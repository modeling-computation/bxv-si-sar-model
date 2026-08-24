rm(list=ls())

set.seed(1234)
libraries = c("deSolve", "dplyr","reshape2","openxlsx","tidyverse","lubridate",
              "readxl","ggplot2","viridisLite","colorspace","ggpubr","gridExtra",
              "png","grid","magrittr","scales","RColorBrewer","fBasics", "mvtnorm")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE) }

datapath <- paste0(getwd(),"/../Data/")
savepath <- paste0(getwd(), "/../Result/")
figurepath <- paste0(getwd(), "/../Figure/")

df <- read_excel(paste0(datapath, "Data.xlsx"))

df <- df %>%
  mutate(serial = case_when(
    serial == "1day" ~ 1,
    serial == "2days" ~ 2,
    serial == "3days" ~ 3,
    serial == "4days" ~ 4,
    serial == "5days" ~ 5,
    serial == "≥6days" ~ 6,
    TRUE ~ 0  ))

df <- df %>% filter(serial > 0)

df <- df %>%
  mutate(delay = case_when(
    dur_onset_med_index == "≤1day" ~ 1,
    dur_onset_med_index == "≤2days" ~ 2,
    dur_onset_med_index == "≤3days" ~ 3,
    dur_onset_med_index == "≤4days" ~ 4,
    dur_onset_med_index == "≤5days" ~ 5,
    dur_onset_med_index == ">5days" ~ 6,
    TRUE ~ NA_real_  )) # NA: no treatment

wotrt <- df %>% filter(anti=="Not visited")
wotrt_serial_freq <- as.data.frame(table(wotrt$serial), stringsAsFactors = FALSE)
names(wotrt_serial_freq) <- c("serial", "Freq")
wotrt_serial_freq$serial <- as.numeric(wotrt_serial_freq$serial)

wtrt <- df %>% filter(anti=="Baloxavir")

treated_serial <- wtrt$serial
treated_delay <- wtrt$delay

# 1. g(s): dist of SI for index cases without treatment
g_pdf <- function(s, g_par1, g_par2, dist_name) {
  density_values <- switch(
    dist_name,
    "Gamma"   = dgamma(s, shape = g_par1, scale = g_par2),
    "Weibull" = dweibull(s, shape = g_par1, scale = g_par2)
  )
  density_values[s < 0] <- 0
  return(density_values)
}

# g'(s)에서 적분하기 위한 g(s) 이산 확률 (P(t-1 < X <= t))
# day=1~5: P(day-1 < X <= day) = F(day) - F(day-1)
# day>=6 (right-censored, "≥6days"): P(X > 5) = 1 - F(5)
prob_g_discrete <- function(t, g_par1, g_par2, dist_name) {
  sapply(t, function(day) {
    if (day >= 6) {
      prob_discrete <- switch(
        dist_name,
        "Gamma"   = 1 - pgamma(5, shape = g_par1, scale = g_par2),
        "Weibull" = 1 - pweibull(5, shape = g_par1, scale = g_par2)
      )
    } else {
      prob_discrete <- switch(
        dist_name,
        "Gamma" = pgamma(day, shape = g_par1, scale = g_par2) - pgamma(day - 1, shape = g_par1, scale = g_par2),
        "Weibull" = pweibull(day, shape = g_par1, scale = g_par2) - pweibull(day - 1, shape = g_par1, scale = g_par2)
      )
    }
      return(prob_discrete)
  })
}


# 2. g'(s): dist of SI for index cases with treatment (Baloxavir)
g_prime_pdf <- function(s, st, epsilon, g_par1, g_par2, dist_name) {
  g_pdf_numerator <- function(x) g_pdf(x, g_par1, g_par2, dist_name)
  g_pdf_reduced_numerator <- function(x) (1 - epsilon) * g_pdf(x, g_par1, g_par2, dist_name)
  
  integral_part1 <- integrate(g_pdf_numerator, lower = 0, upper = st)$value
  integral_part2 <- integrate(g_pdf_reduced_numerator, lower = st, upper = Inf)$value
  denominator <- integral_part1 + integral_part2
  
  if (denominator < 1e-9) return(rep(0, length(s)))
  
  numerator <- ifelse(s < st, g_pdf(s, g_par1, g_par2, dist_name), (1 - epsilon) * g_pdf(s, g_par1, g_par2, dist_name))
  return(numerator / denominator)
}

prob_g_prime_discrete <- function(s_day, st, epsilon, g_par1, g_par2, dist_name = dist) {
  integrand_func <- function(s) g_prime_pdf(s, st, epsilon, g_par1, g_par2, dist_name)
  sapply(s_day, function(day) {
    if (day >= 6) {
      cdf_5 <- integrate(integrand_func, lower = 0, upper = 5)$value
      return(1 - cdf_5)
    } else {
      return(integrate(integrand_func, lower = max(0, day - 1), upper = day)$value)
    }
  })
}

calculate_p_i_part <- function(st, epsilon, g_par1, g_par2, dist_name) {
  g_pdf_numerator <- function(x) g_pdf(x, g_par1, g_par2, dist_name)
  g_pdf_reduced_numerator <- function(x) (1 - epsilon) * g_pdf(x, g_par1, g_par2, dist_name)
  
  integral_part1 <- integrate(g_pdf_numerator, lower = 0, upper = st)$value
  integral_part2 <- integrate(g_pdf_reduced_numerator, lower = st, upper = Inf)$value
  return(integral_part1 + integral_part2)
}

Combined_MLE <- function(params, data_list, dist_name = dist) {
  g_par1  <- params[1]
  g_par2  <- params[2]
  epsilon <- params[3]
  k       <- params[4] 
  
  wotrt_data <- data_list$wotrt_data
  wtrt_data  <- data_list$wtrt_data
  wotrt_data_full <- data_list$wotrt_data_full
  
  probs_L1 <- prob_g_discrete(wotrt_data$serial, g_par1, g_par2, dist_name)
  probs_L1[probs_L1 <= 0] <- 1e-15
  neg_log_L1 <- -sum(wotrt_data$Freq * log(probs_L1))
  
  log_probs_L3 <- sapply(1:nrow(wtrt_data), function(i) {
    prob <- prob_g_prime_discrete(s_day = wtrt_data$serial[i], st = wtrt_data$delay[i], epsilon = epsilon, 
                                  g_par1 = g_par1, g_par2 = g_par2, dist_name = dist_name)
    if (prob <= 0) prob <- 1e-15
    return(log(prob))
  })
  neg_log_L3 <- -sum(log_probs_L3)
  
  log_probs_L2_1 <- sapply(1:nrow(wotrt_data_full), function(i) {
    p_i <- max(1e-9, min(k, 1 - 1e-9)) 
    dbinom(wotrt_data_full$infection_family_count[i], wotrt_data_full$n_member[i], prob = p_i, log = TRUE)
  })
  
  log_probs_L2_2 <- sapply(1:nrow(wtrt_data), function(i) {
     potential_i <- calculate_p_i_part(st = wtrt_data$delay[i], 
                                        epsilon = epsilon, 
                                        g_par1 = g_par1, g_par2 = g_par2, 
                                        dist_name = dist_name)
      p_i <- k * potential_i
      p_i <- max(1e-9, min(p_i, 1 - 1e-9))
      dbinom(wtrt_data$infection_family_count[i], wtrt_data$n_member[i], prob = p_i, log = TRUE)
  })
  
  neg_log_L2 <- -sum(log_probs_L2_1) - sum(log_probs_L2_2)

  # total negative log likelihood
  total_neg_log_L <- neg_log_L1 + neg_log_L2 + neg_log_L3
  
  print(paste("Params:", paste(round(params, 3), collapse=", "), "| -LogL:", round(total_neg_log_L, 2)))
  
  return(total_neg_log_L)
}

gamma_mean_sd <- function(shape, scale) {
  mean_val <- shape * scale
  sd_val <- sqrt(shape) * scale
  return(c(mean = mean_val, sd = sd_val))
}

weibull_mean_sd <- function(shape, scale) {
  mean_val <- scale * gamma(1 + 1 / shape)
  sd_val <- scale * sqrt(gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2)
  return(c(mean = mean_val, sd = sd_val))
}

data_for_mle <- list(
  wotrt_data = wotrt_serial_freq,
  wtrt_data = wtrt,
  wotrt_data_full = wotrt
)


initial_params <- c(2.0, 3.0, 0.5, 0.1) 
lower_bounds <- c(0.01, 0.01, 0, 0)
upper_bounds <- c(Inf, Inf, 1, 1)

dist_candidate <- c("Gamma", "Weibull")
results_list <- list()

for (dist in dist_candidate) {
  cat(paste("\nFitting model with", dist, "distribution. \n"))
  optim_result <- optim( par = initial_params, fn = Combined_MLE, method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds,
                         data_list = data_for_mle, dist_name = dist, 
                         hessian=TRUE)
  
  fisher_info <- solve(optim_result$hessian)
  n_samples <- 10000
  sampled_params <- rmvnorm(n = n_samples, mean = optim_result$par, sigma = fisher_info)
  colnames(sampled_params) <- c("g_par1", "g_par2", "epsilon", "k")
  df_samples <- as.data.frame(sampled_params)
  
  df_samples$g_par1  <- pmax(df_samples$g_par1, 0.01)
  df_samples$g_par2  <- pmax(df_samples$g_par2, 0.01) 
  df_samples$epsilon <- pmax(pmin(df_samples$epsilon, 1), 0)
  df_samples$k       <- pmax(df_samples$k, 0)
  ci_results <- apply(df_samples, 2, function(x) quantile(x, probs = c(0.25, 0.75)))
  final_ci <- data.frame(
    Parameter = colnames(sampled_params),
    Estimate = optim_result$par,
    Lower_CI = ci_results[1, ],
    Upper_CI = ci_results[2, ]
  )
  print("--- Simulation-based 80% Confidence Intervals ---")
  print(final_ci)
  lower_ci <- final_ci$Lower_CI
  upper_ci <- final_ci$Upper_CI
  
  
  if (optim_result$convergence == 0) {
    aic_value <- 2 * length(initial_params) + 2 * optim_result$value
    mean <- if (dist == "Gamma") {
      gamma_mean_sd(optim_result$par[1], optim_result$par[2])[1]
    } else {
      weibull_mean_sd(optim_result$par[1], optim_result$par[2])[1]
    }
    mean_lower <- if (dist == "Gamma") {
      gamma_mean_sd(lower_ci[1], lower_ci[2])[1]
    } else {
      weibull_mean_sd(lower_ci[1], lower_ci[2])[1]
    }
    mean_upper <- if (dist == "Gamma") {
      gamma_mean_sd(upper_ci[1], upper_ci[2])[1]
    } else {
      weibull_mean_sd(upper_ci[1], upper_ci[2])[1]
    }
    sd <- if (dist == "Gamma") {
      gamma_mean_sd(optim_result$par[1], optim_result$par[2])[2]
    } else {
      weibull_mean_sd(optim_result$par[1], optim_result$par[2])[2]
    }
    sd_lower <- if (dist == "Gamma") {
      gamma_mean_sd(lower_ci[1], lower_ci[2])[2]
    } else {
      weibull_mean_sd(lower_ci[1], lower_ci[2])[2]
    }
    sd_upper <- if (dist == "Gamma") {
      gamma_mean_sd(upper_ci[1], upper_ci[2])[2]
    } else {
      weibull_mean_sd(upper_ci[1], upper_ci[2])[2]
    }
    
    results_list[[dist]] <- list(
      distribution = dist,
      estimated_params = optim_result$par,
      # se = se,
      CI_lower = lower_ci,
      CI_upper = upper_ci,
      g_mean = mean,
      g_mean_lower = mean_lower,
      g_mean_upper = mean_upper,
      g_sd = sd,
      g_sd_lower = sd_lower,
      g_sd_upper = sd_upper,
      neg_log_L = optim_result$value,
      aic = aic_value,
      convergence = TRUE
    )
    cat(paste(dist, "model fitting successful. AIC:", round(aic_value, 2), "\n"))
  } else { 
    results_list[[dist]] <- list(
      distribution = dist,
      convergence = FALSE,
      message = optim_result$message
    )
    cat(paste(dist, "model fitting failed.\n"))
  }
}


#### Result
successful_fits <- Filter(function(x) x$convergence, results_list)

if (length(successful_fits) > 0) {
  best_model <- successful_fits[[which.min(sapply(successful_fits, `[[`, "aic"))]]
  estimated_params <- best_model$estimated_params
  names(estimated_params) <- c("g_par1", "g_par2", "epsilon", "k")
  
  model_list <- list()
  for (dist in names(successful_fits)) {
    model_info <- successful_fits[[dist]]
    model_list[[dist]] <- data.frame(
      Distribution = model_info$distribution,
      
      g_par1 = model_info$estimated_params[1],
      g_par1_lower = model_info$CI_lower[1],
      g_par1_upper = model_info$CI_upper[1],
      
      g_par2 = model_info$estimated_params[2],
      g_par2_lower = model_info$CI_lower[2],
      g_par2_upper = model_info$CI_upper[2],
      
      g_mean = model_info$g_mean,
      g_mean_lower = model_info$g_mean_lower,
      g_mean_upper = model_info$g_mean_upper,
      
      g_sd = model_info$g_sd,
      g_sd_lower = model_info$g_sd_lower,
      g_sd_upper = model_info$g_sd_upper,
      
      epsilon = model_info$estimated_params[3],
      epsilon_lower = model_info$CI_lower[3],
      epsilon_upper = model_info$CI_upper[3],
      
      k = model_info$estimated_params[4],
      k_lower = model_info$CI_lower[4],
      k_upper = model_info$CI_upper[4],
      
      AIC = model_info$aic
    )
  }
  best_model_df <- do.call(rbind, model_list)
  write.xlsx(best_model_df, file = paste0(savepath, "All_Model_Fitting_Results.xlsx"), rowNames = FALSE)
  
  
  est <- best_model$estimated_params
  lower <- best_model$CI_lower
  upper <- best_model$CI_upper
  
} else {
  print("All model fittings failed. Cannot select a best model.")
}

write.xlsx(list("Best Model" = data.frame(
  Distribution = best_model$distribution,

  g_par1 = best_model$estimated_params[1],
  g_par1_lower = best_model$CI_lower[1],
  g_par1_upper = best_model$CI_upper[1],

  g_par2 = best_model$estimated_params[2],
  g_par2_lower = best_model$CI_lower[2],
  g_par2_upper = best_model$CI_upper[2],

  g_mean = best_model$g_mean,
  g_mean_lower = best_model$g_mean_lower,
  g_mean_upper = best_model$g_mean_upper,

  g_sd = best_model$g_sd,
  g_sd_lower = best_model$g_sd_lower,
  g_sd_upper = best_model$g_sd_upper,

  epsilon = best_model$estimated_params[3],
  epsilon_lower = best_model$CI_lower[3],
  epsilon_upper = best_model$CI_upper[3],

  k = best_model$estimated_params[4],
  k_lower = best_model$CI_lower[4],
  k_upper = best_model$CI_upper[4],

  AIC = best_model$aic
)), file = paste0(savepath, "Best_Model_Estimation.xlsx"), rowNames = FALSE)


best_model <- read.xlsx(paste0(savepath, "Best_Model_Estimation.xlsx"), sheet = "Best Model")
