# Update: 2026.08.05
# revision용) serial interval>=1인 경우만 모델링

rm(list=ls())

set.seed(1234)
libraries = c("deSolve", "dplyr","reshape2","openxlsx","tidyverse","lubridate",
              "readxl","ggplot2","viridisLite","colorspace","ggpubr","gridExtra",
              "png","grid","magrittr","scales","RColorBrewer","fBasics", "mvtnorm")
for(x in libraries) { library(x,character.only=TRUE,warn.conflicts=FALSE) }

datapath <- paste0(getwd(),"/../Data/")
savepath <- paste0(getwd(), "/../Result/Revision_SI_over0_censoring/")
figurepath <- paste0(getwd(), "/../Figure/Revision_SI_over0_censoring/")

df <- read_excel(paste0(datapath, "Data.xlsx"))

# serial, delay (time from onset to treatment) 변수
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

# index cases without treatment
wotrt <- df %>% filter(anti=="Not visited")
wotrt_serial_freq <- as.data.frame(table(wotrt$serial), stringsAsFactors = FALSE)
names(wotrt_serial_freq) <- c("serial", "Freq")
wotrt_serial_freq$serial <- as.numeric(wotrt_serial_freq$serial)

# index cases with treatment (Baloxavir)
wtrt <- df %>% filter(anti=="Baloxavir")

treated_serial <- wtrt$serial
treated_delay <- wtrt$delay


#### Functions

# 1. g(s): dist of SI for index cases without treatment
g_pdf <- function(s, g_par1, g_par2, dist_name) {
  density_values <- switch(
    dist_name,
    "Gamma"   = dgamma(s, shape = g_par1, scale = g_par2),
    "Weibull" = dweibull(s, shape = g_par1, scale = g_par2)
  )
  density_values[s < 0] <- 0 # 음수 값에 대해서는 0으로 설정
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
  
  if (denominator < 1e-9) return(rep(0, length(s))) # 분모가 너무 작으면 0으로
  
  numerator <- ifelse(s < st, g_pdf(s, g_par1, g_par2, dist_name), (1 - epsilon) * g_pdf(s, g_par1, g_par2, dist_name))
  return(numerator / denominator)
}

# g'(s) 이산 확률 (P(t-1 < X' <= t))
# day=1~5: P(day-1 < X' <= day)
# day>=6 (right-censored, "≥6days"): P(X' > 5) = 1 - integral_0^5 g'(s) ds
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



## p_i 관련 함수
# p_i=k[]의 [] 부분 계산
calculate_p_i_part <- function(st, epsilon, g_par1, g_par2, dist_name) {
  g_pdf_numerator <- function(x) g_pdf(x, g_par1, g_par2, dist_name)
  g_pdf_reduced_numerator <- function(x) (1 - epsilon) * g_pdf(x, g_par1, g_par2, dist_name)
  
  integral_part1 <- integrate(g_pdf_numerator, lower = 0, upper = st)$value
  integral_part2 <- integrate(g_pdf_reduced_numerator, lower = st, upper = Inf)$value
  return(integral_part1 + integral_part2)
}


#### loglikelihood (negative llh)

Combined_MLE <- function(params, data_list, dist_name = dist) {
  g_par1  <- params[1] # g(s) para: shape
  g_par2  <- params[2] # g(s) para: scale
  epsilon <- params[3] # 치료 효과
  k       <- params[4] # 비례 상수 (SAR & SI)
  
  wotrt_data <- data_list$wotrt_data
  wtrt_data  <- data_list$wtrt_data
  wotrt_data_full <- data_list$wotrt_data_full
  
  # L1 계산: 치료받지 않은 index case의 Serial Interval
  probs_L1 <- prob_g_discrete(wotrt_data$serial, g_par1, g_par2, dist_name)
  probs_L1[probs_L1 <= 0] <- 1e-15 # log(0) 방지
  neg_log_L1 <- -sum(wotrt_data$Freq * log(probs_L1))
  
  # L3 계산: Baloxavir 치료받은 index case의 Serial Interval
  log_probs_L3 <- sapply(1:nrow(wtrt_data), function(i) {
    prob <- prob_g_prime_discrete(s_day = wtrt_data$serial[i], st = wtrt_data$delay[i], epsilon = epsilon, 
                                  g_par1 = g_par1, g_par2 = g_par2, dist_name = dist_name)
    if (prob <= 0) prob <- 1e-15
    return(log(prob))
  })
  neg_log_L3 <- -sum(log_probs_L3)
  
  
  
  # L2 계산: 치료 받은/안받은 그룹의 secondary infections
  log_probs_L2_1 <- sapply(1:nrow(wotrt_data_full), function(i) {
    p_i <- max(1e-9, min(k, 1 - 1e-9)) # 확률 p는 0과 1 사이
    dbinom(wotrt_data_full$infection_family_count[i], wotrt_data_full$n_member[i], prob = p_i, log = TRUE)
  })
  
  log_probs_L2_2 <- sapply(1:nrow(wtrt_data), function(i) {
     potential_i <- calculate_p_i_part(st = wtrt_data$delay[i], 
                                        epsilon = epsilon, 
                                        g_par1 = g_par1, g_par2 = g_par2, 
                                        dist_name = dist_name)
      p_i <- k * potential_i
      p_i <- max(1e-9, min(p_i, 1 - 1e-9)) # 확률 p는 0과 1 사이
      dbinom(wtrt_data$infection_family_count[i], wtrt_data$n_member[i], prob = p_i, log = TRUE)
  })
  
  neg_log_L2 <- -sum(log_probs_L2_1) - sum(log_probs_L2_2)

  # total negative log likelihood
  total_neg_log_L <- neg_log_L1 + neg_log_L2 + neg_log_L3
  
  print(paste("Params:", paste(round(params, 3), collapse=", "), "| -LogL:", round(total_neg_log_L, 2)))
  
  return(total_neg_log_L)
}



#### 최적화 실행

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
  
  fisher_info <- solve(optim_result$hessian) # covariance matrix
  n_samples <- 10000 # 샘플 10000개 뽑음
  sampled_params <- rmvnorm(n = n_samples, mean = optim_result$par, sigma = fisher_info)
  colnames(sampled_params) <- c("g_par1", "g_par2", "epsilon", "k")
  df_samples <- as.data.frame(sampled_params)
  
  df_samples$g_par1  <- pmax(df_samples$g_par1, 0.01) # shape > 0
  df_samples$g_par2  <- pmax(df_samples$g_par2, 0.01) # scale > 0
  df_samples$epsilon <- pmax(pmin(df_samples$epsilon, 1), 0) # 0 <= epsilon <= 1
  df_samples$k       <- pmax(df_samples$k, 0)       # 0 <= k <= 1
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
  # AIC가 가장 낮은 모델을 선택
  best_model <- successful_fits[[which.min(sapply(successful_fits, `[[`, "aic"))]]
  estimated_params <- best_model$estimated_params
  names(estimated_params) <- c("g_par1", "g_par2", "epsilon", "k")
  
  # 모든 분포 후보에 대한 피팅 결과 저장
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
  write.xlsx(best_model_df, file = paste0(savepath, "All_Model_Fitting_Results_IQR.xlsx"), rowNames = FALSE)
  
  
  est <- best_model$estimated_params
  lower <- best_model$CI_lower
  upper <- best_model$CI_upper
  
  
  print("--------------------------------------------------")
  print("           Best Model Selection Results           ")
  cat(sprintf("최적 분포 (Best Distribution): %s (AIC: %.2f)\n", best_model$distribution, best_model$aic))
  print("Estimated parameters with 95% CI from best model:")
  cat(sprintf("Shape (g_par1): %.4f (%.4f - %.4f)\n", est[1], lower[1], upper[1]))
  cat(sprintf("Scale (g_par2): %.4f (%.4f - %.4f)\n", est[2], lower[2], upper[2]))
  cat(sprintf("g mean        : %.4f (%.4f - %.4f)\n", best_model$g_mean, best_model$g_mean_lower, best_model$g_mean_upper))
  cat(sprintf("g sd          : %.4f (%.4f - %.4f)\n", best_model$g_sd, best_model$g_sd_lower, best_model$g_sd_upper))
  cat(sprintf("Epsilon       : %.4f (%.4f - %.4f)\n", est[3], lower[3], upper[3]))
  cat(sprintf("k             : %.4f (%.4f - %.4f)\n", est[4], lower[4], upper[4]))
  print("--------------------------------------------------")
  
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
)), file = paste0(savepath, "Best_Model_Estimation_IQR.xlsx"), rowNames = FALSE)


best_model <- read.xlsx(paste0(savepath, "Best_Model_Estimation_IQR.xlsx"), sheet = "Best Model")




######################################################
## Plot
######################################################

# 최적 모델의 g(s)와 g'(s) ggplot 시각화
g_par1 <- best_model$g_par1
g_par2 <- best_model$g_par2
mean <- best_model$g_mean
# sd <- best_model$g_sd
st <- as.numeric(names(sort(table(treated_delay), decreasing=TRUE)[1]))
epsilon <- best_model$epsilon
s_values <- seq(0, 10, length.out = 100)
dist_name <- best_model$Distribution

g_values <- g_pdf(s_values, g_par1, g_par2, dist_name)
g_prime_values <- g_prime_pdf(s_values, st, epsilon, g_par1, g_par2, dist_name)
df_g <- data.frame(s = s_values, density = g_values)
df_g_prime <- data.frame(s = s_values, density = g_prime_values)
df_g$Type <- "g(s)" ; g_prime_name <- "g'(s)"
df_g_prime$Type <- g_prime_name
df_all <- rbind(df_g, df_g_prime)


# max_density <- max(df_all$density, na.rm = TRUE) # 2y축
# max_freq <- max(wotrt_serial_freq$Freq, na.rm = TRUE) # 2y축
# scaling_factor <- max_density / max_freq # 2y축
# ggplot(df_all, aes(x = s, y = density, color = Type)) + geom_line(linewidth = 1) +
#   geom_point(data = wotrt_serial_freq, aes(x = as.numeric(serial), y = Freq * scaling_factor, color = "Population"), size = 2.5, inherit.aes = FALSE) +
#   geom_bar(data = wotrt_serial_freq, aes(x = as.numeric(serial), y = Freq * scaling_factor), stat = "identity", fill = "grey", color = "darkgrey", alpha = 0.3, inherit.aes = FALSE) +
#   labs(title = paste0("g(s) and g'(s)"), x = "Serial interval (days)", color = "") +
#   # geom_vline(xintercept = st, linetype = "dotted", color = "black", linewidth=1) +
#   scale_color_manual( name = "", breaks = c("g(s)", g_prime_name, "Population"),
#                       values = setNames(c("blue", "red", "darkgreen", "black"), c("g(s)", g_prime_name, "Population"))) +
#   scale_y_continuous(name = "Density", sec.axis = sec_axis(trans = ~ . / scaling_factor, name = "Population")) +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), 
#         axis.title.y.right = element_text(color = "black"), axis.text.y.right = element_text(color = "black"))
# 
# # ggsave(filename = paste0(figurepath, "After normalize.png"),
# #        width = 10, height = 6, units = "in", dpi = 300)


## Before normalize 시각화

# g'(s)
g_prime_pdf_2 <- function(s, st, epsilon, g_par1, g_par2, dist_name) {
  g_pdf_numerator <- function(x) g_pdf(x, g_par1, g_par2, dist_name)
  g_pdf_reduced_numerator <- function(x) (1 - epsilon) * g_pdf(x, g_par1, g_par2, dist_name)
  
  integral_part1 <- integrate(g_pdf_numerator, lower = 0, upper = st)$value
  integral_part2 <- integrate(g_pdf_reduced_numerator, lower = st, upper = Inf)$value
  
  numerator <- ifelse(s < st, g_pdf(s, g_par1, g_par2, dist_name), (1 - epsilon) * g_pdf(s, g_par1, g_par2, dist_name))
  return(numerator)
}

g_values <- g_pdf(s_values, g_par1, g_par2, dist_name)
g_prime_values_norm <- g_prime_pdf_2(s_values, st, epsilon, g_par1, g_par2, dist_name)
df_g <- data.frame(s = s_values, density = g_values)
df_g_prime_2 <- data.frame(s = s_values, density = g_prime_values_norm)

df_g$Type <- "g(s)" ; df_g_prime_2$Type <- "g'(s)"
df_all_2 <- rbind(df_g_prime_2, df_g)

# max_density <- max(df_all_2$density, na.rm = TRUE) # 2y축
# max_freq <- max(wotrt_serial_freq$Freq, na.rm = TRUE) # 2y축
# scaling_factor <- max_density / max_freq # 2y축
# ggplot(df_all_2, aes(x = s, y = density, color = Type)) + geom_line(linewidth = 1) +
#   geom_point(data = wotrt_serial_freq, aes(x = as.numeric(serial), y = Freq * scaling_factor, color = "Population"), size = 2.5, inherit.aes = FALSE) +
#   geom_bar(data = wotrt_serial_freq, aes(x = as.numeric(serial), y = Freq * scaling_factor), stat = "identity", fill = "grey", color = "darkgrey", alpha = 0.3, inherit.aes = FALSE) +
#   labs(title = paste0("g(s) and g'(s)"), x = "Serial interval (days)", color = "") +
#   # geom_vline(xintercept = trtdelay_1, linetype = "dotted", color = "black", linewidth=1) +
#   scale_color_manual( name = "", breaks = c("g(s)", g_prime_name, "Population"),
#                       values = setNames(c("blue", "red", "black"), c("g(s)", g_prime_name, "Population"))) +
#   scale_y_continuous(name = "Density", sec.axis = sec_axis(trans = ~ . / scaling_factor, name = "Population")) +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 12), 
#         axis.title.y.right = element_text(color = "black"), axis.text.y.right = element_text(color = "black"))
# 
# # ggsave(filename = paste0(figurepath, "Before normalize.png"),
# # width = 10, height = 6, units = "in", dpi = 300)




######################################################
## SAR Comparison Data Generation (For Python Plotting)
######################################################

#  Observed SAR (Actual Data Aggregation)
# 치료 지연일(delay) 별로 그룹화하여 실제 2차 감염률 계산
sar_obs_df <- wtrt %>%
  group_by(delay) %>%
  summarise(
    total_contacts = sum(n_member),          # 분모: 전체 가족 구성원 수 (Likelihood 식의 n_i)
    total_infections = sum(infection_family_count), # 분자: 실제 2차 감염자 수 (Likelihood 식의 m_i)
    obs_sar = total_infections / total_contacts,    # 관측된 SAR
    count = n() # 해당 delay에 속하는 index case 수
  ) %>%
  ungroup()

# 관측된 SAR의 95% 신뢰구간 계산 (Binomial Confidence Interval using Wilson Score interval)
# Python에서 Error bar를 그리기 위함
sar_ci <- t(sapply(1:nrow(sar_obs_df), function(i) {
  test_result <- prop.test(x = sar_obs_df$total_infections[i],
                           n = sar_obs_df$total_contacts[i],
                           conf.level = 0.95)
  return(test_result$conf.int)
}))

sar_obs_df$sar_lower <- sar_ci[, 1]
sar_obs_df$sar_upper <- sar_ci[, 2]

# 2. Estimated SAR Curve (Model Prediction)
hat_g_par1  <- best_model$g_par1
hat_g_par2  <- best_model$g_par2
hat_epsilon <- best_model$epsilon
hat_k       <- best_model$k
hat_dist    <- best_model$Distribution

# Delay를 0일부터 6일(혹은 데이터 최대값)까지 촘촘하게 생성
max_delay <- max(wtrt$delay, na.rm=TRUE)
delay_seq <- seq(0, max_delay + 1, by = 0.1)

# 각 delay 시점에 대한 p_i (Estimated SAR) 계산
# p_i = k * [Integral of un-normalized modified density]
est_sar_values <- sapply(delay_seq, function(t) {
  integral_val <- calculate_p_i_part(st = t, epsilon = hat_epsilon, g_par1 = hat_g_par1, g_par2 = hat_g_par2, dist_name = hat_dist)
  p_val <- hat_k * integral_val
  return(p_val)
})

sar_est_df <- data.frame(
  delay = delay_seq,
  est_sar = est_sar_values
)

output_filename <- paste0(savepath, "SAR_Comparison.xlsx")
wb <- createWorkbook()
# Sheet 1: Observed Data (Points & Error Bars)
addWorksheet(wb, "Observed_SAR")
writeData(wb, "Observed_SAR", sar_obs_df)

# Sheet 2: Estimated Curve (Line)
addWorksheet(wb, "Estimated_SAR_Curve")
writeData(wb, "Estimated_SAR_Curve", sar_est_df)
# Sheet 3: Metadata (Parameters used)
meta_df <- data.frame(
  Parameter = c("Distribution", "Shape", "Scale", "Epsilon", "k"),
  Value = c(hat_dist, hat_g_par1, hat_g_par2, hat_epsilon, hat_k)
)
addWorksheet(wb, "Model_Metadata")
writeData(wb, "Model_Metadata", meta_df)
saveWorkbook(wb, output_filename, overwrite = TRUE)




######################################################
## Mean
######################################################

# # wotrt에 대한 serial interval 값별 n수, 0인 경우를 제외했을 때의 평균과 표준편차 구하기
# table(wotrt$serial)
# wotrt_serial_over0 <- wotrt%>%filter(serial!=6)
# mean(wotrt_serial_over0$serial)
# sd(wotrt_serial_over0$serial)
# 
# # within 24 hours
# wtrt_24 <- wtrt %>% filter(delay == 1)
# table(wtrt_24$serial)
# wtrt_24_serial_over0 <- wtrt_24 %>% filter(serial!=6)
# mean(wtrt_24_serial_over0$serial)
# sd(wtrt_24_serial_over0$serial)
# 
# # within 48 hours
# wtrt_48 <- wtrt %>% filter(delay == 2)
# table(wtrt_48$serial)
# wtrt_48_serial_over0 <- wtrt_48 %>% filter(serial!=6)
# mean(wtrt_48_serial_over0$serial)
# sd(wtrt_48_serial_over0$serial)
# 
# # within 72 hours
# wtrt_72 <- wtrt %>% filter(delay == 3)
# table(wtrt_72$serial)
# wtrt_72_serial_over0 <- wtrt_72 %>% filter(serial!=6)
# mean(wtrt_72_serial_over0$serial)
# sd(wtrt_72_serial_over0$serial)
# 
# # within 96 hours
# wtrt_96 <- wtrt %>% filter(delay == 4)
# table(wtrt_96$serial)
# wtrt_96_serial_over0 <- wtrt_96 %>% filter(serial!=6)
# mean(wtrt_96_serial_over0$serial)
# sd(wtrt_96_serial_over0$serial)
# 
# # within 120 hours
# wtrt_120 <- wtrt %>% filter(delay == 5)
# table(wtrt_120$serial)
# wtrt_120_serial_over0 <- wtrt_120 %>% filter(serial!=6)
# mean(wtrt_120_serial_over0$serial)
# sd(wtrt_120_serial_over0$serial)
# 
# # over 120 hours
# wtrt_over120 <- wtrt %>% filter(delay == 6)
# table(wtrt_over120$serial)
# wtrt_over120_serial_over0 <- wtrt_over120 %>% filter(serial!=6)
# mean(wtrt_over120_serial_over0$serial)
# sd(wtrt_over120_serial_over0$serial)
# 
# 
# # 전체에 대한 serial interval 평균 (no transmission, 6이상 제외)
# df_summary <- df %>% filter(anti=="Baloxavir" | anti=="Not visited")
# nrow(df_summary)
# table(df_summary$serial, df_summary$delay)
# df_summary_serial_over0 <- df_summary %>% filter(serial!=6)
# mean(df_summary_serial_over0$serial)
# sd(df_summary_serial_over0$serial)


# ==============================================================================
# 6. Model Validation: Observed SAR vs Predicted SAR (Individual Household Level)
# ==============================================================================

# 1. k 파라미터 가져오기 및 개별 가구 SAR 계산
k_est <- best_model$k 

# 데이터프레임에 개별 가구의 SAR (ind_sar) 계산 추가
# infection_family_count: 2차 감염자 수, n_member: 인덱스 케이스 제외한 가족 수
df <- df %>%
  mutate(ind_sar = infection_family_count / n_member)

# 2. 그룹별 관찰된 "평균 SAR" (Mean of Individual SARs) 및 95% 신뢰구간 계산
# 비치료군 (No treatment)
df_wotrt <- df %>% filter(anti == "Not visited")
n_wotrt <- nrow(df_wotrt)
mean_sar_wotrt <- mean(df_wotrt$ind_sar, na.rm = TRUE)
se_sar_wotrt <- sd(df_wotrt$ind_sar, na.rm = TRUE) / sqrt(n_wotrt) # 표준 오차(SE)

# 치료군 (Delay 1~5일)
taus <- 1:5
sar_wtrt_list <- lapply(taus, function(t) {
  df_sub <- df %>% filter(anti == "Baloxavir" & delay == t)
  n_sub <- nrow(df_sub)
  
  if(n_sub == 0) {
    return(data.frame(delay = t, n = 0, obs_sar = NA, lower = NA, upper = NA))
  }
  
  mean_sar <- mean(df_sub$ind_sar, na.rm = TRUE)
  se_sar <- sd(df_sub$ind_sar, na.rm = TRUE) / sqrt(n_sub)
  
  # 95% 신뢰구간 (정규분포 근사: Mean ± 1.96 * SE), 확률이므로 0~1 사이로 제한
  lower_ci <- max(0, mean_sar - 1.96 * se_sar)
  upper_ci <- min(1, mean_sar + 1.96 * se_sar)
  
  data.frame(delay = t, n = n_sub, obs_sar = mean_sar, lower = lower_ci, upper = upper_ci)
})
sar_wtrt_df <- do.call(rbind, sar_wtrt_list)

# 관찰 SAR 데이터프레임 병합
obs_sar_df <- rbind(
  data.frame(Group = "No Tx", delay = 0, n = n_wotrt, 
             obs_sar = mean_sar_wotrt, 
             lower = max(0, mean_sar_wotrt - 1.96 * se_sar_wotrt), 
             upper = min(1, mean_sar_wotrt + 1.96 * se_sar_wotrt)),
  data.frame(Group = paste0("Delay ", sar_wtrt_df$delay), delay = sar_wtrt_df$delay, 
             n = sar_wtrt_df$n, 
             obs_sar = sar_wtrt_df$obs_sar, lower = sar_wtrt_df$lower, upper = sar_wtrt_df$upper)
)

# 3. 그룹별 예측된 SAR (Predicted SAR) 계산
# Pred SAR = k * AUC = k * (1 - AUC_reduction)
compute_reduction_auc_cdf <- function(st, epsilon, g_par1, g_par2, dist_name) {
  tail_prob <- switch(dist_name,
                      "Gamma"   = pgamma(st, shape = g_par1, scale = g_par2, lower.tail = FALSE),
                      "Weibull" = pweibull(st, shape = g_par1, scale = g_par2, lower.tail = FALSE))
  return(epsilon * tail_prob)
}

pred_sar_df <- data.frame(
  delay = c(0, taus), # 0은 No Tx
  AUC_reduction = c(0, sapply(taus, function(st) compute_reduction_auc_cdf(st, epsilon, g_par1, g_par2, dist_name)))
) %>%
  mutate(
    AUC_remaining = 1 - AUC_reduction,
    pred_sar = k_est * AUC_remaining
  )

# 관찰값과 예측값 결합
val_df <- left_join(obs_sar_df, pred_sar_df, by = "delay")
val_df$Group <- factor(val_df$Group, levels = c("No Tx", paste0("Delay ", taus)))

# ==============================================================================
# Figure 1 (Final Revision): Observed Mean SAR (Point/Line) vs Predicted SAR
# ==============================================================================

# 1. No treatment 제외 (Delay 1~5만 남기기)
val_df_sub <- val_df %>%
  filter(delay %in% taus)

# 2. Figure 1 생성
fig1 <- ggplot(val_df_sub, aes(x = Group)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = "Observed SAR"), 
                width = 0.15, linewidth = 1) +
  geom_crossbar(aes(y = obs_sar, ymin = obs_sar, ymax = obs_sar, color = "Observed SAR"), 
                width = 0.3, linewidth = 1.2) +
  
  geom_line(aes(y = pred_sar, color = "Estimated SAR", group = 1), 
            linewidth = 1, linetype = "dashed") + 
  geom_point(aes(y = pred_sar, color = "Predicted SAR"), 
             size = 4, shape = 16) +
  
  scale_color_manual(name = "", 
                     values = c("Observed SAR" = "black", 
                                "Predicted SAR" = "darkred")) +
  
  # y축 설정 (0부터 시작하도록 제한)
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  
  labs(x = "Time from symptom onset to treatment",
       y = "Mean SAR with 95% CI") +
  theme_bw() +
  theme(text = element_text(size = 14), legend.position = "top")

print(fig1)
