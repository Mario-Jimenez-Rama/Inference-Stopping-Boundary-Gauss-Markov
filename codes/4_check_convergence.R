library(goffda)
library(ggplot2)
library(plotly)
library(sde)

simulate_gm_process <- function(n_paths, n_steps, delta, x0 = 1, slope_fun_t, trend_fun_t, vol_fun_t, alpha_gm, beta_gm, sigma_gm, subdiv) {
  T <- delta * n_steps
  time_grid <- seq(0, T, length.out = n_steps + 1)
  paths <- matrix(NA, nrow = n_steps + 1, ncol = n_paths)
  paths[1, ] <- x0
  
  for (j in 1:n_paths) {
    for (i in 2:(n_steps + 1)) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x_ant <- paths[i - 1, j]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t, subdiv)
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1, subdiv)
      vol_value <- vol_int_dt(vol_fun_t, alpha_int_fun, t0, t1, subdiv)
      
      mean_gm <- gm_mean_dt(x_ant, t1, alpha_int_fun, trend_value)
      var_gm <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      if (is.na(mean_gm) || is.infinite(mean_gm) || is.na(var_gm) || var_gm < 0) {
        paths[i, j] <- paths[i - 1, j] # If there's an error, keep the previous value
        next
      }
      
      paths[i, j] <- rnorm(1, mean = mean_gm, sd = sqrt(var_gm))
    }
  }
  return(list(time = time_grid, paths = paths))
}

check_convergence_with_estimated_fim <- function(theta_true, theta_test, N, delta, M, alpha_level = 0.05, 
                                                 slope_fun, trend_fun, vol_fun, 
                                                 dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                                 dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, 
                                                 dvol_dtheta1, dvol_dtheta2, dvol_dtheta3, subdivision = 70) {
  
  T <- N * delta
  
  # --- STEP 1: Simulate M paths using the true parameters ---
  cat("STEP 1: Simulating", M, "paths with the true parameters.\n")
  
  params_true <- list(alpha = theta_true[1], beta = theta_true[2], sigma = theta_true[3])
  slope_fun_true <- function(t) do.call(slope_fun, c(list(t = t), params_true))
  trend_fun_true <- function(t) do.call(trend_fun, c(list(t = t), params_true))
  vol_fun_true   <- function(t) do.call(vol_fun, c(list(t = t), params_true))
  
  # Generate all paths
  all_paths <- simulate_gm_process(n_paths = M, n_steps = N, delta = delta, x0 = 1,
                                   slope_fun_t = slope_fun_true, trend_fun_t = trend_fun_true, vol_fun_t = vol_fun_true,
                                   alpha_gm = theta_true[1], beta_gm = theta_true[2], sigma_gm = theta_true[3], subdiv = subdivision)$paths
  
  # --- STEP 2: Estimate parameters and calculate a FIM for EACH path ---
  cat("STEP 2: Performing", M, "estimations and calculating their individual FIMs.\n")
  
  all_estimates <- list()
  all_fims <- list()
  
  for (i in 1:M) {
    if (i %% 20 == 0) cat("  Processing path:", i, "/", M, "\n")
    
    X <- all_paths[, i, drop = FALSE]
    
    # 2a: Get the maximum likelihood estimate for this path
    result <- tryCatch({
      lk_estimates(subdivision, X, T, theta_test = theta_test, 
                   slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
    }, error = function(e) {
      message("Error in estimation (iteration ", i, "): ", e$message)
      return(NULL)
    })
    
    if (is.null(result)) {
      # If estimation fails, we store NULL to skip it later
      all_estimates[[i]] <- NULL
      all_fims[[i]] <- NULL
      next
    }
    
    theta_estimated <- as.numeric(result@coef)
    all_estimates[[i]] <- theta_estimated
    
    # 2b: Calculate the FIM using the newly estimated parameters for this path
    params_est <- list(alpha = theta_estimated[1], beta = theta_estimated[2], sigma = theta_estimated[3])
    
    # Create functions evaluated at the estimated parameters
    slope_fun_est <- function(t) do.call(slope_fun, c(list(t = t), params_est))
    trend_fun_est <- function(t) do.call(trend_fun, c(list(t = t), params_est))
    vol_fun_est   <- function(t) do.call(vol_fun, c(list(t = t), params_est))
    
    dslope_dtheta1_est <- function(t) do.call(dslope_dtheta1, c(list(t = t), params_est))
    dslope_dtheta2_est <- function(t) do.call(dslope_dtheta2, c(list(t = t), params_est))
    dslope_dtheta3_est <- function(t) do.call(dslope_dtheta3, c(list(t = t), params_est))
    
    dtrend_dtheta1_est <- function(t) do.call(dtrend_dtheta1, c(list(t = t), params_est))
    dtrend_dtheta2_est <- function(t) do.call(dtrend_dtheta2, c(list(t = t), params_est))
    dtrend_dtheta3_est <- function(t) do.call(dtrend_dtheta3, c(list(t = t), params_est))
    
    dvol_dtheta1_est <- function(t) do.call(dvol_dtheta1, c(list(t = t), params_est))
    dvol_dtheta2_est <- function(t) do.call(dvol_dtheta2, c(list(t = t), params_est))
    dvol_dtheta3_est <- function(t) do.call(dvol_dtheta3, c(list(t = t), params_est))
    
    # Calculate the FIM for this path and this estimate
    fim_i <- Fisher_matrix_der_exact_complete(
      T, X, subdivision,
      slope_fun = slope_fun_est, trend_fun = trend_fun_est, vol_fun = vol_fun_est,
      dslope_dtheta1 = dslope_dtheta1_est, dslope_dtheta2 = dslope_dtheta2_est, dslope_dtheta3 = dslope_dtheta3_est,
      dtrend_dtheta1 = dtrend_dtheta1_est, dtrend_dtheta2 = dtrend_dtheta2_est, dtrend_dtheta3 = dtrend_dtheta3_est,
      dvol_dtheta1 = dvol_dtheta1_est, dvol_dtheta2 = dvol_dtheta2_est, dvol_dtheta3 = dvol_dtheta3_est
    )$Fisher_matrix
    
    all_fims[[i]] <- fim_i
  }
  
  # --- STEP 3: Average the Fisher matrices ---
  cat("STEP 3: Averaging the calculated Fisher matrices.\n")
  
  # Filter out any FIM that resulted in NULL due to an estimation error
  valid_fims <- all_fims[!sapply(all_fims, is.null)]
  
  if(length(valid_fims) == 0) {
    stop("Could not calculate any Fisher matrix. All estimations failed.")
  }
  
  # Average the Fisher matrices
  fisher_info_averaged <- Reduce('+', valid_fims) / length(valid_fims)
  
  # Calculate the Cholesky decomposition of the AVERAGED matrix
  R_averaged <- tryCatch({
    chol(fisher_info_averaged)
  }, error = function(e) {
    stop("Error in Cholesky with the averaged FIM. The matrix might not be positive definite. Error: ", e$message)
  })
  
  # --- STEP 4: Calculate convergence statistics using the averaged FIM ---
  cat("STEP 4: Calculating convergence statistics.\n")
  
  chi_squared_stats <- numeric(M)
  vectors <- matrix(NA, nrow = M, ncol = length(theta_true))
  
  valid_estimates_count <- 0
  
  for (i in 1:M){
    theta_estimated <- all_estimates[[i]]
    
    if (is.null(theta_estimated)) {
      # If the estimate was null, the statistic is NA
      chi_squared_stats[i] <- NA
      next
    }
    
    # Calculate the Z-vector using the pre-calculated R_averaged
    z <- sqrt(N) * (R_averaged %*% (theta_estimated - theta_true))
    
    chi_squared_stats[i] <- sum(z^2)
    vectors[i,] <- z
    valid_estimates_count <- valid_estimates_count + 1
  }
  
  chi_squared_critical_value <- qchisq(1 - alpha_level, df = length(theta_true))
  coverage_probability <- mean(chi_squared_stats <= chi_squared_critical_value, na.rm = TRUE)
  
  cat("\nAnalysis complete. Obtained", valid_estimates_count, "valid estimates out of", M, ".\n")
  
  return(list(vectors = vectors, 
              chi_squared_stats = chi_squared_stats,
              critical_value = chi_squared_critical_value,
              coverage_probability = coverage_probability,
              estimates = all_estimates,
              averaged_fim = fisher_info_averaged))
}


# Function to generate sphere (no changes)
generate_sphere <- function(radius, n = 100) {
  theta <- seq(0, 2 * pi, length.out = n)
  phi <- seq(0, pi, length.out = n)
  theta_grid <- rep(theta, each = length(phi))
  phi_grid <- rep(phi, length(theta))
  
  x <- radius * sin(phi_grid) * cos(theta_grid)
  y <- radius * sin(phi_grid) * sin(theta_grid)
  z <- radius * cos(phi_grid)
  
  return(list(x = x, y = y, z = z))
}
# OUB 

# alpha_true <- 1 
# beta_true <- 2
# sigma_true <-2
# slope_fun <- function(t, alpha,...)  -alpha/tanh(alpha*(T-t))
# trend_fun <- function(t, alpha, beta, ...) xf*alpha/sinh(alpha*(T-t)) - beta*tanh(alpha*(T-t)/2) 
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true)
# vol_fun_t <- function(t) vol_fun(t, sigma = sigma_true)
# 
# ## derivadas 
# # derivadas de slope:
# dslope_dtheta1 <- function(t, alpha,...) -coth(alpha*(T - t)) + alpha*(T - t)/sinh(alpha*(T - t))^2 
# dslope_dtheta2 <- function(t,...) rep(0, length(t)) 
# dslope_dtheta3 <- function(t,...) rep(0, length(t)) 
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t, alpha, beta,...) xf*csch(alpha*(T - t)) - alpha*xf*(T - t)*csch(alpha*(T - t))*coth(alpha*(T - t)) - beta*(T - t)/2*sech(alpha*(T - t)/2)^2
# dtrend_dtheta2 <- function(t, alpha,...) -tanh(alpha*(T - t)/2)
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t, sigma,...) 2*sigma

# Hull White model periodic mean difficult
# alpha_true <- 3
# beta_true <- 1
# sigma_true <- 2
# slope_fun <- function(t,...) rep(-1, length(t)) # mean reverting
# trend_fun <- function(t, alpha, beta, sigma) alpha + beta*sin(sigma*t)
# vol_fun <- function(t,...) rep(1, length(t)) # unitary volatility
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas 
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(0, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(1, length(t))
# dtrend_dtheta2 <- function(t, sigma,...) sin(sigma*t)
# dtrend_dtheta3 <- function(t, beta, sigma,...) beta*t*cos(sigma*t)
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(0, length(t))

# OU time-dependent volatility complex
# alpha_true <- 1
# beta_true <- 0.5
# sigma_true <- 0.25
# slope_fun <- function(t,...) rep(-1, length(t)) # cte, mean reversion
# #trend_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t, alpha,...) rep(alpha, length(t))
# vol_fun <- function(t, beta, sigma,...) beta*exp(-sigma*t)
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# # derivadas 
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(0, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(1, length(t))
# dtrend_dtheta2 <- function(t,...) rep(0, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,sigma,...) exp(-sigma*t)
# dvol_dtheta3 <- function(t,beta, sigma,...) -beta*t*exp(-sigma*t)

# Hull White model root squared
# alpha_true <- 1 #mean reversion
# beta_true <- 1
# sigma_true <- 0.5
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta, sigma,...) beta + sigma*sqrt(t)
# vol_fun <- function(t,...) rep(1, length(t))
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivatives
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) sqrt(t)
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(0, length(t))

# OU
# alpha_true <- -1 #mean reversion
# beta_true <- 3
# sigma_true <- 2
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t, beta,...) rep(beta, length(t))
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(1, length(t))

# Hull White periodic mean version 1
alpha_true <- -1 # mean reversion
beta_true <- 3
sigma_true <- 0.5
slope_fun <- function(t, alpha,...) rep(alpha, length(t))
trend_fun <- function(t,beta, sigma,...) beta + sin(sigma*t)
vol_fun <- function(t,...) rep(1, length(t))

# functions with true parameters
slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)

## derivadas
# derivadas de slope:
dslope_dtheta1 <- function(t,...) rep(1, length(t))
dslope_dtheta2 <- function(t,...) rep(0, length(t))
dslope_dtheta3 <- function(t,...) rep(0, length(t))

# derivadas de trend
dtrend_dtheta1 <- function(t,...) rep(0, length(t))
dtrend_dtheta2 <- function(t,...) rep(1, length(t))
dtrend_dtheta3 <- function(t, sigma,...) t*cos(sigma*t)

# derivadas de vol
dvol_dtheta1 <- function(t,...) rep(0, length(t))
dvol_dtheta2 <- function(t,...) rep(0, length(t))
dvol_dtheta3 <- function(t,...) rep(0, length(t))

# Hull White periodic mean version 2
# alpha_true <- 3 
# beta_true <- 3
# sigma_true <- 2
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta, sigma,...) beta + sigma*sin(t) 
# vol_fun <- function(t,...) rep(1, length(t)) 
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) 
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t, sigma,...) sin(t)
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(0, length(t))

# time-dependent decaying volatility OU
# alpha_true <- -0.5
# beta_true <- 3
# sigma_true <- 0.5
# slope_fun <- function(t,alpha,...) rep(alpha, length(t)) # mean reversion
# trend_fun <- function(t,beta,...) rep(beta, length(t))
# #vol_fun <- function(t,sigma,...) 5*exp(-sigma*t)*sin(sigma*t) # 1º
# vol_fun <- function(t,sigma,...) 5*exp(-sigma*t) # 2º
# #vol_fun <- function(t, beta, sigma,...) beta + 5*exp(-sigma*t)
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# #dvol_dtheta3 <- function(t, sigma,...) -t*exp(-sigma*t)*sin(sigma*t) + t*exp(-sigma*t)*cos(sigma*t) # 1º
# dvol_dtheta3 <- function(t, sigma,...) -5*t*exp(-sigma*t) # 2º

# Hull White logarithmic mean 
# alpha_true <- -1
# beta_true <- 2
# sigma_true <- 2
# #slope_fun <- function(t,...) rep(-1, length(t)) #mean reversion
# slope_fun <- function(t, alpha,...) rep(alpha, length(t)) #mean reversion # 3º
# #trend_fun <- function(t,alpha, beta,...) alpha + log(1 + beta*t) # 1º
# #trend_fun <- function(t,alpha, beta,...) alpha + beta*log(1 + t) # 2º
# trend_fun <- function(t,beta,...) beta*log(1 + t) # 3º
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas
# # derivadas de slope:
# #dslope_dtheta1 <- function(t,...) rep(0, length(t)) # 1,2º
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) # 3º
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# #dtrend_dtheta1 <- function(t,...) rep(1, length(t)) # 2º
# dtrend_dtheta1 <- function(t,...) rep(0, length(t)) # 3º
# #dtrend_dtheta2 <- function(t, beta,...) t/(1 + beta*t) # 1º
# dtrend_dtheta2 <- function(t,...) log(1 + t)  # 2º
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(1, length(t))

# time-dependent rational function volatility OU
# alpha_true <- -1
# beta_true <- 3
# sigma_true <- 0.5
# slope_fun <- function(t,alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta,...) rep(beta, length(t))
# vol_fun <- function(t,sigma,...) 1/(1 + sigma*t)
# vol_fun <- function(t,sigma,...) t/(1 + sigma*t^2)
# 
# # functions with true parameters
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivadas
# # derivadas de slope:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de trend
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t, sigma,...) -t/(1 + sigma*t)^2


# --- SIMULATION PARAMETERS ---
M <- 200 # Number of paths and estimations
N <- 100 # Number of time steps
delta <- 1 # Time step size
subdivision <- 70 # Subdivision for numerical integration

# Initial values for the optimization
alpha_test <- -0.5
beta_test <- 1
sigma_test <- 1

alpha_level <- 0.05

convergence_results <- check_convergence_with_estimated_fim(
  theta_true = c(alpha_true, beta_true, sigma_true),
  theta_test = c(alpha_test, beta_test, sigma_test),
  N = N, delta = delta, M = M, alpha_level = alpha_level,
  slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun,
  dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3,
  dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
  dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3,
  subdivision = subdivision
)


# --- VISUALIZATION AND RESULTS ---
vectors <- convergence_results$vectors
critical_value <- convergence_results$critical_value
chi_squared_stats <- convergence_results$chi_squared_stats

# Calculate the critical radius for the confidence sphere
critical_radius <- sqrt(critical_value)

# Generate sphere coordinates
sphere_coords <- generate_sphere(critical_radius)

# Create the 3D plot with Plotly
fig <- plot_ly(x = vectors[,1], y = vectors[,2], z = vectors[,3], 
               type = 'scatter3d', mode = 'markers',
               name = 'Z-statistic',
               marker = list(color = ifelse(chi_squared_stats > critical_value, 'red', 'blue'),
                             size = 4)) %>%
  add_trace(x = sphere_coords$x, y = sphere_coords$y, z = sphere_coords$z,
            type = 'scatter3d', mode = 'markers', 
            name = 'Confidence Region',
            marker = list(color = 'lightblue', size = 2, opacity = 0.2),
            showlegend = FALSE, hoverinfo = 'none') %>%
  layout(title = "Convergence of the Z-statistic",
         scene = list(xaxis = list(title = 'Z1 (alpha)'),
                      yaxis = list(title = 'Z2 (beta)'),
                      zaxis = list(title = 'Z3 (sigma)')))

# Show the plot
fig

# --- COVERAGE REPORT ---
# Calculate the confidence interval for the coverage probability
z_alpha_conf <- qnorm(1 - alpha_level / 2)
standard_error <- sqrt(alpha_level * (1 - alpha_level) / M)
lower_bound <- (1 - alpha_level) - z_alpha_conf * standard_error
upper_bound <- (1 - alpha_level) + z_alpha_conf * standard_error

cat("\n--- Convergence Analysis Results ---\n")
cat("Observed coverage probability:", convergence_results$coverage_probability, "\n")
cat("Theoretical confidence level:", 1 - alpha_level, "\n")
cat(paste0((1-alpha_level)*100, "% Confidence Interval for coverage: [", round(lower_bound, 4), ", ", round(upper_bound, 4), "]\n"))