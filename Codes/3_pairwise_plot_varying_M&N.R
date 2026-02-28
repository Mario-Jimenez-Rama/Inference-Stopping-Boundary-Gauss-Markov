#' Generate Pairwise Asymptotic Normality Plots for Custom Scenarios
#'
#' This function performs a simulation study to check asymptotic normality.
#' It compares different scenarios, each defined by a number of time steps (N)
#' and a number of Monte Carlo paths (M).
#'
#' Key changes in this version:
#' 1. Fisher matrices are calculated for EACH path using the ESTIMATED parameters for that path.
#' 2. These individual Fisher matrices are then averaged to get the final normalization matrix.
#' 3. The function takes a data.frame of scenarios to compare (e.g., high N/M vs. low N/M).
#'
#' @param theta_true Vector of true parameter values (alpha, beta, sigma).
#' @param scenarios_df A data.frame with columns 'N' and 'M' defining the scenarios to run.
#' @param delta Time step size.
#' @param subdivision Number of subdivisions for numerical integration.
#' @param theta_test Initial guess for the estimation procedure.
#' @param ... Functions for slope, trend, vol, and their derivatives.
#' @return A ggpairs plot showing the results.
pairwise_plot_custom_scenarios <- function(theta_true, scenarios_df, delta, subdivision, theta_test, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                           dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3){
  
  # List to hold the final data frames for each scenario
  all_data <- list()
  
  # --- Main loop over the defined scenarios ---
  for (scenario_idx in 1:nrow(scenarios_df)) {
    
    N <- scenarios_df$N[scenario_idx]
    M <- scenarios_df$M[scenario_idx]
    
    cat(paste0("\n--- Processing Scenario ", scenario_idx, ": N = ", N, ", M = ", M, " ---\n"))
    T <- N * delta
    
    # --- Step 1: Data Collection Loop ---
    # In this loop, we run M simulations. For each, we estimate parameters
    # and calculate its individual Fisher matrix based on those estimates.
    
    estimated_coefs_list <- list()
    individual_fisher_matrices <- list()
    
    cat("Step 1: Simulating paths, estimating parameters, and calculating individual Fisher matrices...\n")
    for (i in 1:M) {
      if (i %% 20 == 0) cat("  Processing path", i, "of", M, "\r")
      
      # 1a. Simulate a SINGLE path
      # We use the functions evaluated at the TRUE parameters for simulation
      theta_true_list <- list(alpha = theta_true[1], beta = theta_true[2], sigma = theta_true[3])
      slope_fun_t_true <- function(t) do.call(slope_fun, c(list(t = t), theta_true_list))
      trend_fun_t_true <- function(t) do.call(trend_fun, c(list(t = t), theta_true_list))
      vol_fun_t_true <- function(t) do.call(vol_fun, c(list(t = t), theta_true_list))
      
      sim_result <- simulate_gm_process(n_paths = 1, n_steps = N, delta = delta, x0 = 1, 
                                        slope_fun_t_true, trend_fun_t_true, vol_fun_t_true, 
                                        theta_true[1], theta_true[2], theta_true[3], subdivision)
      
      X <- sim_result$paths
      if (any(is.na(X))) next # Skip if simulation failed
      
      # 1b. Estimate parameters for the current path
      
      # maximization method
      result <- tryCatch({
        lk_estimates(subdivision, X, T, theta_test, slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
      }, error = function(e) {
        message("Error in estimation for path ", i, " (N=", N, "): ", e$message)
        NULL
      })
      
      if (is.null(result)) next # Skip if estimation failed
      
      estimated_coefs <- result@coef
      
      # gradient method
      # result <- tryCatch({
      #   MLE_grad_loglikelihood(initial_guess = c(alpha_test, beta_test, sigma_test),T, X, subdivision, slope_fun = slope_fun,
      #                          trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2,
      #                          dslope_dtheta3 = dslope_dtheta3, dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
      #                          dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3, method = "Newton", global = "cline")
      # }, error = function(e) {
      #   message("Error in estimation for path ", i, " (N=", N, "): ", e$message)
      #   NULL
      # })
      # 
      # 
      # if (is.null(result)) next # Skip if estimation failed
      # 
      # estimated_coefs <- result$x
      
      # 1c. Calculate Fisher matrix for THIS path using ITS estimated parameters
      params_list_est <- list(alpha = estimated_coefs[1], beta = estimated_coefs[2], sigma = estimated_coefs[3])
      
      # Define functions based on the ESTIMATED parameters
      slope_fun_t_est <- function(t) do.call(slope_fun, c(list(t = t), params_list_est))
      trend_fun_t_est <- function(t) do.call(trend_fun, c(list(t = t), params_list_est))
      vol_fun_t_est <- function(t) do.call(vol_fun, c(list(t = t), params_list_est))
      dslope_dtheta1_t_est <- function(t) do.call(dslope_dtheta1, c(list(t = t), params_list_est))
      dslope_dtheta2_t_est <- function(t) do.call(dslope_dtheta2, c(list(t = t), params_list_est))
      dslope_dtheta3_t_est <- function(t) do.call(dslope_dtheta3, c(list(t = t), params_list_est))
      dtrend_dtheta1_t_est <- function(t) do.call(dtrend_dtheta1, c(list(t = t), params_list_est))
      dtrend_dtheta2_t_est <- function(t) do.call(dtrend_dtheta2, c(list(t = t), params_list_est))
      dtrend_dtheta3_t_est <- function(t) do.call(dtrend_dtheta3, c(list(t = t), params_list_est))
      dvol_dtheta1_t_est <- function(t) do.call(dvol_dtheta1, c(list(t = t), params_list_est))
      dvol_dtheta2_t_est <- function(t) do.call(dvol_dtheta2, c(list(t = t), params_list_est))
      dvol_dtheta3_t_est <- function(t) do.call(dvol_dtheta3, c(list(t = t), params_list_est))
      
      fisher_info_i <- Fisher_matrix_der_exact_complete(T, X, subdivision, 
                                                        slope_fun = slope_fun_t_est, trend_fun = trend_fun_t_est, vol_fun = vol_fun_t_est,
                                                        dslope_dtheta1 = dslope_dtheta1_t_est, dslope_dtheta2 = dslope_dtheta2_t_est, dslope_dtheta3 = dslope_dtheta3_t_est,
                                                        dtrend_dtheta1 = dtrend_dtheta1_t_est, dtrend_dtheta2 = dtrend_dtheta2_t_est, dtrend_dtheta3 = dtrend_dtheta3_t_est,
                                                        dvol_dtheta1 = dvol_dtheta1_t_est, dvol_dtheta2 = dvol_dtheta2_t_est, dvol_dtheta3 = dvol_dtheta3_t_est)
      
      # 1d. Store results from this path
      if (!is.null(fisher_info_i) && !any(is.na(fisher_info_i$Fisher_matrix))) {
        estimated_coefs_list[[length(estimated_coefs_list) + 1]] <- estimated_coefs
        individual_fisher_matrices[[length(individual_fisher_matrices) + 1]] <- fisher_info_i$Fisher_matrix
      }
    }
    cat("\nFinished data collection for this scenario.\n")
    
    # --- Step 2: Process collected data for the scenario ---
    
    n_complete <- length(estimated_coefs_list)
    cat("Step 2: Processing results for", n_complete, "successfully completed paths.\n")
    if (n_complete == 0) {
      warning("No paths were successfully processed for this scenario. Skipping.")
      next
    }
    
    # 2a. Average the collected Fisher matrices
    avg_fisher_matrix <- Reduce('+', individual_fisher_matrices) / n_complete
    
    # 2b. Cholesky decomposition of the averaged matrix
    R <- tryCatch({
      chol(avg_fisher_matrix)
    }, error = function(e) {
      message("Cholesky decomposition failed for this scenario: ", e$message)
      message("The averaged Fisher matrix was:\n"); print(avg_fisher_matrix)
      NULL
    })
    
    if (is.null(R)) {
      warning("Skipping this scenario due to non-positive definite Fisher matrix.")
      next
    }
    
    # 2c. Normalize all stored estimates using the single Cholesky matrix R
    data_scenario <- list()
    for (coefs in estimated_coefs_list) {
      point <- sqrt(N) * (R %*% (coefs - theta_true))
      data_scenario[[length(data_scenario) + 1]] <- setNames(as.numeric(point), c("alpha", "beta", "sigma"))
    }
    
    # Convert list to data frame and add scenario identifier
    df_scenario <- do.call(rbind, data_scenario) |> as.data.frame()
    df_scenario$scenario <- paste0("N=", N, ", M=", M)
    all_data[[scenario_idx]] <- df_scenario
    
  } # End of loop over scenarios
  
  # --- Step 3: Combine all data and create plots ---
  if (length(all_data) == 0) {
    message("No data was generated across all scenarios. Cannot create plot.")
    return(invisible(NULL))
  }
  
  combined_data <- do.call(rbind, all_data)
  
  # Custom plot functions for ggpairs
  density_with_normal <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(aes(color = scenario), alpha = 0.7) + 
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linetype = "dashed") +
      theme_minimal() +
      xlim(c(-6, 6))
  }
  
  scatterplots <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(aes(color = scenario), size = 0.8, alpha = 0.7) +
      theme_minimal() +
      xlim(c(-6, 6)) +
      ylim(c(-6, 6))
  }
  
  # Create the final plot
  p <- ggpairs(combined_data, columns = 1:3, aes(color = scenario),
               upper = list(continuous = "cor"),
               lower = list(continuous = wrap(scatterplots)),
               diag = list(continuous = wrap(density_with_normal))
  )
  
  print(p)
  return(p)
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
# ## derivadas ,cuando tiene t esque solo depende de t y no del parametro
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

#Hull White model periodic mean dificil
# alpha_true <- 2
# beta_true <- 0.5
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
# ## derivadas ,cuando tiene t esque solo depende de t y no del parametro
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
# # derivadas ,cuando tiene t esque solo depende de t y no del parametro
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
# alpha_true <- -1 #mean reversion
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
# ## derivadas ,cuando tiene t esque solo depende de t y no del parametro
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
alpha_true <- -1 #mean reversion
beta_true <- 3
sigma_true <- 2
slope_fun <- function(t, alpha,...) rep(alpha, length(t))
trend_fun <- function(t, beta,...) rep(beta, length(t))
vol_fun <- function(t, sigma,...) rep(sigma, length(t))

# functions with true parameters
slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)

## derivadas
# derivadas de slope:
dslope_dtheta1 <- function(t,...) rep(1, length(t)) # cuando tiene t esque solo depende de t y no del parametro
dslope_dtheta2 <- function(t,...) rep(0, length(t))
dslope_dtheta3 <- function(t,...) rep(0, length(t))

# derivadas de trend
dtrend_dtheta1 <- function(t,...) rep(0, length(t))
dtrend_dtheta2 <- function(t,...) rep(1, length(t))
dtrend_dtheta3 <- function(t,...) rep(0, length(t))

# derivadas de vol
dvol_dtheta1 <- function(t,...) rep(0, length(t))
dvol_dtheta2 <- function(t,...) rep(0, length(t))
dvol_dtheta3 <- function(t,...) rep(1, length(t))

# Hull White periodic mean version 1
# alpha_true <- -1 # mean reversion
# beta_true <- 3
# sigma_true <- 1
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta, sigma,...) beta + sin(sigma*t)
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
# dtrend_dtheta3 <- function(t, sigma,...) t*cos(sigma*t)
# 
# # derivadas de vol
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(0, length(t))

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
# # vol_fun <- function(t,sigma,...) t/(1 + sigma*t^2)
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

delta = 1 #time step
subdiv <- 100 # depende de T y de n, para OU no importa ya que las funciones son ctes

# --- 1. Define the scenarios you want to compare ---
# Scenario 1: High N, High M
# Scenario 2: Low N, Low M
scenarios <- data.frame(
  N = c(200, 50),  # Number of time steps
  M = c(200, 75)  # Number of Monte Carlo paths
)

alpha_test <- -0.5 # mean reversion
beta_test <- 0.5
sigma_test <- 1

pairwise_plot_custom_scenarios(
  theta_true = c(alpha_true, beta_true, sigma_true),
  scenarios_df = scenarios,
  delta = delta,
  subdivision = subdiv,
  theta_test = c(alpha_test, beta_test, sigma_test),
  slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun,
  dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3,
  dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
  dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3
)

