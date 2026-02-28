#================================================================================#
# MAIN FUNCTION FOR ROBUST INFERENCE OF BOUNDARIES
#================================================================================#

#' @title Performs robust inference of stopping boundaries for M trajectories.
#' @description Encapsulates a 4-phase workflow: simulation, estimation,
#' Fisher Information Matrix averaging, and inference using the robust FIM.
#' @return An invisible list with the averaged inverse FIM and the results of all inferences.
infer_boundaries_robust <- function(M, N, delta, x0, theta_true, theta_test,
                                    slope_fun, trend_fun, vol_fun,
                                    dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                    dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                    dvol_dtheta1, dvol_dtheta2, dvol_dtheta3,
                                    training_points, partition_length,
                                    strike, discount, z_alpha,
                                    der_tol, subdivision,
                                    save_path, verbose = TRUE) {
  
  # --- PHASE 1: SIMULATION ---
  if (verbose) cat("--- PHASE 1: Simulating", M, "complete trajectories ---\n")
  
  params_true <- list(alpha = theta_true[1], beta = theta_true[2], sigma = theta_true[3])
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params_true))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params_true))
  vol_fun_t   <- function(t) do.call(vol_fun,  c(list(t = t), params_true))
  
  all_paths <- simulate_gm_process(
    n_paths = M, n_steps = N, delta = delta, x0 = x0,
    slope_fun_t, trend_fun_t, vol_fun_t,
    theta_true[1], theta_true[1], theta_true[1], subdivision
  )$paths
  
  # --- PHASE 2: TRAINING (INDIVIDUAL ESTIMATION AND FIM) ---
  if (verbose) cat("--- PHASE 2: Estimating parameters and FIM for each of the", M, "trajectories ---\n")
  
  all_estimates <- vector("list", M)
  all_fims <- vector("list", M)
  training_duration <- training_points * delta
  
  for (i in 1:M) {
    if (verbose && i %% 20 == 0) cat("  Processing training trajectory:", i, "/", M, "\n")
    
    # training_points + 1 because the first point x0 is included
    X_training <- all_paths[1:(training_points + 1), i, drop = FALSE]
    
    # MLE method
    result <- tryCatch(
      lk_estimates(subdivision, X_training, training_duration, theta_test = theta_test,
                   slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun),
      error = function(e) { message("Error in estimation (trajectory ", i, "): ", e$message); NULL }
    )
    
    if (is.null(result)) next
    
    theta_estimated <- as.numeric(result@coef)
    
    # gradient method
    # result <- tryCatch({
    #   MLE_grad_loglikelihood(initial_guess = c(alpha_test, beta_test, sigma_test), training_duration, X_training, subdivision,
    #                          slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun,
    #                          dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3,
    #                          dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
    #                          dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3,
    #                          method = "Newton", global = "cline")
    # }, error = function(e) {
    #   message("Error in estimation for path ", i, " (N=", N, "): ", e$message)
    #   NULL
    # })
    #
    # if (is.null(result)) next # Skip if estimation failed
    #
    # theta_estimated <- as.numeric(result$x)
    
    all_estimates[[i]] <- theta_estimated
    
    params_est <- list(alpha = theta_estimated[1], beta = theta_estimated[2], sigma = theta_estimated[3])
    
    # Calculate the FIM for this trajectory
    fim_i <- Fisher_matrix_der_exact_complete(
      training_duration, X_training, subdivision,
      slope_fun = function(t) do.call(slope_fun, c(list(t = t), params_est)),
      trend_fun = function(t) do.call(trend_fun, c(list(t = t), params_est)),
      vol_fun   = function(t) do.call(vol_fun,   c(list(t = t), params_est)),
      dslope_dtheta1 = function(t) do.call(dslope_dtheta1, c(list(t = t), params_est)),
      dslope_dtheta2 = function(t) do.call(dslope_dtheta2, c(list(t = t), params_est)),
      dslope_dtheta3 = function(t) do.call(dslope_dtheta3, c(list(t = t), params_est)),
      dtrend_dtheta1 = function(t) do.call(dtrend_dtheta1, c(list(t = t), params_est)),
      dtrend_dtheta2 = function(t) do.call(dtrend_dtheta2, c(list(t = t), params_est)),
      dtrend_dtheta3 = function(t) do.call(dtrend_dtheta3, c(list(t = t), params_est)),
      dvol_dtheta1   = function(t) do.call(dvol_dtheta1,   c(list(t = t), params_est)),
      dvol_dtheta2   = function(t) do.call(dvol_dtheta2,   c(list(t = t), params_est)),
      dvol_dtheta3   = function(t) do.call(dvol_dtheta3,   c(list(t = t), params_est))
    )$Fisher_matrix
    
    all_fims[[i]] <- fim_i
  }
  
  # --- PHASE 3: AGGREGATION (AVERAGE FIM) ---
  if (verbose) cat("--- PHASE 3: Averaging the FIMs to obtain a robust matrix ---\n")
  
  valid_fims <- all_fims[!sapply(all_fims, is.null)]
  if (length(valid_fims) == 0) stop("Could not calculate any Fisher Matrix. All estimations failed.")
  
  fisher_info_averaged <- Reduce('+', valid_fims) / length(valid_fims)
  
  fisher_info_inv_averaged <- tryCatch(
    solve(fisher_info_averaged),
    error = function(e) stop("Error inverting the averaged FIM. Error: ", e$message)
  )
  
  # --- PHASE 4: INFERENCE (USING AVERAGED FIM) ---
  if (verbose) cat("--- PHASE 4: Inferring boundaries and their CIs for each trajectory ---\n")
  
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  all_results_list <- vector("list", M)
  
  for (i in 1:M) {
    if (verbose && i %% 20 == 0) cat("  Processing inference trajectory:", i, "/", M, "\n")
    
    theta_estimated <- all_estimates[[i]]
    if (is.null(theta_estimated)) next
    
    expiration <- delta * partition_length
    boundary <- boundary_wrapper(theta = theta_estimated, strike = strike, partition_length = partition_length - 1,
                                 discount = discount, expiration = expiration, slope_fun, trend_fun, vol_fun)
    
    h_alpha <- der_tol * abs(theta_estimated[1])
    h_beta  <- der_tol * abs(theta_estimated[2])
    h_sigma <- der_tol * abs(theta_estimated[3])
    
    g1 <- (boundary_wrapper(c(theta_estimated[1] + h_alpha, theta_estimated[2:3]), strike, expiration,
                            partition_length - 1, discount, slope_fun, trend_fun, vol_fun) -
             boundary_wrapper(c(theta_estimated[1] - h_alpha, theta_estimated[2:3]), strike, expiration,
                              partition_length - 1, discount, slope_fun, trend_fun, vol_fun)) / (2 * h_alpha)
    
    g2 <- (boundary_wrapper(c(theta_estimated[1], theta_estimated[2] + h_beta, theta_estimated[3]), strike, expiration,
                            partition_length - 1, discount, slope_fun, trend_fun, vol_fun) -
             boundary_wrapper(c(theta_estimated[1], theta_estimated[2] - h_beta, theta_estimated[3]), strike, expiration,
                              partition_length - 1, discount, slope_fun, trend_fun, vol_fun)) / (2 * h_beta)
    
    g3 <- (boundary_wrapper(c(theta_estimated[1:2], theta_estimated[3] + h_sigma), strike, expiration,
                            partition_length - 1, discount, slope_fun, trend_fun, vol_fun) -
             boundary_wrapper(c(theta_estimated[1:2], theta_estimated[3] - h_sigma), strike, expiration,
                              partition_length - 1, discount, slope_fun, trend_fun, vol_fun)) / (2 * h_sigma)
    
    grad_matrix <- as.matrix(data.frame(g1, g2, g3))
    
    # averaged FIM
    # variances <- rowSums((grad_matrix %*% fisher_info_inv_averaged) * grad_matrix) / training_points
    
    # individual FIM
    fim_i_inv <- solve(all_fims[[i]])
    variances <- rowSums((grad_matrix %*% fim_i_inv) * grad_matrix) / training_points
    
    std_dev <- sqrt(variances)
    
    z_critical <- qnorm(1 - z_alpha / 2)
    upper_bound <- boundary + z_critical * std_dev
    lower_bound <- boundary - z_critical * std_dev
    
    inference_result <- list(lower_bound = lower_bound,
                             boundary_est = boundary,
                             upper_bound = upper_bound,
                             est_theta = theta_estimated)
    all_results_list[[i]] <- inference_result
    
    # Save the specified file index
    file_index <- i
    folder_name <- "rds_files"
    
    # Define file paths
    x_path_file <- file.path(folder_name, paste0("path_", save_path, "_", file_index, ".rds"))
    inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_", save_path, "_", file_index, ".rds"))
    
    # Save the data to .rds files
    saveRDS(all_paths[, i], x_path_file)
    saveRDS(inference_result, inferred_boundary_file)
  }
  
  if (verbose) cat("--- Process completed. Results saved in '", save_path, "' ---\n")
  
  return(invisible(list(
    averaged_fim_inverse = fisher_info_inv_averaged,
    all_inference_results = all_results_list,
    all_paths = all_paths
  )))
}

# OUB (commented out)

# time-dependent rational function volatility OU
alpha_true <- -1
beta_true <- 3
sigma_true <- 0.5
slope_fun <- function(t, alpha, ...) rep(alpha, length(t))
trend_fun <- function(t, beta, ...) rep(beta, length(t))
vol_fun   <- function(t, sigma, ...) 1/(1 + sigma * t)
# vol_fun <- function(t, sigma, ...) t/(1 + sigma * t^2)

# functions with true parameters
slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
vol_fun_t   <- function(t) vol_fun(t,   alpha = alpha_true, beta = beta_true, sigma = sigma_true)

## derivatives
# slope derivatives:
dslope_dtheta1 <- function(t, ...) rep(1, length(t))
dslope_dtheta2 <- function(t, ...) rep(0, length(t))
dslope_dtheta3 <- function(t, ...) rep(0, length(t))

# trend derivatives
dtrend_dtheta1 <- function(t, ...) rep(0, length(t))
dtrend_dtheta2 <- function(t, ...) rep(1, length(t))
dtrend_dtheta3 <- function(t, ...) rep(0, length(t))

# volatility derivatives
dvol_dtheta1 <- function(t, ...) rep(0, length(t))
dvol_dtheta2 <- function(t, ...) rep(0, length(t))
dvol_dtheta3 <- function(t, sigma, ...) -t/(1 + sigma * t)^2

# --- 2. Define all parameters for the run ---
run_params <- list(
  M = 200,
  N = 100,
  delta = 1,
  x0 = 10,
  theta_true = c(alpha = alpha_true, beta = beta_true, sigma = sigma_true),
  theta_test = c(alpha = -0.5, beta = 1, sigma = 0.5),
  slope_fun = slope_fun,
  trend_fun = trend_fun,
  vol_fun = vol_fun,
  dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3,
  dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
  dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3,
  training_points = 75,
  partition_length = 25,
  strike = 9,
  discount = 0,
  z_alpha = 0.1,
  der_tol = 1e-2,
  subdivision = 100,
  save_path = "rational", # OU, rational, irrational, logarithmic, or periodic
  verbose = TRUE
)

# --- 3. Call the main function with all parameters ---
results <- do.call(infer_boundaries_robust, run_params)