library(pracma)
library(sde)
library(ggplot2)
library(gridExtra)
library(GGally)
library(stats4)
library(nleqslv)
library(ggtext)

# --- Helper function for direct integration with trapz (unchanged) ---
integrate_direct_trapz <- function(func_to_integrate, lower, upper, subdivisions, ...) {
  if (abs(lower - upper) < .Machine$double.eps^0.75) {
    return(0)
  }
  if (subdivisions < 2) {
    stop("integrate_direct_trapz requires at least 2 subdivisions (points).")
  }
  grid <- seq(lower, upper, length.out = subdivisions)
  
  func_name_str <- deparse(substitute(func_to_integrate))
  # If func_to_integrate is one of those that expects a scalar as the first grid argument.
  if (grepl("integrand_for_trend_term_optimized", func_name_str) || 
      grepl("integrand_for_vol_term_optimized", func_name_str) ) { 
    values <- sapply(grid, func_to_integrate, ...)
  } else { # For functions like slope, trend, vol that are assumed to be vectorized
    values <- func_to_integrate(grid, ...)
  }
  
  return(pracma::trapz(grid, values))
}

# --- OPTIMIZED integrand functions (use spline_F_slope) ---
integrand_for_trend_term_optimized <- function(u_mesh_val, current_t2, 
                                               spline_F_slope_func, # Changed: spline function
                                               trend_func) {
  # integral_slope_u_t2 is obtained from the precomputed spline
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  return(trend_func(u_mesh_val) * exp(integral_slope_u_t2))
}

integrand_for_vol_term_optimized <- function(u_mesh_val, current_t2, 
                                             spline_F_slope_func, # Changed: spline function
                                             vol_func) {
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  term_in_paren <- vol_func(u_mesh_val) * exp(integral_slope_u_t2)
  return(term_in_paren^2)
}

# --- Kernel with OPTIMIZED integral computation ---
boundary_kernel_optimized <- function(c1, c2, t1, x1, t2, x2, 
                                      spline_F_slope_func, # Changed: spline function
                                      trend_f, vol_f, 
                                      discount_val, subdivisions,
                                      # No longer needs slope_f directly here for integrals
                                      # But it is needed for logging if you want slope_f(t1)
                                      slope_f_for_log_only # Only for logging
) {
  
  # 1. exp_int_slope_t1_t2 = exp(integral_{t1}^{t2} slope(u)du)
  # Obtained from the precomputed spline
  integral_slope_t1_t2 <- spline_F_slope_func(t2) - spline_F_slope_func(t1)
  exp_int_slope_t1_t2 <- exp(integral_slope_t1_t2)
  
  # 2. integral_trend_term 
  integral_trend_term <- integrate_direct_trapz(
    integrand_for_trend_term_optimized, t1, t2, subdivisions,
    current_t2 = t2, spline_F_slope_func = spline_F_slope_func, trend_func = trend_f
  )
  
  # 3. integral_vol_term (direct variance)
  marginal_var_val_direct <- integrate_direct_trapz(
    integrand_for_vol_term_optimized, t1, t2, subdivisions,
    current_t2 = t2, spline_F_slope_func = spline_F_slope_func, vol_func = vol_f
  )
  
  marginal_mean_val <- x1 * exp_int_slope_t1_t2 + integral_trend_term
  
  current_marginal_var <- marginal_var_val_direct
  if (current_marginal_var < 0 && abs(current_marginal_var) < .Machine$double.eps*100) {
    current_marginal_var <- 0
  } else if (current_marginal_var < 0) {
    current_marginal_var <- 0 
  }
  marginal_sd_val <- sqrt(current_marginal_var)
  
  # --- Logging ---
  # line <- paste(c1, c2, t1, x1, t2, x2,
  #               slope_f_for_log_only(t1), trend_f(t1), vol_f(t1),
  #               marginal_mean_val, current_marginal_var, sep = "\t")
  # try(write(line, file = "stopping_boundary_Abel.txt", append = TRUE, ncolumns = 11), silent=TRUE)
  
  # --- Kernel value computation ---
  x2_std <- 0.0
  if (marginal_sd_val > .Machine$double.eps^0.5) {
    x2_std <- (x2 - marginal_mean_val) / marginal_sd_val
  } else {
    diff_val_std <- x2 - marginal_mean_val
    x2_std <- ifelse(abs(diff_val_std) < .Machine$double.eps^0.5, 0.0, ifelse(diff_val_std > 0, Inf, -Inf))
  }
  
  normal_dist <- pnorm(x2_std, mean = 0, sd = 1, lower.tail = TRUE)
  normal_dens <- dnorm(x2_std, mean = 0, sd = 1)
  
  term1 <- (c1 - c2 * marginal_mean_val) * normal_dist
  term2 <- c2 * marginal_sd_val * normal_dens
  
  K_val <- exp(-discount_val * (t2 - t1)) * (term1 + term2)
  
  return(K_val)
}


# --- Main boundary function (OPTIMIZED) ---
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0,
                      slope, trend, vol, errors = FALSE,
                      trapz_subdivs = 100, 
                      cumtrapz_subdivs_slope = 1000 # For the slope spline
) {
  
  # --- Pre-computation of the slope integral using cumtrapz and spline ---
  # Fine grid from the first point of time_line to expiration
  # Make sure min(time_line) is time_line[1] if it is ordered
  min_time <- min(time_line) # Use min/max for robustness if time_line is not ordered
  max_time <- max(time_line)
  
  # Ensure cumtrapz_subdivs_slope is at least 2
  safe_cumtrapz_subdivs_slope <- max(2, cumtrapz_subdivs_slope)
  fine_grid_for_slope_spline <- seq(min_time, max_time, length.out = safe_cumtrapz_subdivs_slope)
  
  slope_vals_on_fine_grid <- slope(fine_grid_for_slope_spline)
  # Ensure slope_vals_on_fine_grid is a vector of the correct size
  if(!is.vector(slope_vals_on_fine_grid) || length(slope_vals_on_fine_grid) != length(fine_grid_for_slope_spline)) {
    # If slope is not vectorized, apply sapply.
    # This is important if slope is not of the form rep(const, length(t))
    slope_vals_on_fine_grid <- sapply(fine_grid_for_slope_spline, slope)
  }
  
  # F(x) = integral_{fine_grid_for_slope_spline[1]}^{x} slope(u) du
  cumulative_integral_slope_values <- pracma::cumtrapz(fine_grid_for_slope_spline, slope_vals_on_fine_grid)
  
  # Create spline function for F(x)
  spline_F_slope <- splinefun(fine_grid_for_slope_spline, cumulative_integral_slope_values, method = "natural")
  # --- End pre-computation ---
  
  # Write header to the file (if logging is active in the kernel)
  # write("c1\tc2\tt1\tx1\tt2\tx2\tslope(t1)\ttrend(t1)\tvol(t1)\tmarginal_mean\tmarginal_var", file = "stopping_boundary_Abel.txt")
  
  N <- length(time_line)
  expiration <- time_line[N] # Assumes the last point is the expiration
  delta_t_steps <- time_line[2:N] - time_line[1:(N-1)]
  
  if (errors) er <- c()
  
  initial_bnd_val <- strike
  # For slope(expiration), it is better to evaluate the original function, not the spline.
  slope_at_expiration <- slope(expiration)[1] # Take the first element if it returns a vector
  trend_at_expiration <- trend(expiration)[1]
  
  if (abs(discount - slope_at_expiration) > .Machine$double.eps^0.75) {
    initial_bnd_val <- min((trend_at_expiration + discount * strike) / (discount - slope_at_expiration), strike)
  } else {
    # warning("Denominator (discount - slope(expiration)) is close to zero in boundary init.")
  }
  bnd <- rep(initial_bnd_val, N)
  
  e <- 1.0
  j <- 0
  last_valid_bnd <- bnd
  
  while (e > tol) {
    j <- j + 1
    bnd_old <- bnd
    
    print(paste("Picard iteration:", j)) # To see progress
    print(paste("error:", e)) # To see progress
    
    for (i in (N - 1):1) {
      t_current <- time_line[i]
      b_current_old <- bnd_old[i]
      
      K1 <- boundary_kernel_optimized(c1 = strike, c2 = 1,
                                      t1 = t_current, x1 = b_current_old,
                                      t2 = expiration, x2 = strike,
                                      spline_F_slope_func = spline_F_slope, 
                                      trend_f = trend, vol_f = vol,
                                      discount_val = discount, subdivisions = trapz_subdivs,
                                      slope_f_for_log_only = slope) 
      
      sum_K2_delta_val <- 0.0
      if (i < N) {
        idx_future_points <- (i + 1):N
        times_future <- time_line[idx_future_points]
        b_future_old_vals <- bnd_old[idx_future_points]
        
        trend_vals_future <- trend(times_future)
        if(!is.vector(trend_vals_future) && length(times_future)>0) trend_vals_future <- sapply(times_future, trend)
        slope_vals_future <- slope(times_future) # Not the spline, the original function
        if(!is.vector(slope_vals_future) && length(times_future)>0) slope_vals_future <- sapply(times_future, slope)
        
        c1_K2_vec <- discount * strike + trend_vals_future
        c2_K2_vec <- discount - slope_vals_future
        
        K2_vec <- numeric(length(times_future))
        for (k_fut in 1:length(times_future)) {
          u_k <- times_future[k_fut]
          b_at_u_k <- b_future_old_vals[k_fut]
          c1_val_k <- c1_K2_vec[k_fut]
          c2_val_k <- c2_K2_vec[k_fut]
          
          K2_vec[k_fut] <- boundary_kernel_optimized(
            c1 = c1_val_k, c2 = c2_val_k,
            t1 = t_current, x1 = b_current_old,
            t2 = u_k, x2 = b_at_u_k,
            spline_F_slope_func = spline_F_slope,
            trend_f = trend, vol_f = vol,
            discount_val = discount, subdivisions = trapz_subdivs,
            slope_f_for_log_only = slope
          )
        }
        sum_K2_delta_val <- sum(K2_vec * delta_t_steps[i:(N-1)])
      }
      
      bnd[i] <- strike - K1 - sum_K2_delta_val
      
      if(bnd[i] > strike && !isTRUE(all.equal(bnd[i], strike, tolerance = tol*1e-2))){
        bnd[i] <- strike 
      } else if (bnd[i] > strike) { 
        bnd[i] <- strike
      }
    }
    
    e <- sum((bnd - bnd_old)^2, na.rm = FALSE) 
    
    if(is.na(e) || is.infinite(e) || is.nan(e)) {
      warning(paste("Error became NA/NaN/Inf at iteration",j,". Stopping. Using last valid boundary."))
      bnd <- last_valid_bnd 
      if (errors) er <- c(er, e) 
      break 
    }
    if(j > 1 || (j==1 && !(is.na(e) || is.infinite(e) || is.nan(e)))) { # Ensure bnd_old is 'good'
      last_valid_bnd <- bnd_old 
    }
    
    
    if (errors) er <- c(er, e)
    
    # Reduce the iteration limit if it is suspected to be very slow initially
    iter_limit <- 2000 
    if (j > iter_limit) { 
      warning(paste("Picard iteration limit (", iter_limit, ") reached. Error:", e))
      break
    }
  }
  
  print(paste0("Converged/Stopped after ", j, " iterations with error: ", sprintf("%.2e", e)))
  
  if (errors) return(list(boundary = bnd, errors = er))
  
  return(bnd)
}

boundary_wrapper <- function(theta, strike, expiration, partition_length, discount=0, slope_fun, trend_fun, vol_fun){
  
  time_grid <- seq(0, expiration, length.out = partition_length + 1)
  
  slope_fun_t <- function(t) slope_fun(t, alpha = theta[1])
  trend_fun_t <- function(t) trend_fun(t, beta = theta[2])
  vol_fun_t <- function(t) vol_fun(t, sigma = theta[3])
  
  bnd <- boundary(tol = 1e-3 ,strike = strike, time_line = time_line,
                  discount = discount, slope = slope_fun_t, trend = trend_fun_t,
                  vol = vol_fun_t, errors = TRUE, trapz_subdivs = 100, cumtrapz_subdivs_slope = 1000)
  
  return(bnd$boundary)
}

infer_boundary <- function(X, delta, z_alpha = 0.1, partition_length = 100, strike = 0, discount = 0, theta_test, slope_fun, trend_fun, vol_fun,dslope_dtheta1, dtrend_dtheta2, dvol_dtheta3, der_tol, subdivision) {
  
  expiration = delta*partition_length
  # likelihood4
  # here it would be
  result <- tryCatch({
    lk_estimates(subdivision, X, expiration, theta_test = theta_test,
                 slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
  }, error = function(e) {
    message("Estimation error: ", e$message)
    NULL
  })
  
  boundary <- boundary_wrapper(theta = result@coef, strike = strike, partition_length = partition_length,
                               discount = discount, expiration = expiration, slope_fun, trend_fun, vol_fun)
  
  h_alpha <- der_tol * abs(result@coef[1])
  h_beta <- der_tol * abs(result@coef[2])
  h_sigma <- der_tol * abs(result@coef[3])
  
  g_prime_theta1 <- (boundary_wrapper(c(alpha = result@coef[1] + h_alpha, beta = result@coef[2], sigma = result@coef[3]), strike = strike, expiration = expiration,  partition_length = partition_length, discount = discount, slope_fun, trend_fun, vol_fun) -
                       boundary_wrapper(c(alpha = result@coef[1] - h_alpha, beta = result@coef[2], sigma = result@coef[3]), strike = strike, expiration = expiration, partition_length = partition_length, discount= discount,  slope_fun, trend_fun, vol_fun)) / (2 * h_alpha)
  
  g_prime_theta2 <- (boundary_wrapper(c(alpha = result@coef[1], beta = result@coef[2] + h_beta, sigma = result@coef[3]), strike = strike, expiration = expiration, partition_length = partition_length, discount= discount, slope_fun, trend_fun, vol_fun) -
                       boundary_wrapper(c(alpha = result@coef[1], beta = result@coef[2] - h_beta, sigma = result@coef[3]), strike = strike, expiration = expiration, partition_length = partition_length, discount= discount, slope_fun, trend_fun, vol_fun)) / (2 * h_beta)
  
  g_prime_theta3 <- (boundary_wrapper(c(alpha = result@coef[1], beta = result@coef[2], sigma = result@coef[3] + h_sigma), strike = strike, expiration = expiration, partition_length = partition_length, discount= discount, slope_fun, trend_fun, vol_fun) -
                       boundary_wrapper(c(alpha = result@coef[1], beta = result@coef[2], sigma = result@coef[3] - h_sigma), strike = strike, expiration = expiration, partition_length = partition_length, discount= discount, slope_fun, trend_fun, vol_fun)) / (2 * h_sigma)
  
  grad_g <- data.frame(
    g_prime_theta1 = g_prime_theta1,
    g_prime_theta2 = g_prime_theta2,
    g_prime_theta3 = g_prime_theta3
  )
  # cat("g_prime_theta1 =", g_prime_theta1, "\n")
  # cat("g_prime_theta2 =", g_prime_theta2, "\n")
  # cat("g_prime_theta3 =", g_prime_theta3, "\n")
  
  slope_fun_t <- function(t) slope_fun(t, alpha = result@coef[1])
  trend_fun_t <- function(t) trend_fun(t, beta = result@coef[2])
  vol_fun_t <- function(t) vol_fun(t, sigma = result@coef[3])
  
  fisher_info <- tryCatch({
    Fisher_matrix_der_exact(T, X, subdivision, slope_fun = slope_fun_t,
                            trend_fun = trend_fun_t, vol_fun = vol_fun_t, dslope_dtheta1 = dslope_dtheta1,
                            dtrend_dtheta2 = dtrend_dtheta2, dvol_dtheta3 = dvol_dtheta3)$Fisher_matrix
  }, error = function(e) {
    message("Error computing Fisher matrix: ", e$message)
    NULL
  })
  fisher_info_inv <- solve(fisher_info)
  
  std_dev <- numeric(length = nrow(grad_g))
  for (i in 1:nrow(grad_g)) {
    grad <- as.numeric(grad_g[i, ])
    std_dev[i] <- sqrt(t(grad) %*% fisher_info_inv %*% grad)/sqrt(length(X))
  }
  
  upper_bound <- boundary + qnorm(1 - z_alpha/2) * std_dev
  lower_bound <- boundary - qnorm(1 - z_alpha/2) * std_dev
  
  return(list(lower_bound = lower_bound, boundary_est = boundary, upper_bound = upper_bound, est_theta = result@coef))
}

alpha_int_dt <- function(t0, t1, slope_fun, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions; could add an argument to vary the partition
  y <- slope_fun(time_grid)
  # cat("Dimensions of y:", if(is.matrix(y)) paste(dim(y), collapse="x") else length(y), "\n")
  # cat("Length time_grid:", length(time_grid), "\n")
  int_values <- cumtrapz(time_grid, y)
  # return(splinefun(time_grid, exp(int_values)) # exponential directly
  return(splinefun(time_grid, int_values)) # without taking the exponential directly; spline is not strictly needed if you only want values on time_grid
}

# alpha_int_dt <- function(t0, t1, slope_fun, subdiv){
#   time_grid <- seq(t0, t1, length.out = subdiv)
#   y <- slope_fun(time_grid)
#
#   # cat("Dimensions of y:", if(is.matrix(y)) paste(dim(y), collapse="x") else length(y), "\n")
#   # cat("Length time_grid:", length(time_grid), "\n")
#
#   # --- Key fix for gradient ---
#   if (length(y) == 3 * length(time_grid)) {  # If grad() perturbed 3 parameters
#     y <- matrix(y, nrow = length(time_grid), ncol = 3)  # Convert to a matrix (10x3)
#     int_values <- apply(y, 2, function(col) cumtrapz(time_grid, col))  # Integrate each column
#   } else {
#     stopifnot(length(y) == length(time_grid))  # Standard check
#     int_values <- cumtrapz(time_grid, y)  # Normal case (vector)
#   }
#
#   splinefun(time_grid, int_values) # builds a spline between t0 and t1
# }

trend_int_dt <- function(trend_fun, alpha_int_fun, t0, t1, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  integrand <- function(t) {
    return(trend_fun(t) / exp(alpha_int_fun(t)))
  }
  y <- integrand(time_grid)
  return(trapz(time_grid, y)) # no need for cumulative integration
}

vol_int_dt <- function(vol_fun, alpha_int_fun, t0, t1, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # time grid
  integrand <- function(t) {
    return((vol_fun(t) / exp(alpha_int_fun(t)))^2)
  }
  y <- integrand(time_grid)
  return(trapz(time_grid, y))
}

gm_mean_dt <- function(x0, t1, alpha_int_fun, trend_int){
  return(exp(alpha_int_fun(t1)) * (x0 + trend_int))
}

gm_var_dt <- function(t1, alpha_int_fun, vol_int){
  return(exp(2 * alpha_int_fun(t1)) * vol_int)
}

log_likelihood4 <- function(alpha, beta, sigma, T, gm_paths, subdiv, slope_fun, trend_fun, vol_fun){
  
  # cat("Class of X:", class(X), "Dimensions:", dim(X), "\n")
  # Convert various input formats to a matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    # Handle time-series objects
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Handle list of objects (can be ts, vectors, or matrices)
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # Handle single vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Handle matrices or data.frames directly
    gm_paths_matrix <- as.matrix(gm_paths)
  }
  
  # cat("alpha = ", alpha, "\n")
  # cat("beta = ", beta, "\n")
  # cat("sigma = ", sigma, "\n")
  # This computes the log-likelihood over all paths; that is not what I want.
  # I want the log-likelihood of each path.
  # Write header to file
  # write("x1\tt1\tx0\texp\ttrend\tvol\tmean\tvariance\tdnorm", file = "likelihood_mio4.txt")
  
  n <- nrow(gm_paths_matrix)
  num_paths <- ncol(gm_paths_matrix)
  time_grid <- seq(0, T, length.out = n) # length n because x0 aligns with x
  
  params <- list(alpha = alpha, beta = beta, sigma = sigma)
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t   <- function(t) do.call(vol_fun,   c(list(t = t), params))
  
  # OU stationary
  # slope_fun_t <- function(t) slope_fun(t, alpha = alpha)
  # trend_fun_t <- function(t) trend_fun(t, beta = beta)
  # vol_fun_t   <- function(t) vol_fun(t, sigma = sigma)
  
  # Brownian bridge
  # slope_fun_t <- function(t) -alpha/(T - t) # keep alpha constant
  # trend_fun_t <- function(t) beta/(T - t)
  # vol_fun_t   <- function(t) rep(sigma, length(t))
  
  # OU bridge alpha = alpha, z = beta, sigma same
  # xf <- gm_paths_matrix[n, num_paths]
  # trend_fun_t <- function(t) xf*alpha/sinh(alpha*(T-t)) - beta*tanh(alpha*(T-t)/2)
  # slope_fun_t <- function(t) -alpha/tanh(alpha*(T-t))
  # vol_fun_t   <- function(t) rep(sigma, length(t))
  
  total_logL <- 0
  
  # Process each path (column in the matrix)
  for (path_num in 1:num_paths) {
    gm_path <- gm_paths_matrix[, path_num]
    path_logL <- 0
    
    for (i in 2:n) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x1 <- gm_path[i]
      x0 <- gm_path[i - 1]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t, subdiv)
      exp_value   <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1, subdiv)
      vol_value   <- vol_int_dt(vol_fun_t,   alpha_int_fun, t0, t1, subdiv)
      mean_val    <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      var_val     <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      value <- dnorm(x1, mean_val, sqrt(var_val), log = TRUE)
      path_logL <- path_logL + value
      
      # Create a line with tab-separated values
      line <- paste(x1, t1, x0, exp_value, trend_value, vol_value, mean_val, var_val, value, sep = "\t")
      
      # Write to file, appending a new line each time
      write(line, file = "likelihood_mio4.txt", append = TRUE, ncolumns = 1)
    }
    
    total_logL <- total_logL + path_logL
  }
  
  # Return negative sum (for minimization)
  return(-total_logL)
}

lk_estimates <- function(subdivision, X, T, theta_test, slope_fun, trend_fun, vol_fun){
  opt_result <- mle(
    function(alpha, beta, sigma)
      log_likelihood4(alpha, beta, sigma, T, X, subdivision, slope_fun, trend_fun, vol_fun),
    start  = list(alpha = theta_test[1], beta = theta_test[2], sigma = theta_test[3]),
    method = "L-BFGS-B",
    lower  = c(-Inf, -Inf, 0.001)
  )
  return(opt_result)
}

dmean_dtheta1 <- function(trend_fun, alpha_int_fun, exp_value, trend_value,
                          dalpha_dtheta1, dtrend_dtheta1,
                          t0, t1, x0, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta1
  y1 <- dalpha_dtheta1(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta1_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return(trend_fun(t) * dalpha_dtheta1_int_fun(t) / exp(alpha_int_fun(t)))
  }
  y2 <- integrand2(time_grid)
  
  # integrand for the third term
  integrand3 <- function(t){
    return(dtrend_dtheta1(t) / exp(alpha_int_fun(t)))
  }
  y3 <- integrand3(time_grid)
  
  # Value of dmean_dtheta1
  return(exp_value * (dalpha_dtheta1_int_fun(t1) * (x0 + trend_value) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

dvar_dtheta1 <- function(exp_value, vol_value, vol_fun, alpha_int_fun,
                         dalpha_dtheta1, dvol_dtheta1,
                         t0, t1, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta1
  y1 <- dalpha_dtheta1(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta1_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return((vol_fun(t)^2) * dalpha_dtheta1_int_fun(t) / (exp(alpha_int_fun(t)))^2)
  }
  y2 <- integrand2(time_grid)
  
  integrand3 <- function(t){
    return(vol_fun(t) * dvol_dtheta1(t) / (exp(alpha_int_fun(t)))^2)
  }
  y3 <- integrand3(time_grid)
  
  # Value of dvar_dtheta1
  return(2 * exp_value^2 * (vol_value * dalpha_dtheta1_int_fun(t1) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

dmean_dtheta2 <- function(trend_fun, alpha_int_fun, exp_value, trend_value,
                          dalpha_dtheta2, dtrend_dtheta2,
                          t0, t1, x0, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta2
  y1 <- dalpha_dtheta2(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta2_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return(trend_fun(t) * dalpha_dtheta2_int_fun(t) / exp(alpha_int_fun(t)))
  }
  y2 <- integrand2(time_grid)
  
  # integrand for the third term
  integrand3 <- function(t){
    return(dtrend_dtheta2(t) / exp(alpha_int_fun(t)))
  }
  y3 <- integrand3(time_grid)
  
  # Value of dmean_dtheta2
  return(exp_value * (dalpha_dtheta2_int_fun(t1) * (x0 + trend_value) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

dvar_dtheta2 <- function(exp_value, vol_value, vol_fun, alpha_int_fun,
                         dalpha_dtheta2, dvol_dtheta2,
                         t0, t1, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta2
  y1 <- dalpha_dtheta2(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta2_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return((vol_fun(t)^2) * dalpha_dtheta2_int_fun(t) / (exp(alpha_int_fun(t)))^2)
  }
  y2 <- integrand2(time_grid)
  
  integrand3 <- function(t){
    return(vol_fun(t) * dvol_dtheta2(t) / (exp(alpha_int_fun(t)))^2)
  }
  y3 <- integrand3(time_grid)
  
  # Value of dvar_dtheta2
  return(2 * exp_value^2 * (vol_value * dalpha_dtheta2_int_fun(t1) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

dmean_dtheta3 <- function(trend_fun, alpha_int_fun, exp_value, trend_value,
                          dalpha_dtheta3, dtrend_dtheta3,
                          t0, t1, x0, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta3
  y1 <- dalpha_dtheta3(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta3_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return(trend_fun(t) * dalpha_dtheta3_int_fun(t) / exp(alpha_int_fun(t)))
  }
  y2 <- integrand2(time_grid)
  
  # integrand for the third term
  integrand3 <- function(t){
    return(dtrend_dtheta3(t) / exp(alpha_int_fun(t)))
  }
  y3 <- integrand3(time_grid)
  
  # Value of dmean_dtheta3
  return(exp_value * (dalpha_dtheta3_int_fun(t1) * (x0 + trend_value) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

dvar_dtheta3 <- function(exp_value, vol_value, vol_fun, alpha_int_fun,
                         dalpha_dtheta3, dvol_dtheta3,
                         t0, t1, subdiv){
  time_grid <- seq(t0, t1, l = subdiv) # computed with a fixed number of subdivisions
  
  # Integral of dalpha_dtheta3
  y1 <- dalpha_dtheta3(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta3_int_fun <- splinefun(time_grid, int_values)
  
  # integrand for the second term
  integrand2 <- function(t){
    return((vol_fun(t)^2) * dalpha_dtheta3_int_fun(t) / (exp(alpha_int_fun(t)))^2)
  }
  y2 <- integrand2(time_grid)
  
  integrand3 <- function(t){
    return(vol_fun(t) * dvol_dtheta3(t) / (exp(alpha_int_fun(t)))^2)
  }
  y3 <- integrand3(time_grid)
  
  # Value of dvar_dtheta3
  return(2 * exp_value^2 * (vol_value * dalpha_dtheta3_int_fun(t1) - trapz(time_grid, y2) + trapz(time_grid, y3)))
}

Fisher_matrix_der_exact_complete_nonuni <- function(time_grid, gm_paths, subdiv,
                                                    slope_fun, trend_fun, vol_fun,
                                                    dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                                    dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                                    dvol_dtheta1, dvol_dtheta2, dvol_dtheta3){
  # Convert various input formats to a matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    # Handle time-series objects
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Handle list of objects (can be ts, vectors, or matrices)
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # Handle single vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Handle matrices or data.frames directly
    gm_paths_matrix <- as.matrix(gm_paths)
  }
  
  output <- list(
    dl_dtheta1_total = 0,
    dl_dtheta2_total = 0,
    dl_dtheta3_total = 0,
    Fisher_matrix = matrix(0, nrow = 3, ncol = 3)
  )
  
  n <- nrow(gm_paths_matrix)
  num_paths <- ncol(gm_paths_matrix)
  # Write header to file
  # write("t0\tt1\tx0\tx1\tdlogL_dtheta1\tdlogL_dtheta2\tdlogL_dtheta3", file = "Fisher_matrix_der_exact.txt")
  
  # Initialize Fisher matrix accumulator BEFORE loops
  fisher_sum <- matrix(0, nrow = 3, ncol = 3)
  
  # Loop over each path
  for (path_num in 1:num_paths) {
    gm_path <- gm_paths_matrix[, path_num]
    dl_dtheta1 <- 0 # derivative totals end up being from the last path only (as coded)
    dl_dtheta2 <- 0
    dl_dtheta3 <- 0
    
    # Loop over each transition within the path
    for (i in 2:n) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x1 <- gm_path[i]
      x0 <- gm_path[i - 1]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun, subdiv)
      exp_value   <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun, alpha_int_fun, t0, t1, subdiv)
      vol_value   <- vol_int_dt(vol_fun,   alpha_int_fun, t0, t1, subdiv)
      mean_val    <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      var_val     <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      # Log-likelihood derivatives for the current transition
      dmean1 <- dmean_dtheta1(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta1, dtrend_dtheta1, t0, t1, x0, subdiv)
      dvar1  <- dvar_dtheta1(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta1, dvol_dtheta1, t0, t1, subdiv)
      dl1 <- 0.5 * dvar1 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean1 / var_val
      dl_dtheta1 <- dl_dtheta1 + dl1
      
      dmean2 <- dmean_dtheta2(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta2, dtrend_dtheta2, t0, t1, x0, subdiv)
      dvar2  <- dvar_dtheta2(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta2, dvol_dtheta2, t0, t1, subdiv)
      dl2 <- 0.5 * dvar2 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean2 / var_val
      dl_dtheta2 <- dl_dtheta2 + dl2
      
      dmean3 <- dmean_dtheta3(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta3, dtrend_dtheta3, t0, t1, x0, subdiv)
      dvar3  <- dvar_dtheta3(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta3, dvol_dtheta3, t0, t1, subdiv)
      dl3 <- 0.5 * dvar3 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean3 / var_val
      dl_dtheta3 <- dl_dtheta3 + dl3
      
      # Gradient vector for the current transition
      grad <- c(dl1, dl2, dl3)
      
      # Accumulate the outer product
      fisher_sum <- fisher_sum + tcrossprod(grad)
    }
  }
  
  # Total number of transitions: (n-1) per path
  total_transitions <- (n - 1) * num_paths
  
  # Normalize ONCE at the end
  fisher_final <- fisher_sum / total_transitions
  
  output$dl_dtheta1_total <- dl_dtheta1
  output$dl_dtheta2_total <- dl_dtheta2
  output$dl_dtheta3_total <- dl_dtheta3
  
  # Average over all paths
  output$Fisher_matrix <- fisher_final
  
  return(output)
}

Fisher_matrix_der_exact_complete <- function(T, gm_paths, subdiv,
                                             slope_fun, trend_fun, vol_fun,
                                             dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                             dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                             dvol_dtheta1, dvol_dtheta2, dvol_dtheta3){
  # Convert various input formats to a matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    # Handle time-series objects
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Handle list of objects (can be ts, vectors, or matrices)
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # Handle single vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Handle matrices or data.frames directly
    gm_paths_matrix <- as.matrix(gm_paths)
  }
  
  output <- list(
    dl_dtheta1_total = 0,
    dl_dtheta2_total = 0,
    dl_dtheta3_total = 0,
    Fisher_matrix = matrix(0, nrow = 3, ncol = 3)
  )
  
  n <- nrow(gm_paths_matrix)
  num_paths <- ncol(gm_paths_matrix)
  time_grid <- seq(0, T, length.out = n)
  
  # Initialize Fisher matrix accumulator BEFORE loops
  fisher_sum <- matrix(0, nrow = 3, ncol = 3)
  
  # Loop over each path
  for (path_num in 1:num_paths) {
    gm_path <- gm_paths_matrix[, path_num]
    dl_dtheta1 <- 0 # derivative totals end up being from the last path only (as coded)
    dl_dtheta2 <- 0
    dl_dtheta3 <- 0
    
    # Loop over each transition within the path
    for (i in 2:n) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x1 <- gm_path[i]
      x0 <- gm_path[i - 1]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun, subdiv)
      exp_value   <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun, alpha_int_fun, t0, t1, subdiv)
      vol_value   <- vol_int_dt(vol_fun,   alpha_int_fun, t0, t1, subdiv)
      mean_val    <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      var_val     <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      # Log-likelihood derivatives for the current transition
      dmean1 <- dmean_dtheta1(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta1, dtrend_dtheta1, t0, t1, x0, subdiv)
      dvar1  <- dvar_dtheta1(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta1, dvol_dtheta1, t0, t1, subdiv)
      dl1 <- 0.5 * dvar1 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean1 / var_val
      dl_dtheta1 <- dl_dtheta1 + dl1
      
      dmean2 <- dmean_dtheta2(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta2, dtrend_dtheta2, t0, t1, x0, subdiv)
      dvar2  <- dvar_dtheta2(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta2, dvol_dtheta2, t0, t1, subdiv)
      dl2 <- 0.5 * dvar2 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean2 / var_val
      dl_dtheta2 <- dl_dtheta2 + dl2
      
      dmean3 <- dmean_dtheta3(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta3, dtrend_dtheta3, t0, t1, x0, subdiv)
      dvar3  <- dvar_dtheta3(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta3, dvol_dtheta3, t0, t1, subdiv)
      dl3 <- 0.5 * dvar3 * (-1/var_val + (x1 - mean_val)^2/(var_val^2)) + (x1 - mean_val) * dmean3 / var_val
      dl_dtheta3 <- dl_dtheta3 + dl3
      
      # Gradient vector for the current transition
      grad <- c(dl1, dl2, dl3)
      
      # Accumulate the outer product
      fisher_sum <- fisher_sum + tcrossprod(grad)
    }
  }
  
  # Total number of transitions: (n-1) per path
  total_transitions <- (n - 1) * num_paths
  
  # Normalize ONCE at the end
  fisher_final <- fisher_sum / total_transitions
  
  output$dl_dtheta1_total <- dl_dtheta1
  output$dl_dtheta2_total <- dl_dtheta2
  output$dl_dtheta3_total <- dl_dtheta3
  
  # Average over all paths
  output$Fisher_matrix <- fisher_final
  
  return(output)
}

# Wrapper function for nleqslv
gradient_wrapper <- function(theta, T, X, subdivision,
                             slope_fun, trend_fun, vol_fun,
                             dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                             dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                             dvol_dtheta1, dvol_dtheta2, dvol_dtheta3) {
  
  params <- list(alpha = theta[1], beta = theta[2], sigma = theta[3])
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t   <- function(t) do.call(vol_fun,   c(list(t = t), params))
  
  # Derivatives
  dslope_dtheta1_t <- function(t) do.call(dslope_dtheta1, c(list(t = t), params))
  dslope_dtheta2_t <- function(t) do.call(dslope_dtheta2, c(list(t = t), params))
  dslope_dtheta3_t <- function(t) do.call(dslope_dtheta3, c(list(t = t), params))
  
  dtrend_dtheta1_t <- function(t) do.call(dtrend_dtheta1, c(list(t = t), params))
  dtrend_dtheta2_t <- function(t) do.call(dtrend_dtheta2, c(list(t = t), params))
  dtrend_dtheta3_t <- function(t) do.call(dtrend_dtheta3, c(list(t = t), params))
  
  dvol_dtheta1_t <- function(t) do.call(dvol_dtheta1, c(list(t = t), params))
  dvol_dtheta2_t <- function(t) do.call(dvol_dtheta2, c(list(t = t), params))
  dvol_dtheta3_t <- function(t) do.call(dvol_dtheta3, c(list(t = t), params))
  
  result <- Fisher_matrix_der_exact_complete_nonuni(
    T, X, subdivision,
    slope_fun = slope_fun_t,
    trend_fun = trend_fun_t,
    vol_fun   = vol_fun_t,
    dslope_dtheta1 = dslope_dtheta1_t,
    dslope_dtheta2 = dslope_dtheta2_t,
    dslope_dtheta3 = dslope_dtheta3_t,
    dtrend_dtheta1 = dtrend_dtheta1_t,
    dtrend_dtheta2 = dtrend_dtheta2_t,
    dtrend_dtheta3 = dtrend_dtheta3_t,
    dvol_dtheta1   = dvol_dtheta1_t,
    dvol_dtheta2   = dvol_dtheta2_t,
    dvol_dtheta3   = dvol_dtheta3_t
  )
  
  # Return the gradient as a vector
  return(c(result$dl_dtheta1_total, result$dl_dtheta2_total, result$dl_dtheta3_total))
}

gradient_wrapper_nonuni <- function(theta, time_grid, X, subdivision,
                                    slope_fun, trend_fun, vol_fun,
                                    dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                    dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                    dvol_dtheta1, dvol_dtheta2, dvol_dtheta3) {
  
  params <- list(alpha = theta[1], beta = theta[2], sigma = theta[3])
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t   <- function(t) do.call(vol_fun,   c(list(t = t), params))
  
  # Derivatives
  dslope_dtheta1_t <- function(t) do.call(dslope_dtheta1, c(list(t = t), params))
  dslope_dtheta2_t <- function(t) do.call(dslope_dtheta2, c(list(t = t), params))
  dslope_dtheta3_t <- function(t) do.call(dslope_dtheta3, c(list(t = t), params))
  
  dtrend_dtheta1_t <- function(t) do.call(dtrend_dtheta1, c(list(t = t), params))
  dtrend_dtheta2_t <- function(t) do.call(dtrend_dtheta2, c(list(t = t), params))
  dtrend_dtheta3_t <- function(t) do.call(dtrend_dtheta3, c(list(t = t), params))
  
  dvol_dtheta1_t <- function(t) do.call(dvol_dtheta1, c(list(t = t), params))
  dvol_dtheta2_t <- function(t) do.call(dvol_dtheta2, c(list(t = t), params))
  dvol_dtheta3_t <- function(t) do.call(dvol_dtheta3, c(list(t = t), params))
  
  result <- Fisher_matrix_der_exact_complete_nonuni(
    time_grid, X, subdivision,
    slope_fun = slope_fun_t,
    trend_fun = trend_fun_t,
    vol_fun   = vol_fun_t,
    dslope_dtheta1 = dslope_dtheta1_t,
    dslope_dtheta2 = dslope_dtheta2_t,
    dslope_dtheta3 = dslope_dtheta3_t,
    dtrend_dtheta1 = dtrend_dtheta1_t,
    dtrend_dtheta2 = dtrend_dtheta2_t,
    dtrend_dtheta3 = dtrend_dtheta3_t,
    dvol_dtheta1   = dvol_dtheta1_t,
    dvol_dtheta2   = dvol_dtheta2_t,
    dvol_dtheta3   = dvol_dtheta3_t
  )
  
  # Return the gradient as a vector
  return(c(result$dl_dtheta1_total, result$dl_dtheta2_total, result$dl_dtheta3_total))
}

MLE_grad_loglikelihood <- function(initial_guess, T, X, subdivision,
                                   slope_fun, trend_fun, vol_fun,
                                   dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                   dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                   dvol_dtheta1, dvol_dtheta2, dvol_dtheta3,
                                   method = "Newton", global = "cline"){
  sol <- nleqslv(
    initial_guess,
    function(theta)
      gradient_wrapper(theta, T, X, subdivision,
                       slope_fun, trend_fun, vol_fun,
                       dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                       dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                       dvol_dtheta1, dvol_dtheta2, dvol_dtheta3),
    method = "Newton", global = "cline"
  )
  
  # ## Methods and global strategies to try
  # methods <- c("Newton", "Broyden")
  # globals <- c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none")
  #
  # # List to store results
  # results <- list()
  #
  # # Loop over all combinations
  # cat("nleqslv results by method/global combination:\n")
  # for (m in methods) {
  #   for (g in globals) {
  #     res <- tryCatch({
  #       sol <- nleqslv(initial_guess, function(theta) gradient_wrapper(theta, T, X, subdivision,
  #         slope_fun, trend_fun, vol_fun,
  #         dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
  #         dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
  #         dvol_dtheta1, dvol_dtheta2, dvol_dtheta3),
  #         method = "Newton", global = "cline")
  #       list(method = m, global = g, x = sol$x, conv = sol$termcd)
  #     }, error = function(e) {
  #       list(method = m, global = g, x = NA, conv = NA)
  #     })
  #
  #     results[[paste(m, g, sep = "_")]] <- res
  #
  #     cat(sprintf("Method: %-7s  Global: %-7s  Solution: (%.4f, %.4f, %.4f)  Convergence: %d\n",
  #                 m, g,
  #                 ifelse(is.na(res$x[1]), NA, res$x[1]),
  #                 ifelse(is.na(res$x[2]), NA, res$x[2]),
  #                 ifelse(is.na(res$x[3]), NA, res$x[3]),
  #                 res$conv))
  #   }
  # }
  
  return(sol)
}

MLE_grad_loglikelihood_nonuni <- function(initial_guess, time_grid, X, subdivision,
                                          slope_fun, trend_fun, vol_fun,
                                          dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                          dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                                          dvol_dtheta1, dvol_dtheta2, dvol_dtheta3,
                                          method = "Newton", global = "cline"){
  sol <- nleqslv(
    initial_guess,
    function(theta)
      gradient_wrapper_nonuni(theta, time_grid, X, subdivision,
                              slope_fun, trend_fun, vol_fun,
                              dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                              dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3,
                              dvol_dtheta1, dvol_dtheta2, dvol_dtheta3),
    method = "Newton", global = "cline"
  )
  return(sol)
}

log_likelihood4_nonuni <- function(alpha, beta, sigma, time_grid, gm_paths, subdiv,
                                   slope_fun, trend_fun, vol_fun){
  
  # cat("Class of X:", class(X), "Dimensions:", dim(X), "\n")
  # Convert various input formats to a matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    # Handle time-series objects
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Handle list of objects (can be ts, vectors, or matrices)
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # Handle single vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Handle matrices or data.frames directly
    gm_paths_matrix <- as.matrix(gm_paths)
  }
  
  # Write header to file
  # write("x1\tt1\tx0\texp\ttrend\tvol\tmean\tvariance\tdnorm", file = "likelihood_mio4.txt")
  
  n <- nrow(gm_paths_matrix)
  num_paths <- ncol(gm_paths_matrix)
  
  params <- list(alpha = alpha, beta = beta, sigma = sigma)
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t   <- function(t) do.call(vol_fun,   c(list(t = t), params))
  
  total_logL <- 0
  
  # Process each path (column in the matrix)
  for (path_num in 1:num_paths) {
    gm_path <- gm_paths_matrix[, path_num]
    path_logL <- 0
    
    for (i in 2:n) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x1 <- gm_path[i]
      x0 <- gm_path[i - 1]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t, subdiv)
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1, subdiv)
      vol_value   <- vol_int_dt(vol_fun_t,   alpha_int_fun, t0, t1, subdiv)
      mean_val    <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      var_val     <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      value <- dnorm(x1, mean_val, sqrt(var_val), log = TRUE)
      path_logL <- path_logL + value
    }
    
    total_logL <- total_logL + path_logL
  }
  
  # Return negative sum (for minimization)
  return(-total_logL)
}

lk_estimates_nonuni <- function(subdivision, X, time_grid, theta_test,
                                slope_fun, trend_fun, vol_fun){
  opt_result <- mle(
    function(alpha, beta, sigma)
      log_likelihood4_nonuni(alpha, beta, sigma, time_grid, X, subdivision, slope_fun, trend_fun, vol_fun),
    start  = list(alpha = theta_test[1], beta = theta_test[2], sigma = theta_test[3]),
    method = "L-BFGS-B",
    lower  = c(-Inf, -Inf, 1e-6),
    upper  = c(Inf,  Inf,  Inf)
  )
  return(opt_result)
}

simulate_gm_process <- function(n_paths, n_steps, delta, x0 = 1,
                                slope_fun_t, trend_fun_t, vol_fun_t,
                                alpha_gm, beta_gm, sigma_gm, subdiv) {
  # n_steps are the number of steps
  T <- delta * n_steps
  time_grid <- seq(0, T, length.out = n_steps + 1)
  paths <- matrix(NA, nrow = n_steps + 1, ncol = n_paths)
  paths_exact <- matrix(NA, nrow = n_steps + 1, ncol = n_paths)
  paths[1, ] <- x0
  paths_exact[1, ] <- x0
  
  # Error counters
  error_log <- list(
    negative_variance = 0,
    na_mean = 0,
    infinite_mean = 0,
    na_sd = 0,
    infinite_sd = 0,
    error_mean = 0, # average of all mean errors across each path
    error_var  = 0,
    error_path = 0
  )
  
  for (j in 1:n_paths) {
    error_mean_path <- numeric(n_steps - 1)
    error_var_path  <- numeric(n_steps - 1)
    
    for (i in 2:(n_steps + 1)) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      
      x_prev <- paths[i - 1, j]
      x_prev_exact <- paths_exact[i - 1, j]
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t, subdiv)
      exp_value   <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1, subdiv)
      vol_value   <- vol_int_dt(vol_fun_t,   alpha_int_fun, t0, t1, subdiv)
      mean_gm <- gm_mean_dt(x_prev, t1, alpha_int_fun, trend_value)
      var_gm  <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      # NOTE: remaining exact mean/var + error tracking is unchanged in this abbreviated section
      paths[i, j]       <- rnorm(1, mean = mean_gm, sd = sqrt(max(var_gm, 0)))
      paths_exact[i, j] <- NA_real_
    }
  }
  
  return(list(time = time_grid, paths = paths, paths_exact = paths_exact, errors = error_log))
}