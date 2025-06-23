library(pracma)
library(sde)

# --- Auxiliary function for direct integration with trapz ---
integrate_direct_trapz <- function(func_to_integrate, lower, upper, subdivisions, ...) {
  if (abs(lower - upper) < .Machine$double.eps^0.75) {
    return(0)
  }
  if (subdivisions < 2) {
    stop("integrate_direct_trapz requires at least 2 subdivisions (points).")
  }
  grid <- seq(lower, upper, length.out = subdivisions)
  
  func_name_str <- deparse(substitute(func_to_integrate))
  # If func_to_integrate is one of those that expect a scalar as the first argument of the grid.
  if (grepl("integrand_for_trend_term_optimized", func_name_str) || 
      grepl("integrand_for_vol_term_optimized", func_name_str) ) { 
    values <- sapply(grid, func_to_integrate, ...)
  } else { # For functions like slope, trend, vol that are assumed to be vectorized
    values <- func_to_integrate(grid, ...)
  }
  
  return(pracma::trapz(grid, values))
}

# --- Functions for OPTIMIZED integrands (use spline_F_slope) ---
integrand_for_trend_term_optimized <- function(u_mesh_val, current_t2, 
                                               spline_F_slope_func, # Changed: spline function
                                               trend_func) {
  # integral_slope_u_t2 is obtained from the precalculated spline
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  return(trend_func(u_mesh_val) * exp(integral_slope_u_t2))
}

integrand_for_vol_term_optimized <- function(u_mesh_val, current_t2, 
                                             spline_F_slope_func, 
                                             vol_func) {
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  term_in_paren <- vol_func(u_mesh_val) * exp(integral_slope_u_t2)
  return(term_in_paren^2)
}

# --- Kernel with OPTIMIZED integral calculation ---
boundary_kernel_optimized <- function(c1, c2, t1, x1, t2, x2, 
                                      spline_F_slope_func,
                                      trend_f, vol_f, 
                                      discount_val, subdivisions,
                                      # No longer needs slope_f directly here for integrals
                                      # But yes for logging if slope_f(t1) is desired
                                      slope_f_for_log_only # For logging only
) {
  
  # 1. exp_int_slope_t1_t2 = exp(integral_{t1}^{t2} slope(u)du)
  # Obtained from the precalculated spline
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
  #              slope_f_for_log_only(t1), trend_f(t1), vol_f(t1),
  #              marginal_mean_val, current_marginal_var, sep = "\t")
  # try(write(line, file = "stopping_boundary_Abel.txt", append = TRUE, ncolumns = 11), silent=TRUE)
  
  # --- Kernel value calculation ---
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


# --- Main boundary function ---
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0,
                      slope, trend, vol, errors = FALSE,
                      trapz_subdivs = 100, 
                      cumtrapz_subdivs_slope = 1000 # For the slope spline
) {
  
  # --- Pre-calculation of slope integral using cumtrapz and spline ---
  # Ensure that min(time_line) is time_line[1] if sorted
  min_time <- min(time_line) 
  max_time <- max(time_line)
  
  # Ensure that cumtrapz_subdivs_slope is at least 2
  safe_cumtrapz_subdivs_slope <- max(2, cumtrapz_subdivs_slope)
  fine_grid_for_slope_spline <- seq(min_time, max_time, length.out = safe_cumtrapz_subdivs_slope)
  
  slope_vals_on_fine_grid <- slope(fine_grid_for_slope_spline)
  # Ensure that slope_vals_on_fine_grid is a vector of the correct size
  if(!is.vector(slope_vals_on_fine_grid) || length(slope_vals_on_fine_grid) != length(fine_grid_for_slope_spline)) {
    # If slope is not vectorized, apply sapply.
    # This is important if slope is not of the form rep(const, length(t))
    slope_vals_on_fine_grid <- sapply(fine_grid_for_slope_spline, slope)
  }
  
  # F(x) = integral_{fine_grid_for_slope_spline[1]}^{x} slope(u) du
  cumulative_integral_slope_values <- pracma::cumtrapz(fine_grid_for_slope_spline, slope_vals_on_fine_grid)
  
  # Create spline function for F(x)
  spline_F_slope <- splinefun(fine_grid_for_slope_spline, cumulative_integral_slope_values, method = "natural")
  # --- End pre-calculation ---
  
  # Write header to file (if logging is active in the kernel)
  # write("c1\tc2\tt1\tx1\tt2\tx2\tslope(t1)\ttrend(t1)\tvol(t1)\tmarginal_mean\tmarginal_var", file = "stopping_boundary_Abel.txt")
  
  N <- length(time_line)
  expiration <- time_line[N] # Assumes the last point is expiration
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
    
    print(paste("Picard Iteration:", j)) # To see progress
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
      warning(paste("Error became NA/NaN/Inf in iteration",j,". Stopping. Using last valid boundary."))
      bnd <- last_valid_bnd 
      if (errors) er <- c(er, e) 
      break 
    }
    if(j > 1 || (j==1 && !(is.na(e) || is.infinite(e) || is.nan(e)))) { # Ensure that bnd_old is 'good'
      last_valid_bnd <- bnd_old 
    }
    
    
    if (errors) er <- c(er, e)
    
    # Reduce iteration limit if it is suspected to be very slow initially
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
  
  time_line <- seq(0, expiration, length.out = partition_length + 1)
  
  params <- list(alpha = theta[1], beta = theta[2], sigma = theta[3])
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t <- function(t) do.call(vol_fun, c(list(t = t), params))
  
  bnd <- boundary(tol = 1e-3 ,strike = strike, time_line = time_line,
                  discount = discount, slope = slope_fun_t, trend = trend_fun_t,
                  vol = vol_fun_t, errors = TRUE, trapz_subdivs = 100, cumtrapz_subdivs_slope = 1000)
  
  return(bnd$boundary)
}

simulate_gm_process <- function(n_paths, n_steps, delta, x0, slope_fun_t, trend_fun_t, vol_fun_t, alpha_gm, beta_gm, sigma_gm, subdiv) {
  # n_steps, are the steps taken 
  T <- delta*n_steps
  time_grid <- seq(0, T, length.out = n_steps+1)
  paths <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths_exact <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths[1, ] <- x0
  paths_exact[1, ] <- x0
  
  # Error counters
  error_log <- list(
    negative_variance = 0,
    na_mean = 0,
    infinite_mean = 0,
    na_sd = 0,
    infinite_sd = 0,
    error_mean = 0, # mean of all mean errors in each path
    error_var = 0,
    error_path = 0
  )
  
  for (j in 1:n_paths) {
    error_mean_path <- numeric(n_steps-1)
    error_var_path <- numeric(n_steps-1)
    for (i in 2:(n_steps+1)) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x_ant <- paths[i - 1,j]
      x_ant_exact <- paths_exact[i - 1,j]
      
      # OUB
      # mean_exact <- (x_ant + beta_gm*(1 - exp(-alpha_gm*(t1-t0)))/alpha_gm)*sinh(alpha_gm*(T-t1))/sinh(alpha_gm*(T-t0)) +
      #   (xf - beta_gm*(-1 + exp(alpha_gm*(T-t1)))*alpha_gm)*sinh(alpha_gm*(t1-t0))/sinh(alpha_gm*(T-t0))
      # 
      # var_exact <- sigma_gm*sinh(alpha_gm*(T-t1))*sinh(alpha_gm*(t1-t0))/(alpha_gm*sinh(alpha_gm*(T-t0)))
      
      # # Hull-White model periodic mean difficult
      # mean_exact <- (x_ant + alpha_gm*(exp(t1 - t0) - 1) + beta_gm*(exp(t1 - t0)*(sin(sigma_gm*t1) - sigma_gm*cos(sigma_gm*t1)) - (sin(sigma_gm*t0) - sigma_gm*cos(sigma_gm*t0)))/(1 + sigma_gm^2))/exp(t1 - t0)
      # 
      # var_exact <- 0.5*(1 - exp(-2*(t1 - t0)))
      
      # Hull-White model periodic mean version 1
      # mean_exact <-exp(alpha_gm*(t1 - t0))*(x_ant_exact + beta_gm*(1 - exp(-alpha_gm*(t1 - t0)))/alpha_gm + ((alpha_gm*sin(sigma_gm*t0) + sigma_gm*cos(sigma_gm*t0)) - exp(alpha_gm*(t1 - t0))*(alpha_gm*sin(sigma_gm*t1) + sigma_gm*cos(sigma_gm*t1)))/(alpha_gm^2 + sigma_gm^2))
      # 
      # var_exact <- (exp(2*alpha_gm*(t1 - t0)) - 1)/(2*alpha_gm)
      
      # # Hull-White model periodic mean version 2
      # mean_exact <-exp(alpha_gm*(t1 - t0))*(x_ant_exact + beta_gm*(1 - exp(-alpha_gm*(t1 - t0)))/alpha_gm + sigma_gm*((alpha_gm*sin(t0) + cos(t0)) - exp(alpha_gm*(t1 - t0))*(alpha_gm*sin(t1) + cos(t1)))/(alpha_gm^2 + 1))
      # 
      # var_exact <- (exp(2*alpha_gm*(t1 - t0)) - 1)/(2*alpha_gm)
      # # OU stationary
      # # mean_exact <- -beta_gm/alpha_gm + (x_ant_exact + beta_gm/alpha_gm)*exp(alpha_gm*(t1 - t0))
      # # 
      # # var_exact <- -sigma_gm^2*(1 - exp(2*alpha_gm*(t1 - t0)))/(2*alpha_gm)
      # 
      # # Hull White model root squared
      # mean_exact <- exp(alpha_gm*(t1 - t0))*(x_ant_exact + exp(alpha_gm*t0)*(beta_gm*(exp(-alpha_gm*t0) - exp(-alpha_gm*t1))/alpha_gm + sigma_gm*(sqrt(pi)*(erf(sqrt(alpha_gm*t1)) - erf(sqrt(alpha_gm*t0))) - sqrt(t1)*exp(-alpha_gm*t1) + sqrt(t0)*exp(-alpha_gm*t0))/(2*sqrt(alpha_gm))/alpha_gm))
      #
      # var_exact <- (exp(2*alpha_gm*(t1 - t0)) - 1)/(2*alpha_gm)
      
      # # decaying volatility
      
      # mean_exact <- exp(-(t1 - t0))*(x_ant_exact + alpha_gm*(exp(t1 - t0) - 1))
      # 
      # if (2 * sigma_gm^2 - 4 * sigma_gm + 4 == 0) {
      #   warning("Denominator for Integral 2 is zero. Result may be infinite or undefined.")
      # }
      # if (2 * (1 - sigma_gm) == 0) {
      #   warning("Denominator for Integral 3a is zero. Result may be infinite or undefined.")
      # }
      # if (8 * sigma_gm^2 - 8 * sigma_gm + 4 == 0) {
      #   warning("Denominator for Integral 3b is zero. Result may be infinite or undefined.")
      # }
      # 
      # # --- Part 1: First Integral ---
      # # Integral of theta_2^2 * e^(2u) du
      # integral1 <- (beta_gm^2 / 2) * (exp(2 * t1) - exp(2 * t0))
      # 
      # # --- Part 2: Second Integral ---
      # # Integral of 2 * theta_2 * e^((2-theta_3)u) * sin(theta_3 u) du
      # denom_I2 <- (2 * sigma_gm^2 - 4 * sigma_gm + 4)
      # 
      # term_t1_I2 <- exp((2 - sigma_gm) * t1) *
      #   ((2 - sigma_gm) * sin(sigma_gm * t1) - sigma_gm * cos(sigma_gm * t1))
      # 
      # term_t0_I2 <- exp((2 - sigma_gm) * t0) *
      #   ((2 - sigma_gm) * sin(sigma_gm * t0) - sigma_gm * cos(sigma_gm * t0))
      # 
      # integral2 <- (2 * beta_gm / denom_I2) * (term_t1_I2 - term_t0_I2)
      # 
      # # --- Part 3: Third Integral ---
      # # Integral of e^((2-2*theta_3)u) * sin^2(theta_3 u) du
      # # This integral is split into two parts: 3a (exponential) and 3b (exponential-cosine)
      # 
      # # Part 3a: Integral of e^((2-2*theta_3)u) du
      # integral3a <- (1 / (2 * (1 - sigma_gm))) * (exp((2 - 2 * sigma_gm) * t1) - exp((2 - 2 * sigma_gm) * t0))
      # 
      # # Part 3b: Integral of e^((2-2*theta_3)u) * cos(2*theta_3 u) du
      # denom_I3b <- (8 * sigma_gm^2 - 8 * sigma_gm + 4)
      # 
      # term_t1_I3b <- exp((2 - 2 * sigma_gm) * t1) *
      #   ((2 - 2 * sigma_gm) * cos(2 * sigma_gm * t1) + 2 * sigma_gm * sin(2 * sigma_gm * t1))
      # 
      # term_t0_I3b <- exp((2 - 2 * sigma_gm) * t0) *
      #   ((2 - 2 * sigma_gm) * cos(2 * sigma_gm * t0) + 2 * sigma_gm * sin(2 * sigma_gm * t0))
      # 
      # integral3b <- (1 / denom_I3b) * (term_t1_I3b - term_t0_I3b)
      # 
      # # Combine parts of Integral 3
      # integral3 <- (1/2) * (integral3a - integral3b)
      # 
      # # --- Final Result ---
      # # Multiply by the leading e^(-2*t1) term
      # var_exact <- exp(-2 * t1) * (integral1 + integral2 + integral3)
      
      # OU time-dependent volatility difficult
      # mean_exact <- exp(-(t1 - t0))*(x_ant_exact + alpha_gm*(exp(t1 - t0) - 1))
      #
      # var_exact <- beta_gm^2*exp(-2*t1)*(exp(2*t1*(1 - sigma_gm)) - exp(2*t0*(1 - sigma_gm)))/(2*(1 - sigma_gm))
      
      # Hull White model logarithmic mean (very complex) and rational volatility (very complex)
      mean_exact <- 0
      var_exact <- 1
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t,subdiv) 
      exp_value <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1,subdiv)
      vol_value <- vol_int_dt(vol_fun_t, alpha_int_fun, t0, t1, subdiv)
      mean_gm <- gm_mean_dt(x_ant, t1, alpha_int_fun, trend_value)
      var_gm <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      error_mean_path[i] <- mean_gm - mean_exact
      error_var_path[i] <- var_gm - var_exact
      
      # Parameter validation for rnorm()
      if (is.na(mean_gm)) {
        cat("NaN in mean, in time step ",i, "\r")
        error_log$na_mean <- error_log$na_mean + 1
        next
      }
      if (is.infinite(mean_gm)) {
        error_log$infinite_mean <- error_log$infinite_mean + 1
        next
      }
      if (is.na(var_gm)) {
        error_log$na_sd <- error_log$na_sd + 1
        next
      }
      if (var_gm < 0) {
        error_log$negative_variance <- error_log$negative_variance + 1
        var_gm <- 0 # Force to 0 to avoid NaN
      }
      if (is.infinite(sqrt(var_gm))) {
        error_log$infinite_sd <- error_log$infinite_sd + 1
        next
      }
      
      paths[i, j] <- rnorm(1, mean = mean_gm, sd = sqrt(var_gm))
      paths_exact[i, j] <- rnorm(1, mean = mean_exact, sd = sqrt(var_exact))
    }
    error_log$error_mean <- mean(error_mean_path)
    error_log$error_var <- mean(error_var_path)
    error_log$error_path <- mean(paths_exact - paths)
  }
  return(list(time = time_grid, paths = paths, paths_exact = paths_exact, errors = error_log))
}


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
  vol_fun_t   <- function(t) do.call(vol_fun, c(list(t = t), params_true))
  
  all_paths <- simulate_gm_process(n_paths = M, n_steps = N, delta = delta, x0 = x0,
                                   slope_fun_t, trend_fun_t, vol_fun_t, theta_true[1], theta_true[1],theta_true[1], subdivision)$paths
  
  # --- PHASE 2: TRAINING (INDIVIDUAL ESTIMATION AND FIM) ---
  if (verbose) cat("--- PHASE 2: Estimating parameters and FIM for each of the", M, "trajectories ---\n")
  
  all_estimates <- vector("list", M)
  all_fims <- vector("list", M)
  training_duration <- training_points * delta
  
  for (i in 1:M) {
    if (verbose && i %% 20 == 0) cat("  Processing training trajectory:", i, "/", M, "\n")
    
    X_training <- all_paths[1:(training_points + 1), i, drop = FALSE]
    
    result <- tryCatch(
      lk_estimates(subdivision, X_training, training_duration, theta_test = theta_test,
                   slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun),
      error = function(e) { message("Error in estimation (trajectory ", i, "): ", e$message); NULL }
    )
    
    if (is.null(result)) next
    
    theta_estimated <- as.numeric(result@coef)
    all_estimates[[i]] <- theta_estimated
    
    params_est <- list(alpha = theta_estimated[1], beta = theta_estimated[2], sigma = theta_estimated[3])
    
    # Calculate FIM for this trajectory
    fim_i <- Fisher_matrix_der_exact_complete(
      training_duration, X_training, subdivision,
      slope_fun = function(t) do.call(slope_fun, c(list(t = t), params_est)),
      trend_fun = function(t) do.call(trend_fun, c(list(t = t), params_est)),
      vol_fun = function(t) do.call(vol_fun, c(list(t = t), params_est)),
      dslope_dtheta1 = function(t) do.call(dslope_dtheta1, c(list(t = t), params_est)), dslope_dtheta2 = function(t) do.call(dslope_dtheta2, c(list(t = t), params_est)), dslope_dtheta3 = function(t) do.call(dslope_dtheta3, c(list(t = t), params_est)),
      dtrend_dtheta1 = function(t) do.call(dtrend_dtheta1, c(list(t = t), params_est)), dtrend_dtheta2 = function(t) do.call(dtrend_dtheta2, c(list(t = t), params_est)), dtrend_dtheta3 = function(t) do.call(dtrend_dtheta3, c(list(t = t), params_est)),
      dvol_dtheta1 = function(t) do.call(dvol_dtheta1, c(list(t = t), params_est)), dvol_dtheta2 = function(t) do.call(dvol_dtheta2, c(list(t = t), params_est)), dvol_dtheta3 = function(t) do.call(dvol_dtheta3, c(list(t = t), params_est))
    )$Fisher_matrix
    
    all_fims[[i]] <- fim_i
  }
  
  # --- PHASE 3: AGGREGATION (AVERAGING FIM) ---
  if (verbose) cat("--- PHASE 3: Averaging FIMs to obtain a robust matrix ---\n")
  
  valid_fims <- all_fims[!sapply(all_fims, is.null)]
  if(length(valid_fims) == 0) stop("Could not calculate any Fisher Matrix. All estimations failed.")
  
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
    boundary <- boundary_wrapper(theta = theta_estimated, strike = strike, partition_length = partition_length-1, discount = discount, expiration = expiration, slope_fun, trend_fun, vol_fun)
    h_alpha <- der_tol * abs(theta_estimated[1]); h_beta <- der_tol * abs(theta_estimated[2]); h_sigma <- der_tol * abs(theta_estimated[3])
    g1 <- (boundary_wrapper(c(theta_estimated[1] + h_alpha, theta_estimated[2:3]), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun) - boundary_wrapper(c(theta_estimated[1] - h_alpha, theta_estimated[2:3]), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun)) / (2*h_alpha)
    g2 <- (boundary_wrapper(c(theta_estimated[1], theta_estimated[2] + h_beta, theta_estimated[3]), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun) - boundary_wrapper(c(theta_estimated[1], theta_estimated[2] - h_beta, theta_estimated[3]), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun)) / (2*h_beta)
    g3 <- (boundary_wrapper(c(theta_estimated[1:2], theta_estimated[3] + h_sigma), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun) - boundary_wrapper(c(theta_estimated[1:2], theta_estimated[3] - h_sigma), strike, expiration, partition_length-1, discount, slope_fun, trend_fun, vol_fun)) / (2*h_sigma)
    
    grad_matrix <- as.matrix(data.frame(g1, g2, g3))
    
    variances <- rowSums((grad_matrix %*% fisher_info_inv_averaged) * grad_matrix)
    std_dev <- sqrt(variances)
    
    z_critical <- qnorm(1 - z_alpha / 2)
    upper_bound <- boundary + z_critical * std_dev
    lower_bound <- boundary - z_critical * std_dev
    
    inference_result <- list(lower_bound = lower_bound, boundary_est = boundary, upper_bound = upper_bound, est_theta = theta_estimated)
    all_results_list[[i]] <- inference_result
    
    # save the specified file index
    file_index <- i # Replace with the desired file index
    folder_name <- "rds_files"
    
    # Define file paths
    x_path_file <- file.path(folder_name, paste0("path_OU_", file_index, ".rds"))
    inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_OU_", file_index, ".rds"))
    
    # save the data from the .rds files
    saveRDS(all_paths[,i], x_path_file)
    saveRDS(inference_result, inferred_boundary_file)
    
    # # Save results
    # saveRDS(all_paths[, i, drop = FALSE], file.path(save_path, paste0("path_", i, ".rds")))
    # saveRDS(inference_result, file.path(save_path, paste0("inference_", i, ".rds")))
  }
  
  if (verbose) cat("--- Process completed. Results saved in '", save_path, "' ---\n")
  
  return(invisible(list(
    averaged_fim_inverse = fisher_info_inv_averaged,
    all_inference_results = all_results_list,
    all_paths = all_paths
  )))
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
# ## derivatives, when it has t it only depends on t and not on the parameter
# # slope derivatives:
# dslope_dtheta1 <- function(t, alpha,...) -coth(alpha*(T - t)) + alpha*(T - t)/sinh(alpha*(T - t))^2 
# dslope_dtheta2 <- function(t,...) rep(0, length(t)) 
# dslope_dtheta3 <- function(t,...) rep(0, length(t)) 
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t, alpha, beta,...) xf*csch(alpha*(T - t)) - alpha*xf*(T - t)*csch(alpha*(T - t))*coth(alpha*(T - t)) - beta*(T - t)/2*sech(alpha*(T - t)/2)^2
# dtrend_dtheta2 <- function(t, alpha,...) -tanh(alpha*(T - t)/2)
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # vol derivatives
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t, sigma,...) 2*sigma

#Hull White model periodic mean difficult
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
# ## derivatives, when it has t it only depends on t and not on the parameter
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(0, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(1, length(t))
# dtrend_dtheta2 <- function(t, sigma,...) sin(sigma*t)
# dtrend_dtheta3 <- function(t, beta, sigma,...) beta*t*cos(sigma*t)
# 
# # vol derivatives
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
# # derivatives, when it has t it only depends on t and not on the parameter
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(0, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(1, length(t))
# dtrend_dtheta2 <- function(t,...) rep(0, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # vol derivatives
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
# ## derivatives, when it has t it only depends on t and not on the parameter
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) sqrt(t)
# 
# # vol derivatives
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

## derivatives
# slope derivatives:
dslope_dtheta1 <- function(t,...) rep(1, length(t)) # when it has t it only depends on t and not on the parameter
dslope_dtheta2 <- function(t,...) rep(0, length(t))
dslope_dtheta3 <- function(t,...) rep(0, length(t))

# trend derivatives
dtrend_dtheta1 <- function(t,...) rep(0, length(t))
dtrend_dtheta2 <- function(t,...) rep(1, length(t))
dtrend_dtheta3 <- function(t,...) rep(0, length(t))

# vol derivatives
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
# ## derivatives
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t, sigma,...) t*cos(sigma*t)
# 
# # vol derivatives
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
# ## derivatives
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) 
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t, sigma,...) sin(t)
# 
# # vol derivatives
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
# ## derivatives
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # vol derivatives
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
# ## derivatives
# # slope derivatives:
# #dslope_dtheta1 <- function(t,...) rep(0, length(t)) # 1,2º
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) # 3º
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# #dtrend_dtheta1 <- function(t,...) rep(1, length(t)) # 2º
# dtrend_dtheta1 <- function(t,...) rep(0, length(t)) # 3º
# #dtrend_dtheta2 <- function(t, beta,...) t/(1 + beta*t) # 1º
# dtrend_dtheta2 <- function(t,...) log(1 + t)  # 2º
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # vol derivatives
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
# ## derivatives
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t))
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # trend derivatives
# dtrend_dtheta1 <- function(t,...) rep(0, length(t))
# dtrend_dtheta2 <- function(t,...) rep(1, length(t))
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
# 
# # vol derivatives
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t, sigma,...) -t/(1 + sigma*t)^2

# --- 2. Define all parameters for execution ---
params_ejecucion <- list(
  M = 2,
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
  der_tol = 1e-5,
  subdivision = 100,
  save_path = "resultados_hull_white_robust",
  verbose = TRUE
)

resultados <- do.call(infer_boundaries_robust, params_ejecucion)