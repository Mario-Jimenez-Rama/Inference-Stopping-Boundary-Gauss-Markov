library(pracma)
library(sde)

####
## - The function "boundary" computes the optimal stopping boundary of the (exp.
## discounted optimal stopping problem with finite horizon and gain function 
## G(x) = (strike - x)^+. The underlying process with parameters is a 
## time-dependent Ornstein-Uhlenbeck given by the SDE
##         dX_t = (slope(t)X_t + trend(t)) + vol(t)dW_t.
## - The boundary is computed by running a Picard iteration algorithm,
## that stops when the L2 distance between consecutive boundaries is less than
## tol, and solves the free-boundary equation.
## - The boundary is computed at the time points provided in time_line
## - If errors == TRUE, a vector of the errors of the Picard scheme is provided
## alongside the boundary
####
# --- Auxiliary function for direct integration with trapz ---
integrate_direct_trapz <- function(func_to_integrate, lower, upper, subdivisions, ...) {
  if (abs(lower - upper) < .Machine$double.eps^0.75) { # If numerically equal
    return(0)
  }
  # seq() behaves well if lower > upper, it creates a descending sequence.
  # trapz() should also handle it (or give an error if not ordered, depending on implementation).
  # For safety, if one wanted to force ascending order and sign:
  # if (lower > upper) {
  #   return(-integrate_direct_trapz(func_to_integrate, upper, lower, subdivisions, ...))
  # }
  if (subdivisions < 2) {
    stop("integrate_direct_trapz requires at least 2 subdivisions (points).")
  }
  grid <- seq(lower, upper, length.out = subdivisions)
  
  # Vectorize the call to func_to_integrate if it is not inherently vectorized
  # This is crucial for functions like integrand_for_trend_term_direct
  # which take additional arguments (...).
  # We will use sapply if func_to_integrate does not seem to be vectorized over its first argument.
  # A simple test: func_to_integrate(grid[1:2], ...)
  # However, it is safer to assume it must be vectorized by the user or use sapply.
  
  # In this case, the integrand_for_xxx_term_direct functions are designed
  # to take a scalar 'u_mesh_val' as the first argument. So sapply is necessary.
  
  # Determine if the function is one of our internal functions that needs sapply
  func_name_str <- deparse(substitute(func_to_integrate))
  if (grepl("integrand_for_trend_term_direct", func_name_str) || 
      grepl("integrand_for_vol_term_direct", func_name_str) ||
      !is.vector(try(func_to_integrate(grid, ...), silent = TRUE)) # Generic attempt to detect non-vectorization
  ) { 
    values <- sapply(grid, func_to_integrate, ...)
  } else { # For functions like slope, trend, vol that are assumed to be vectorized
    values <- func_to_integrate(grid, ...)
  }
  
  return(pracma::trapz(grid, values))
}

# --- Functions for the integrands of the mean and variance terms (direct calculation) ---
integrand_for_trend_term_direct <- function(u_mesh_val, current_t2, slope_func, trend_func, subdivs_inner) {
  integral_slope_u_t2 <- integrate_direct_trapz(slope_func, u_mesh_val, current_t2, subdivs_inner)
  return(trend_func(u_mesh_val) * exp(integral_slope_u_t2))
}

integrand_for_vol_term_direct <- function(u_mesh_val, current_t2, slope_func, vol_func, subdivs_inner) {
  integral_slope_u_t2 <- integrate_direct_trapz(slope_func, u_mesh_val, current_t2, subdivs_inner)
  term_in_paren <- vol_func(u_mesh_val) * exp(integral_slope_u_t2)
  return(term_in_paren^2)
}

# --- Kernel with direct integral calculation ---
boundary_kernel_direct <- function(c1, c2, t1, x1, t2, x2, 
                                   slope_f, trend_f, vol_f, 
                                   discount_val, subdivisions, # subdivisions for all trapz here
                                   strike_price_for_log) { # strike_price_for_log only for writing
  
  # 1. exp_int_slope_t1_t2 = exp(integral_{t1}^{t2} slope(u)du)
  integral_slope_t1_t2 <- integrate_direct_trapz(slope_f, t1, t2, subdivisions)
  exp_int_slope_t1_t2 <- exp(integral_slope_t1_t2)
  
  # 2. integral_trend_term = integral_{t1}^{t2} trend(u) * exp(integral_{u}^{t2} slope(v)dv) du
  integral_trend_term <- integrate_direct_trapz(
    integrand_for_trend_term_direct, t1, t2, subdivisions,
    current_t2 = t2, slope_func = slope_f, trend_func = trend_f, subdivs_inner = subdivisions
  )
  
  # 3. integral_vol_term (this is the variance directly)
  # integral_{t1}^{t2} (vol(u) * exp(integral_{u}^{t2} slope(v)dv))^2 du
  marginal_var_val_direct <- integrate_direct_trapz(
    integrand_for_vol_term_direct, t1, t2, subdivisions,
    current_t2 = t2, slope_func = slope_f, vol_func = vol_f, subdivs_inner = subdivisions
  )
  
  # Calculation of marginal mean and standard deviation
  marginal_mean_val <- x1 * exp_int_slope_t1_t2 + integral_trend_term
  
  current_marginal_var <- marginal_var_val_direct
  if (current_marginal_var < 0 && abs(current_marginal_var) < .Machine$double.eps*100) {
    current_marginal_var <- 0
  } else if (current_marginal_var < 0) {
    # warning(paste("Direct marginal_var_val < 0 :", current_marginal_var, "at t1,t2",t1,t2))
    current_marginal_var <- 0 # Force non-negativity
  }
  marginal_sd_val <- sqrt(current_marginal_var)
  
  # --- Logging (similar to before, now with direct t1, t2) ---
  # Ensure c1, c2, x1, x2 are scalars for this log section
  # (the function will be called with scalars in the K2 loop)
  line <- paste(c1, c2, t1, x1, t2, x2,
                slope_f(t1), trend_f(t1), vol_f(t1), # Assumes slope_f etc. can take a scalar
                marginal_mean_val, current_marginal_var, sep = "\t")
  try(write(line, file = "stopping_boundary_Abel.txt", append = TRUE, ncolumns = 11), silent=TRUE)
  
  # --- Kernel value calculation (as before) ---
  x2_std <- 0 
  if (marginal_sd_val > .Machine$double.eps^0.5) {
    x2_std <- (x2 - marginal_mean_val) / marginal_sd_val
  } else {
    diff_val_std <- x2 - marginal_mean_val
    x2_std <- ifelse(abs(diff_val_std) < .Machine$double.eps^0.5, 0, ifelse(diff_val_std > 0, Inf, -Inf))
  }
  
  normal_dist <- pnorm(x2_std, mean = 0, sd = 1, lower.tail = TRUE)
  normal_dens <- dnorm(x2_std, mean = 0, sd = 1)
  
  # If sd is effectively zero, the density term (c2 * sd * dens) will be zero.
  # No special adjustment for normal_dens is needed if marginal_sd_val is used below.
  
  term1 <- (c1 - c2 * marginal_mean_val) * normal_dist
  term2 <- c2 * marginal_sd_val * normal_dens # This term becomes zero if sd is very small
  
  time_diff <- t2 - t1 # time_diff = t2 - t1
  if (abs(time_diff) < .Machine$double.eps^0.75 && discount_val == 0) { # Avoid exp(0) * (Inf) if t1=t2 and discount=0
    # If t1=t2, the Kernel should be zero if there is no immediate payment at that instant.
    # The original Kernel formula has e^(-lambda(t2-t1)). If t1=t2, this is 1.
    # The marginal mean is x1. The marginal variance is 0. Sd is 0.
    # x2_std is (x2-x1)/0 -> Inf, -Inf, or 0.
    # If x2 > x1, dist=1, dens=0. K = (c1-c2*x1)*1 = c1-c2*x1
    # If x2 < x1, dist=0, dens=0. K = 0
    # If x2 = x1, dist=0.5, dens=dnorm(0). K = (c1-c2*x1)*0.5
    # This is for the special case t1=t2.
    # The problem definition generally assumes t2 > t1.
    # If t1=t2, the probability 1{Xt2 <= x2} is 1 if x1 <= x2, and 0 if x1 > x2.
    # This does not match the above. Better to trust that t2 > t1 for the kernel.
    # If they are equal, the kernel value must be defined by the problem context (often 0).
    # For now, the general formula is allowed to act.
  }
  
  K_val <- exp(-discount_val * (t2 - t1)) * (term1 + term2)
  
  return(K_val)
}


# --- Main boundary function ---
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0,
                      slope, trend, vol, errors = FALSE,
                      trapz_subdivs = 100 # Only subdivision needed now
                      # cumtrapz_subdivs_slope is no longer necessary
) {
  
  # Write header to file (updated for direct t1, x1, t2, x2)
  write("c1\tc2\tt1\tx1\tt2\tx2\tslope(t1)\ttrend(t1)\tvol(t1)\tmarginal_mean\tmarginal_var", file = "stopping_boundary_Abel.txt")
  
  N <- length(time_line)
  expiration <- time_line[N]
  delta <- time_line[2:N] - time_line[1:(N-1)] # For the K2 sum
  
  if (errors) er <- c()
  
  initial_bnd_val <- strike
  slope_at_expiration <- slope(expiration)
  trend_at_expiration <- trend(expiration)
  
  if (abs(discount - slope_at_expiration) > .Machine$double.eps^0.75) {
    initial_bnd_val <- min((trend_at_expiration + discount * strike) / (discount - slope_at_expiration), strike)
  } else {
    # warning("Denominator (discount - slope(expiration)) is close to zero in boundary init.")
  }
  bnd <- rep(initial_bnd_val, N)
  
  # No pre-calculations of I1, I2, I3
  
  e <- 1
  j <- 0
  
  while (e > tol) {
    j <- j + 1
    bnd_old <- bnd
    
    print(paste("Picard Iteration:", j)) # To see progress
    print(paste("error:", e)) # To see progress
    
    for (i in (N - 1):1) { # From t_{N-1} to t_1 (index 1 = time_line[1])
      t_current <- time_line[i]
      b_current_old <- bnd_old[i] # x1 for the kernel
      
      # K1: Kernel(strike, 1, t_current, b_current_old, expiration, strike)
      K1 <- boundary_kernel_direct(c1 = strike, c2 = 1,
                                   t1 = t_current, x1 = b_current_old,
                                   t2 = expiration, x2 = strike,
                                   slope_f = slope, trend_f = trend, vol_f = vol,
                                   discount_val = discount, subdivisions = trapz_subdivs,
                                   strike_price_for_log = strike) # last arg only for log if needed
      
      # K2: Sum of Kernels for future times
      sum_K2_delta_val <- 0
      if (i < N) { # Only if there are future time points
        idx_future_points <- (i + 1):N
        times_future <- time_line[idx_future_points]
        b_future_old_vals <- bnd_old[idx_future_points] # x2 for the kernel at each u_k
        
        # c1 and c2 for K2 depend on times_future (u_k)
        # Ensure vectorization of trend and slope or use sapply
        trend_vals_future <- trend(times_future)
        if(!is.vector(trend_vals_future) && length(times_future)>0) trend_vals_future <- sapply(times_future, trend)
        slope_vals_future <- slope(times_future)
        if(!is.vector(slope_vals_future) && length(times_future)>0) slope_vals_future <- sapply(times_future, slope)
        
        c1_K2_vec <- discount * strike + trend_vals_future
        c2_K2_vec <- discount - slope_vals_future
        
        K2_vec <- numeric(length(times_future))
        for (k_fut in 1:length(times_future)) {
          u_k <- times_future[k_fut]
          b_at_u_k <- b_future_old_vals[k_fut]
          c1_val_k <- c1_K2_vec[k_fut]
          c2_val_k <- c2_K2_vec[k_fut]
          
          K2_vec[k_fut] <- boundary_kernel_direct(
            c1 = c1_val_k, c2 = c2_val_k,
            t1 = t_current, x1 = b_current_old,
            t2 = u_k, x2 = b_at_u_k,
            slope_f = slope, trend_f = trend, vol_f = vol,
            discount_val = discount, subdivisions = trapz_subdivs,
            strike_price_for_log = strike
          )
        }
        # delta is of length N-1. delta[k] = time_line[k+1]-time_line[k]
        # The sum is K2(u_j) * (t_{j+1}-t_j) or similar.
        # The original delta: delta <- time_line[2:N] - time_line[1:(N-1)]
        # If we sum from i to N-1 for delta: delta[i:(N-1)]
        # This corresponds to the segments (t_i, t_{i+1}), ..., (t_{N-1}, t_N)
        # The length of K2_vec is length((i+1):N) = N-i.
        # The length of delta[i:(N-1)] is (N-1)-i+1 = N-i. They match.
        sum_K2_delta_val <- sum(K2_vec * delta[i:(N-1)])
      }
      
      bnd[i] <- strike - K1 - sum_K2_delta_val
      
      if(bnd[i] > strike && !isTRUE(all.equal(bnd[i], strike, tolerance = 1e-9))){
        # print(paste0("Clipping bnd[",i,"]=",bnd[i]," at iter ", j))
        bnd[i] <- strike 
      } else if (bnd[i] > strike) { 
        bnd[i] <- strike
      }
    } # End for loop i
    
    e <- sum((bnd - bnd_old)^2, na.rm=TRUE) # na.rm just in case
    # print(paste0("Iteration: ", j, ", error: ", sprintf("%.2e", e))) # More verbose
    if(is.na(e) || is.infinite(e)) {
      warning(paste("Error became NA/Inf in iteration",j,". Stopping."))
      # Recover the last valid boundary if possible
      if(exists("last_valid_bnd")) bnd <- last_valid_bnd else bnd <- bnd_old 
      # You might also want to save bnd_old as last_valid_bnd in each iteration
      break 
    }
    if (j>1) last_valid_bnd <- bnd_old # Save the last stable boundary
    
    if (errors) er <- c(er, e)
    
    if (j > 200) { # Lower iteration limit for this slow version
      warning(paste("Picard iteration limit (200) reached. Error:", e))
      break
    }
  } # End while e > tol
  
  print(paste0("Converged/Stopped after ", j, " iterations with error: ", sprintf("%.2e", e)))
  
  if (errors) return(list(boundary = bnd, errors = er))
  
  return(bnd)
}