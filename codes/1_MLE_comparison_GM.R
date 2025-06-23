library(nleqslv)
library(pracma)
library(stats4)
library(sde)

alpha_int_dt <- function(t0,t1,slope_fun,subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  y <- slope_fun(time_grid)
  int_values <- cumtrapz(time_grid, y)
  return(splinefun(time_grid, int_values)) # no exponential
}

trend_int_dt <- function(trend_fun, alpha_int_fun, t0, t1,subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  func2 <- function(t) {
    return(trend_fun(t) /exp(alpha_int_fun(t)))  
  }
  y <- func2(time_grid)
  return(trapz(time_grid,y))
}

vol_int_dt <- function(vol_fun, alpha_int_fun, t0, t1,subdiv){
  time_grid <- seq(t0, t1, l = subdiv)  # sub-partition
  func3 <- function(t) {
    return((vol_fun(t) /exp(alpha_int_fun(t)))^2)  
  }
  y = func3(time_grid)
  return(trapz(time_grid,y))
}

gm_mean_dt <- function(x0, t1, alpha_int_fun, trend_int){
  return(exp(alpha_int_fun(t1))*(x0 + trend_int))
}

gm_var_dt <- function(t1, alpha_int_fun, vol_int){
  return(exp(2*alpha_int_fun(t1))*vol_int)
}

log_likelihood4 <- function(alpha, beta, sigma, T, gm_paths, subdiv, slope_fun, trend_fun, vol_fun){
  
  # convert different many formats to matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Other formats
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # one vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Matrix and data.frames 
    gm_paths_matrix <- as.matrix(gm_paths)
  }
  # Oper file (optional)
  # write("x1\tt1\tx0\texp\ttrend\tvol\tmean\tvariance\tdnorm", file = "likelihood_mio4.txt")
  
  n <- nrow(gm_paths_matrix) 
  num_paths <- ncol(gm_paths_matrix)
  time_grid <- seq(0, T, length.out = n) 
  
  params <- list(alpha = alpha, beta = beta, sigma = sigma)
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t <- function(t) do.call(vol_fun, c(list(t = t), params))
  
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
      
      alpha_int_fun <- alpha_int_dt(t0,t1,slope_fun_t,subdiv) 
      exp_value <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1,subdiv)
      vol_value <- vol_int_dt(vol_fun_t, alpha_int_fun, t0, t1,subdiv)
      mean <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      variance <- gm_var_dt(t1, alpha_int_fun, vol_value)
      value <- dnorm(x1, mean, sqrt(variance), log = TRUE)
      path_logL <- path_logL + value
      
      # line for the .txt file
      # linea <- paste(x1, t1, x0,exp_value, trend_value, vol_value, mean, variance, value, sep = "\t")
     # write(linea, file = "likelihood_mio4.txt", append = TRUE, ncolumns = 1)
    }
    
    total_logL <- total_logL + path_logL
  }
  
  # Return negative sum (for minimization)
  return(-total_logL)
}

lk_estimates <-  function(subdivision, X, T, theta_test, slope_fun, trend_fun, vol_fun){
  opt_result <- mle(function(alpha,beta,sigma) log_likelihood4(alpha, beta, sigma, T, X, subdivision, slope_fun, trend_fun, vol_fun),
                    start=list(alpha=theta_test[1],beta= theta_test[2], sigma= theta_test[3]),
                    method="L-BFGS-B", lower=c(-Inf,-Inf,0.1))
  return(opt_result)
}

dmean_dtheta1 <- function(trend_fun, alpha_int_fun, exp_value, trend_value, dalpha_dtheta1, dtrend_dtheta1, t0, t1, x0, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  
  # integral function of dalpha_dtheta1
  
  y1 <- dalpha_dtheta1(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta1_int_fun <- splinefun(time_grid,int_values)
  
  # integrand of second term dmean_dtheta1
  
  func2 <- function(t){
    return(trend_fun(t)*dalpha_dtheta1_int_fun(t)/exp(alpha_int_fun(t)))
  }
  y2 <- func2(time_grid)
  
  # integrand of third term dmean_dtheta1
  func3 <- function(t){
    return(dtrend_dtheta1(t)/exp(alpha_int_fun(t)))
  }
  y3 <- func3(time_grid)
  
  # dmean_dtheta1
  return(exp_value*(dalpha_dtheta1_int_fun(t1)*(x0 + trend_value) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

dvar_dtheta1 <- function(exp_value, vol_value, vol_fun, alpha_int_fun, dalpha_dtheta1, dvol_dtheta1, t0, t1, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  
  # integral function of dalpha_dtheta1
  
  y1 <- dalpha_dtheta1(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta1_int_fun <- splinefun(time_grid,int_values) 
  
  # integrand second term of dvar_dtheta1
  
  func2 <- function(t){ 
    return((vol_fun(t)^2)*dalpha_dtheta1_int_fun(t)/(exp(alpha_int_fun(t)))^2)
  }
  y2 <- func2(time_grid)
  
  # integrand third term of dvar_dtheta1 
  func3 <- function(t){
    return(vol_fun(t)*dvol_dtheta1(t)/(exp(alpha_int_fun(t)))^2)
  }
  y3 <- func3(time_grid)
  
  return(2*exp_value^2*(vol_value*dalpha_dtheta1_int_fun(t1) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

dmean_dtheta2 <- function(trend_fun, alpha_int_fun, exp_value, trend_value, dalpha_dtheta2, dtrend_dtheta2, t0, t1, x0, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  
  # integral function of dalpha_dtheta2
  
  y1 <- dalpha_dtheta2(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta2_int_fun <- splinefun(time_grid,int_values)
  
  # integrand second term of dmean_dtheta2
  
  func2 <- function(t){
    return(trend_fun(t)*dalpha_dtheta2_int_fun(t)/exp(alpha_int_fun(t)))
  }
  y2 <- func2(time_grid)
  
  # integrand third termn of dmean_dtheta2
  func3 <- function(t){
    return(dtrend_dtheta2(t)/exp(alpha_int_fun(t)))
  }
  y3 <- func3(time_grid)
  
  return(exp_value*(dalpha_dtheta2_int_fun(t1)*(x0 + trend_value) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

dvar_dtheta2 <- function(exp_value, vol_value, vol_fun, alpha_int_fun, dalpha_dtheta2, dvol_dtheta2, t0, t1, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # sub-partition
  
  # integrand function of dalpha_dtheta2
  
  y1 <- dalpha_dtheta2(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta2_int_fun <- splinefun(time_grid,int_values) 
  
  # # integrand of second term dvar_dtheta2
  
  func2 <- function(t){ 
    return((vol_fun(t)^2)*dalpha_dtheta2_int_fun(t)/(exp(alpha_int_fun(t)))^2)
  }
  y2 <- func2(time_grid)
  
  # integrand of third term dvar_dtheta2
  
  func3 <- function(t){
    return(vol_fun(t)*dvol_dtheta2(t)/(exp(alpha_int_fun(t)))^2)
  }
  y3 <- func3(time_grid)
  
  return(2*exp_value^2*(vol_value*dalpha_dtheta2_int_fun(t1) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

dmean_dtheta3 <- function(trend_fun, alpha_int_fun, exp_value, trend_value, dalpha_dtheta3, dtrend_dtheta3, t0, t1, x0, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # la calculo con 10 iteraciones
  
  # # integral function of dalpha_dtheta3
  
  y1 <- dalpha_dtheta3(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta3_int_fun <- splinefun(time_grid,int_values)
  
  # integrand of second term dmean_dtheta3
  
  func2 <- function(t){
    return(trend_fun(t)*dalpha_dtheta3_int_fun(t)/exp(alpha_int_fun(t)))
  }
  y2 <- func2(time_grid)
  
  # integrand of third term dmean_dtheta3
  
  func3 <- function(t){
    return(dtrend_dtheta3(t)/exp(alpha_int_fun(t)))
  }
  y3 <- func3(time_grid)

  return(exp_value*(dalpha_dtheta3_int_fun(t1)*(x0 + trend_value) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

dvar_dtheta3 <- function(exp_value, vol_value, vol_fun, alpha_int_fun, dalpha_dtheta3, dvol_dtheta3, t0, t1, subdiv){
  time_grid <- seq(t0,t1,l=subdiv) # la calculo con 10 iteraciones
  
  # # integrand function of dalpha_dtheta3
  
  y1 <- dalpha_dtheta3(time_grid)
  int_values <- cumtrapz(time_grid, y1)
  dalpha_dtheta3_int_fun <- splinefun(time_grid,int_values) 
  
  # # integrand of second term dvar_dtheta3
  
  func2 <- function(t){ 
    return((vol_fun(t)^2)*dalpha_dtheta3_int_fun(t)/(exp(alpha_int_fun(t)))^2)
  }
  y2 <- func2(time_grid)
  
  # integrand of third term dvar_dtheta3
  
  func3 <- function(t){
    return(vol_fun(t)*dvol_dtheta3(t)/(exp(alpha_int_fun(t)))^2)
  }
  y3 <- func3(time_grid)

  return(2*exp_value^2*(vol_value*dalpha_dtheta3_int_fun(t1) - trapz(time_grid,y2) + trapz(time_grid,y3)))
}

Fisher_matrix_der_exact_complete <- function(T, gm_paths, subdiv, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                             dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3){
  # convert different many formats to matrix
  if (inherits(gm_paths, c("ts", "mts", "xts", "zoo"))) {
    gm_paths_matrix <- as.matrix(gm_paths)
  } else if (is.list(gm_paths)) {
    # Other formats
    gm_paths_matrix <- do.call(cbind, lapply(gm_paths, function(x) {
      if (inherits(x, c("ts", "mts", "xts", "zoo"))) {
        as.vector(as.matrix(x))
      } else {
        as.vector(x)
      }
    }))
  } else if (is.vector(gm_paths)) {
    # one vector
    gm_paths_matrix <- matrix(gm_paths, ncol = 1)
  } else {
    # Matrix and data.frames 
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
  
  # open file (optional) 
  #write("t0\tt1\tx0\tx1\tdlogL_dtheta1\tdlogL_dtheta2\tdlogL_dtheta3", file = "Fisher_matrix_der_exacta.txt")
  
  matrix_Fisher_der_exact <- matrix(0, nrow = 3, ncol = 3)
  
  for (path_num in 1:num_paths) {
    gm_path <- gm_paths_matrix[, path_num]
    dlogL_dtheta1 <- 0
    dlogL_dtheta2 <- 0
    dlogL_dtheta3 <- 0
    
    for (i in 2:n) {
      t0 <- time_grid[i - 1]
      t1 <- time_grid[i]
      x1 <- gm_path[i]
      x0 <- gm_path[i - 1]
      
      ## (log_likelihood4) 
      alpha_int_fun <- alpha_int_dt(t0,t1,slope_fun,subdiv) 
      exp_value <- exp(alpha_int_fun(t1)) # phi(t0,t1)
      trend_value <- trend_int_dt(trend_fun, alpha_int_fun, t0, t1,subdiv)
      vol_value <- vol_int_dt(vol_fun, alpha_int_fun, t0, t1,subdiv)
      mean <- gm_mean_dt(x0, t1, alpha_int_fun, trend_value)
      variance <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      ## New
      # deriv theta1
      dmean_dtheta1_value <- dmean_dtheta1(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta1, dtrend_dtheta1, t0, t1, x0, subdiv)
      dvar_dtheta1_value <- dvar_dtheta1(exp_value, vol_value,vol_fun, alpha_int_fun, dslope_dtheta1, dvol_dtheta1, t0, t1, subdiv)
      value_dl_dtheta1 <- 0.5*dvar_dtheta1_value*(-1/variance + (x1 - mean)^2/(variance^2)) + (x1 - mean)*dmean_dtheta1_value/variance
      dlogL_dtheta1 <- dlogL_dtheta1 + value_dl_dtheta1 
      
      # deriv theta2
      dmean_dtheta2_value <- dmean_dtheta2(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta2, dtrend_dtheta2, t0, t1, x0, subdiv)
      dvar_dtheta2_value <- dvar_dtheta2(exp_value, vol_value,vol_fun, alpha_int_fun, dslope_dtheta2, dvol_dtheta2, t0, t1, subdiv)
      value_dl_dtheta2 <- 0.5*dvar_dtheta2_value*(-1/variance + (x1 - mean)^2/(variance^2)) + (x1 - mean)*dmean_dtheta2_value/variance
      dlogL_dtheta2 <- dlogL_dtheta2 + value_dl_dtheta2 
      
      # deriv theta3 
      dmean_dtheta3_value <- dmean_dtheta3(trend_fun, alpha_int_fun, exp_value, trend_value, dslope_dtheta3, dtrend_dtheta3, t0, t1, x0, subdiv)
      dvar_dtheta3_value <- dvar_dtheta3(exp_value, vol_value, vol_fun, alpha_int_fun, dslope_dtheta3, dvol_dtheta3, t0, t1, subdiv)
      value_dl_dtheta3 <- 0.5*dvar_dtheta3_value*(-1/variance + (x1 - mean)^2/(variance^2)) + (x1 - mean)*dmean_dtheta3_value/variance
      dlogL_dtheta3 <- dlogL_dtheta3 + value_dl_dtheta3 
      
      # Copy line
      # linea <- paste(t0, t1, x0, x1, value_dl_dtheta1, value_dl_dtheta2, value_dl_dtheta3, sep = "\t")
      # write(linea, file = "Fisher_matrix_der_exacta.txt", append = TRUE, ncolumns = 1)
      
      grad <- c(value_dl_dtheta1, value_dl_dtheta2, value_dl_dtheta3)
      matrix_Fisher_der_exact <- matrix_Fisher_der_exact + tcrossprod(grad)
      
    }
    
    output$dl_dtheta1_total <- dlogL_dtheta1
    output$dl_dtheta2_total <- dlogL_dtheta2
    output$dl_dtheta3_total <- dlogL_dtheta3
    
    output$Fisher_matrix = matrix_Fisher_der_exact/n
  }
  
  return(output)
}
# wrapper function for nleqslv
gradient_wrapper <- function(theta, T, X, subdivision, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                             dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3) {
  
  params <- list(alpha = theta[1], beta = theta[2], sigma = theta[3])
  
  slope_fun_t <- function(t) do.call(slope_fun, c(list(t = t), params))
  trend_fun_t <- function(t) do.call(trend_fun, c(list(t = t), params))
  vol_fun_t <- function(t) do.call(vol_fun, c(list(t = t), params))
  
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
  
  
  result <- Fisher_matrix_der_exact_complete(T, X, subdivision, slope_fun = slope_fun_t,
                                             trend_fun = trend_fun_t, vol_fun = vol_fun_t, dslope_dtheta1 = dslope_dtheta1_t,
                                             dslope_dtheta2 = dslope_dtheta2_t, dslope_dtheta3 = dslope_dtheta3_t,
                                             dtrend_dtheta1 = dtrend_dtheta1_t, dtrend_dtheta2 = dtrend_dtheta2_t,
                                             dtrend_dtheta3 = dtrend_dtheta3_t, dvol_dtheta1 = dvol_dtheta1_t, 
                                             dvol_dtheta2 = dvol_dtheta2_t, dvol_dtheta3 = dvol_dtheta3_t)

  return(c(result$dl_dtheta1_total, result$dl_dtheta2_total, result$dl_dtheta3_total))
}

MLE_grad_loglikelihood <- function(initial_guess,T, X, subdivision, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3, 
                                   dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3, method = "Newton", global = "cline"){
  sol <- nleqslv(initial_guess, function(theta) gradient_wrapper(theta, T, X, subdivision, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                                                                 dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3), method = "Newton", global = "cline")
  # ## Methods y global strategies 
  # methods <- c("Newton", "Broyden")
  # globals <- c("cline", "qline", "gline", "pwldog", "dbldog", "hook", "none")
  # 
  # results <- list()
  # 
  # # loop
  # cat("Resultados de nleqslv por combinación método/global:\n")
  # for (m in methods) {
  #   for (g in globals) {
  #     res <- tryCatch({
  #       sol <- nleqslv(initial_guess, function(theta) gradient_wrapper(theta, T, X, subdivision, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
  #             dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3), method = "Newton", global = "cline")
  #       list(method = m, global = g, x = sol$x, conv = sol$termcd)
  #     }, error = function(e) {
  #       list(method = m, global = g, x = NA, conv = NA)
  #     })
  #     
  #     results[[paste(m, g, sep = "_")]] <- res
  #     
  #     cat(sprintf("Método: %-7s  Global: %-7s  Solución: (%.4f, %.4f, %.4f)  Convergencia: %d\n",
  #                 m, g,
  #                 ifelse(is.na(res$x[1]), NA, res$x[1]),
  #                 ifelse(is.na(res$x[2]), NA, res$x[2]),
  #                 ifelse(is.na(res$x[3]), NA, res$x[3]),
  #                 res$conv))
  #   }
  # }
  return(sol)
}

simulate_gm_process <- function(n_paths, n_steps, delta, x0, slope_fun_t, trend_fun_t, vol_fun_t, alpha_gm, beta_gm, sigma_gm, subdiv) {
  # n_steps, in the number of steps it takes, NOT including x0 
  T <- delta*n_steps
  time_grid <- seq(0, T, length.out = n_steps+1)
  paths <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths_exact <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths[1, ] <- x0
  paths_exact[1, ] <- x0
  
  # errors
  error_log <- list(
    negative_variance = 0,
    na_mean = 0,
    infinite_mean = 0,
    na_sd = 0,
    infinite_sd = 0,
    error_mean = 0, 
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
      
      # # Hull-White model periodic mean dificil
      # mean_exact <- (x_ant + alpha_gm*(exp(t1 - t0) - 1) + beta_gm*(exp(t1 - t0)*(sin(sigma_gm*t1) - sigma_gm*cos(sigma_gm*t1)) - (sin(sigma_gm*t0) - sigma_gm*cos(sigma_gm*t0)))/(1 + sigma_gm^2))/exp(t1 - t0)
      # 
      # var_exact <- 0.5*(1 - exp(-2*(t1 - t0)))
      
      # Hull-White model periodic mean version 1
      mean_exact <-exp(alpha_gm*(t1 - t0))*(x_ant_exact + beta_gm*(1 - exp(-alpha_gm*(t1 - t0)))/alpha_gm + ((alpha_gm*sin(sigma_gm*t0) + sigma_gm*cos(sigma_gm*t0)) - exp(alpha_gm*(t1 - t0))*(alpha_gm*sin(sigma_gm*t1) + sigma_gm*cos(sigma_gm*t1)))/(alpha_gm^2 + sigma_gm^2))

      var_exact <- (exp(2*alpha_gm*(t1 - t0)) - 1)/(2*alpha_gm)
      
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

      # OU time-dependent volatility dificil
      # mean_exact <- exp(-(t1 - t0))*(x_ant_exact + alpha_gm*(exp(t1 - t0) - 1))
      #
      # var_exact <- beta_gm^2*exp(-2*t1)*(exp(2*t1*(1 - sigma_gm)) - exp(2*t0*(1 - sigma_gm)))/(2*(1 - sigma_gm))

      # Hull White model logarithmic mean (very complex) and rational volatility (very complex)
      # mean_exact <- 0
      # var_exact <- 1
      
      alpha_int_fun <- alpha_int_dt(t0, t1, slope_fun_t,subdiv) 
      exp_value <- exp(alpha_int_fun(t1))
      trend_value <- trend_int_dt(trend_fun_t, alpha_int_fun, t0, t1,subdiv)
      vol_value <- vol_int_dt(vol_fun_t, alpha_int_fun, t0, t1, subdiv)
      mean_gm <- gm_mean_dt(x_ant, t1, alpha_int_fun, trend_value)
      var_gm <- gm_var_dt(t1, alpha_int_fun, vol_value)
      
      error_mean_path[i] <- mean_gm - mean_exact
      error_var_path[i] <- var_gm - var_exact
      
      # checks
      if (is.na(mean_gm)) {
        cat("NaN en mean, en la time step ",i, "\r")
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
        var_gm <- 0  # Force 0
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
# Plot GM 
plot_gm_paths <- function(simulated_paths, 
                          title = "GM Process Paths", 
                          plot_exact = FALSE,
                          ...) {
  n_paths <- ncol(simulated_paths$paths)
  
  show_exact_only <- !plot_exact || n_paths > 1
  show_both <- plot_exact && n_paths == 1
  
  # Calculate ranges
  y_range <- if(show_exact_only) {
    range(simulated_paths$paths_exact, na.rm = TRUE)
  } else if(show_both) {
    range(c(simulated_paths$paths, simulated_paths$paths_exact), na.rm = TRUE)
  } else {
    range(simulated_paths$paths, na.rm = TRUE)
  }
  
  # empty plot
  plot(simulated_paths$time, 
       if(show_exact_only) simulated_paths$paths_exact else simulated_paths$paths[,1], 
       type = "n", 
       ylim = y_range,
       xlab = "Time", ylab = "OU Bridge", main = title,
       ...)
  
  # Case 1: only paths_exact (default)
  if(show_exact_only && !is.null(simulated_paths$paths_exact)) {
    lines(simulated_paths$time, simulated_paths$paths_exact,
          col = "black", lwd = 3, lty = 1)
  }
  
  # Case 2: both(plot_exact=TRUE and n_paths=1)
  if(show_both) {
    # simulated trayectory
    lines(simulated_paths$time, simulated_paths$paths[,1], 
          col = "red", lwd = 2, lty = 1)
    
    # exact trayectory
    lines(simulated_paths$time, simulated_paths$paths_exact,
          col = "black", lwd = 3, lty = 2)
    
  }
  
  # Case 3: multiple simulated paths (n_paths > 1)
  if(n_paths > 1 && !show_exact_only) {
    for (j in 1:n_paths) {
      lines(simulated_paths$time, simulated_paths$paths[,j], 
            col = rainbow(n_paths)[j], lwd = 2)
    }
  }
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
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) # cuando tiene t esque solo depende de t y no del parametro
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
sigma_true <- 1
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

delta <- 1
n_steps <- 200 
T <- delta*n_steps
n_paths <- 1  
x0 <- 1
subdivision <- 100

# Test parameters
alpha_test <- -0.5 # mean reversion
beta_test <- 0.5
sigma_test <- 1
theta_test <- c(alpha_test, beta_test, sigma_test)

simulated_path <- simulate_gm_process(n_paths, n_steps, delta, x0, slope_fun_t, trend_fun_t, vol_fun_t, alpha_true, beta_true, sigma_true, subdivision)

# Show simulation errors
print(simulated_path$errors)

plot_gm_paths(simulated_path, plot_exact = TRUE)
# Uncomment for Hull White squared root, Hull White with logarithmic mean and rational volatility 
X <- simulated_path$paths
#X <- simulated_path$paths_exact
#X <- sde.sim(model = "OU", theta = c(beta_true, -alpha_true, sigma_true), N = n_steps, delta = delta, M = 1) # (beta, -alpha, sigma)
# First method
# opt_result <- tryCatch({
#   lk_estimates(subdivision, X, T, theta_test = theta_test,
#                slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
# }, error = function(e) {
#   message("Error in first method, iteration: ", e$message)
#   NULL
# })
# summary(opt_result@coef)
# 
# # estimated parameters
# alpha_est <- opt_result@coef[1]
# beta_est <- opt_result@coef[2]
# sigma_est <- opt_result@coef[3]
# 
# print("value theta1 optim")
# print(alpha_est)
# print("value theta2 optim")
# print(beta_est)
# print("Value sigma optim")
# print(sigma_est)

# Second method
result <- tryCatch({
  MLE_grad_loglikelihood(initial_guess = c(alpha_test, beta_test, sigma_test),T, X, subdivision, slope_fun = slope_fun,
                         trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2,
                         dslope_dtheta3 = dslope_dtheta3, dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
                         dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3, method = "Newton", global = "cline")
}, error = function(e) {
  message("Error in second method, iteration: ", e$message)
  NULL
})
print(result$x)

# Estimated parameters
alpha_est <- result$x[1]
beta_est <- result$x[2]
sigma_est <- result$x[3]

print("value theta1 grad")
print(alpha_est)
print("value theta2 grad")
print(beta_est)
print("Value sigma grad")
print(sigma_est)