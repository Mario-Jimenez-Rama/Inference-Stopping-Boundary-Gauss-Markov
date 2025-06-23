# It takes too long 
#source("C:/Users/Mario Jimenez/Desktop/Máster Universitario en Matemática Aplicada y Computacional, UC3M/TFM/Codigos/likelihoodentero4.R")
library(pracma)
library(sde)
library(ggplot2)
library(gridExtra)
library(GGally)

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

pairwise_plot <- function(theta_true, N, delta, subdivision, theta_test, slope_fun, trend_fun, vol_fun, dslope_dtheta1, dslope_dtheta2, dslope_dtheta3,
                          dtrend_dtheta1, dtrend_dtheta2, dtrend_dtheta3, dvol_dtheta1, dvol_dtheta2, dvol_dtheta3){
  
  T <- N*delta
  data_100 <- list()

  #for (i in 1:100){
  for (i in 1:3){
    X <- sde.sim(model = "OU", theta = c(theta_true[2],-theta_true[1],theta_true[3]), N = N, delta = delta) # (beta, -alpha, sigma)
    
    # likelihood4
    result <- tryCatch({
      lk_estimates(subdivision, X, T, theta_test, slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
    }, error = function(e) {
      message("Error in estimation (iteration ", i, ") in data_100: ", e$message)
      NULL
    })
    
    params <- list(alpha = result@coef[1], beta = result@coef[2], sigma = result@coef[3])
    
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

    fisher_info_matrix <- tryCatch({
      Fisher_matrix_der_exact_complete(T, X, subdivision, slope_fun = slope_fun_t,
                                       trend_fun = trend_fun_t, vol_fun = vol_fun_t, dslope_dtheta1 = dslope_dtheta1_t,
                                       dslope_dtheta2 = dslope_dtheta2_t, dslope_dtheta3 = dslope_dtheta3_t,
                                       dtrend_dtheta1 = dtrend_dtheta1_t, dtrend_dtheta2 = dtrend_dtheta2_t, dtrend_dtheta3 = dtrend_dtheta3_t,
                                       dvol_dtheta1 = dvol_dtheta1_t, dvol_dtheta2 = dvol_dtheta2_t, dvol_dtheta3 = dvol_dtheta3_t)$Fisher_matrix
    }, error = function(e) {
      message("Error in Fisher matrix (iteration ", i, ") in data_100: ", e$message)
      NULL
    })
    
    if (is.null(fisher_info_matrix)) next  
    
    R <- tryCatch({
      chol(fisher_info_matrix)
    }, error = function(e) {
      message("Error in Cholesky descomposition (iteration ", i, ") in data_100: ", e$message)
      NULL
    })
    
    if (is.null(R)) next  
    
    point <- sqrt(N) * (R %*% (c(result@coef[1], result@coef[2], result@coef[3]) - theta_true))
    data_100[[length(data_100) + 1]] <- setNames(as.numeric(point), c("alpha", "beta", "sigma"))
    
  }
  
  # Generate data for sample_size = 1000
  data_1000 <- list()
  #for (i in 1:1000){
  for (i in 1:4){
    X <- sde.sim(model = "OU", theta = c(theta_true[2],-theta_true[1],theta_true[3]), N = N, delta = delta) # (beta, -alpha, sigma)
    
    # likelihood4
    result <- tryCatch({
      lk_estimates(subdivision, X, T, theta_test, slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
    }, error = function(e) {
      message("Error in estimation (iteration ", i, ") in data_1000: ", e$message)
      NULL
    })
    
    params <- list(alpha = result@coef[1], beta = result@coef[2], sigma = result@coef[3])
    
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
    
    fisher_info_matrix <- tryCatch({
      Fisher_matrix_der_exact_complete(T, X, subdivision, slope_fun = slope_fun_t,
                                       trend_fun = trend_fun_t, vol_fun = vol_fun_t, dslope_dtheta1 = dslope_dtheta1_t,
                                       dslope_dtheta2 = dslope_dtheta2_t, dslope_dtheta3 = dslope_dtheta3_t,
                                       dtrend_dtheta1 = dtrend_dtheta1_t, dtrend_dtheta2 = dtrend_dtheta2_t, dtrend_dtheta3 = dtrend_dtheta3_t,
                                       dvol_dtheta1 = dvol_dtheta1_t, dvol_dtheta2 = dvol_dtheta2_t, dvol_dtheta3 = dvol_dtheta3_t)$Fisher_matrix
    }, error = function(e) {
      message("Error Fisher matrix (iteration ", i, ") in data_1000: ", e$message)
      NULL
    })
    
    if (is.null(fisher_info_matrix)) next
    
    R <- tryCatch({
      chol(fisher_info_matrix)
    }, error = function(e) {
      message("Error in Cholesky descomposition (iteration ", i, ")in data_1000: ", e$message)
      NULL
    })
    
    if (is.null(R)) next  
    
    point <- sqrt(N) * (R %*% (c(result@coef[1], result@coef[2], result@coef[3]) - theta_true))
    data_1000[[length(data_1000) + 1]] <- setNames(as.numeric(point), c("alpha", "beta", "sigma"))
    
  }
  # Change to dataframe as list are more efficient when we dont know the size of the columns
  data_100 <- do.call(rbind, data_100) |> as.data.frame()
  data_1000 <- do.call(rbind, data_1000) |> as.data.frame()
  
  # Combine data
  data_100$sample_size <- "100"
  data_1000$sample_size <- "1000"
  combined_data <- rbind(data_100, data_1000)
  
  # Compute correlations
  cor_100 <- cor(data_100[,1:3])
  cor_1000 <- cor(data_1000[,1:3])
  
  cat("Correlations for sample_size = 100:\n")
  print(cor_100)
  cat("\nCorrelations for sample_size = 1000:\n")
  print(cor_1000)
  
  # Plot functions
  density_with_normal <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(aes(color = sample_size), alpha = 0.7) + 
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) +
      theme_minimal() +
      xlim(c(-6, 6))
  }
  
  scatterplots <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(aes(color = sample_size),size = 0.8, alpha = 0.7) +
      theme_minimal() +
      xlim(c(-6, 6)) +
      ylim(c(-6, 6))
  }
  
  upper_corr <- function(data, mapping, ...) {
    x_var <- as.character(mapping$x[2])
    y_var <- as.character(mapping$y[2])
    
    corr_100 <- round(cor_100[x_var, y_var], 4)
    corr_1000 <- round(cor_1000[x_var, y_var], 4)
    
    ggplot() + 
      annotate("text", x = 0.5, y = c(0.55, 0.5, 0.49, 0.45), 
               label = c("", paste(corr_100), paste(corr_1000), ""), 
               color = c("white", "red", "turquoise", "white"), size = 5, hjust = 0.5) +
      theme_minimal()
  }
  
  ggpairs(combined_data, columns = 1:3, 
          upper = list(continuous = wrap(upper_corr)),
          lower = list(continuous = scatterplots),
          diag = list(continuous = density_with_normal)
  )
}

N = 100 # number of time steps
M = 100 # number of paths
delta = 1 #time step
subdivision <- 70 

alpha_test <- -0.5
beta_test <- 1
sigma_test <- 1

alpha_true <- -1
beta_true <- 3
sigma_true <- 2

# OU stacionary
slope_fun <- function(t, alpha, ...) rep(alpha, length(t))
trend_fun <- function(t, beta,...) rep(beta, length(t))
vol_fun <- function(t, sigma,...) rep(sigma, length(t))

## derivadas
# derivadas de slope:
dslope_dtheta1 <- function(t,...) rep(1, length(t)) 
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


plot_combined <- pairwise_plot(c(alpha_true, beta_true, sigma_true), N, M, delta, subdivision, c(alpha_test, beta_test, sigma_test),
                               slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1,
                               dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3, dtrend_dtheta1 = dtrend_dtheta1,
                               dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3, dvol_dtheta1 = dvol_dtheta1,
                               dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3)
                               
print(plot_combined)