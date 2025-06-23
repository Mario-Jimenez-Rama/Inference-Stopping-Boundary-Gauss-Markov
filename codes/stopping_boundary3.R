library(pracma)

# --- Función auxiliar para integración directa con trapz (sin cambios) ---
integrate_direct_trapz <- function(func_to_integrate, lower, upper, subdivisions, ...) {
  if (abs(lower - upper) < .Machine$double.eps^0.75) {
    return(0)
  }
  if (subdivisions < 2) {
    stop("integrate_direct_trapz requiere al menos 2 subdivisiones (puntos).")
  }
  grid <- seq(lower, upper, length.out = subdivisions)
  
  func_name_str <- deparse(substitute(func_to_integrate))
  # Si func_to_integrate es una de las que esperan un escalar como primer argumento de la grid.
  if (grepl("integrand_for_trend_term_optimized", func_name_str) || 
      grepl("integrand_for_vol_term_optimized", func_name_str) ) { 
    values <- sapply(grid, func_to_integrate, ...)
  } else { # Para funciones como slope, trend, vol que se asume son vectorizadas
    values <- func_to_integrate(grid, ...)
  }
  
  return(pracma::trapz(grid, values))
}

# --- Funciones para los integrandos OPTIMIZADAS (usan spline_F_slope) ---
integrand_for_trend_term_optimized <- function(u_mesh_val, current_t2, 
                                               spline_F_slope_func, # Cambiado: función spline
                                               trend_func) {
  # integral_slope_u_t2 se obtiene de la spline precalculada
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  return(trend_func(u_mesh_val) * exp(integral_slope_u_t2))
}

integrand_for_vol_term_optimized <- function(u_mesh_val, current_t2, 
                                             spline_F_slope_func, # Cambiado: función spline
                                             vol_func) {
  integral_slope_u_t2 <- spline_F_slope_func(current_t2) - spline_F_slope_func(u_mesh_val)
  term_in_paren <- vol_func(u_mesh_val) * exp(integral_slope_u_t2)
  return(term_in_paren^2)
}

# --- Kernel con cálculo OPTIMIZADO de integrales ---
boundary_kernel_optimized <- function(c1, c2, t1, x1, t2, x2, 
                                      spline_F_slope_func, # Cambiado: función spline
                                      trend_f, vol_f, 
                                      discount_val, subdivisions,
                                      # Ya no necesita slope_f directamente aquí para integrales
                                      # Pero sí para el logging si se quiere slope_f(t1)
                                      slope_f_for_log_only # Solo para logging
) {
  
  # 1. exp_int_slope_t1_t2 = exp(integral_{t1}^{t2} slope(u)du)
  # Se obtiene de la spline precalculada
  integral_slope_t1_t2 <- spline_F_slope_func(t2) - spline_F_slope_func(t1)
  exp_int_slope_t1_t2 <- exp(integral_slope_t1_t2)
  
  # 2. integral_trend_term 
  integral_trend_term <- integrate_direct_trapz(
    integrand_for_trend_term_optimized, t1, t2, subdivisions,
    current_t2 = t2, spline_F_slope_func = spline_F_slope_func, trend_func = trend_f
  )
  
  # 3. integral_vol_term (varianza directa)
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
  # linea <- paste(c1, c2, t1, x1, t2, x2,
  #                slope_f_for_log_only(t1), trend_f(t1), vol_f(t1),
  #                marginal_mean_val, current_marginal_var, sep = "\t")
  # try(write(linea, file = "stopping_boundary_Abel.txt", append = TRUE, ncolumns = 11), silent=TRUE)
  
  # --- Cálculo del valor del Kernel ---
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


# --- Función principal boundary (OPTIMIZADA) ---
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0,
                      slope, trend, vol, errors = FALSE,
                      trapz_subdivs = 100, 
                      cumtrapz_subdivs_slope = 1000 # Para la spline de slope
) {
  
  # --- Pre-cálculo de la integral de slope usando cumtrapz y spline ---
  # Rejilla fina desde el primer punto de time_line hasta expiration
  # Asegurarse de que min(time_line) es time_line[1] si está ordenado
  min_time <- min(time_line) # Usar min/max para robustez si time_line no está ordenado
  max_time <- max(time_line)
  
  # Asegurar que cumtrapz_subdivs_slope sea al menos 2
  safe_cumtrapz_subdivs_slope <- max(2, cumtrapz_subdivs_slope)
  fine_grid_for_slope_spline <- seq(min_time, max_time, length.out = safe_cumtrapz_subdivs_slope)
  
  slope_vals_on_fine_grid <- slope(fine_grid_for_slope_spline)
  # Asegurar que slope_vals_on_fine_grid es un vector del tamaño correcto
  if(!is.vector(slope_vals_on_fine_grid) || length(slope_vals_on_fine_grid) != length(fine_grid_for_slope_spline)) {
    # Si slope no es vectorizada, aplicar sapply.
    # Esto es importante si slope no es de la forma rep(const, length(t))
    slope_vals_on_fine_grid <- sapply(fine_grid_for_slope_spline, slope)
  }
  
  # F(x) = integral_{fine_grid_for_slope_spline[1]}^{x} slope(u) du
  cumulative_integral_slope_values <- pracma::cumtrapz(fine_grid_for_slope_spline, slope_vals_on_fine_grid)
  
  # Crear función spline para F(x)
  spline_F_slope <- splinefun(fine_grid_for_slope_spline, cumulative_integral_slope_values, method = "natural")
  # --- Fin pre-cálculo ---
  
  # Escribir encabezado en el archivo (si el logging está activo en el kernel)
  # write("c1\tc2\tt1\tx1\tt2\tx2\tslope(t1)\ttrend(t1)\tvol(t1)\tmarginal_mean\tmarginal_var", file = "stopping_boundary_Abel.txt")
  
  N <- length(time_line)
  expiration <- time_line[N] # Asume que el último punto es la expiración
  delta_t_steps <- time_line[2:N] - time_line[1:(N-1)]
  
  if (errors) er <- c()
  
  initial_bnd_val <- strike
  # Para slope(expiration), es mejor evaluar la función original, no la spline.
  slope_at_expiration <- slope(expiration)[1] # Tomar el primer elemento si devuelve vector
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
    
    print(paste("Iteración de Picard:", j)) # Para ver el progreso
    print(paste("error:", e)) # Para ver el progreso
    
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
        slope_vals_future <- slope(times_future) # No la spline, la función original
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
      warning(paste("Error se volvió NA/NaN/Inf en iteración",j,". Deteniendo. Usando última frontera válida."))
      bnd <- last_valid_bnd 
      if (errors) er <- c(er, e) 
      break 
    }
    if(j > 1 || (j==1 && !(is.na(e) || is.infinite(e) || is.nan(e)))) { # Asegurar que bnd_old es 'bueno'
      last_valid_bnd <- bnd_old 
    }
    
    
    if (errors) er <- c(er, e)
    
    # Reducir el límite de iteraciones si se sospecha que es muy lento inicialmente
    iter_limit <- 2000 
    if (j > iter_limit) { 
      warning(paste("Picard iteration limit (", iter_limit, ") reached. Error:", e))
      break
    }
  }
  
  print(paste0("Convergido/Detenido tras ", j, " iteraciones con error: ", sprintf("%.2e", e)))
  
  if (errors) return(list(boundary = bnd, errors = er))
  
  return(bnd)
}

simulate_gm_process <- function(n_paths, n_steps, delta, x0, slope_fun_t, trend_fun_t, vol_fun_t, alpha_gm, beta_gm, sigma_gm, subdiv) {
  # n_steps, son los pasos que da 
  T <- delta*n_steps
  time_grid <- seq(0, T, length.out = n_steps+1)
  paths <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths_exact <- matrix(NA, nrow = n_steps+1, ncol = n_paths)
  paths[1, ] <- x0
  paths_exact[1, ] <- x0
  
  # Contadores de errores
  error_log <- list(
    negative_variance = 0,
    na_mean = 0,
    infinite_mean = 0,
    na_sd = 0,
    infinite_sd = 0,
    error_mean = 0, # media de todos los errores medios en cada path
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
      
      # Validación de parámetros para rnorm()
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
        var_gm <- 0  # Forzar a 0 para evitar NaN
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
# Función mejorada para visualización
plot_gm_paths <- function(simulated_paths, 
                          title = "GM Process Paths", 
                          plot_exact = FALSE,
                          ...) {
  n_paths <- ncol(simulated_paths$paths)
  
  # Determinar qué curvas mostrar y calcular el rango Y adecuado
  show_exact_only <- !plot_exact || n_paths > 1
  show_both <- plot_exact && n_paths == 1
  
  # Calcular rango Y considerando las curvas relevantes
  y_range <- if(show_exact_only) {
    range(simulated_paths$paths_exact, na.rm = TRUE)
  } else if(show_both) {
    range(c(simulated_paths$paths, simulated_paths$paths_exact), na.rm = TRUE)
  } else {
    range(simulated_paths$paths, na.rm = TRUE)
  }
  
  # Configurar el plot vacío
  plot(simulated_paths$time, 
       if(show_exact_only) simulated_paths$paths_exact else simulated_paths$paths[,1], 
       type = "n", 
       ylim = y_range,
       xlab = "Time", ylab = "OU Bridge", main = title,
       ...)
  
  # Caso 1: Mostrar solo paths_exact (por defecto)
  if(show_exact_only && !is.null(simulated_paths$paths_exact)) {
    lines(simulated_paths$time, simulated_paths$paths_exact,
          col = "black", lwd = 3, lty = 1)
  }
  
  # Caso 2: Mostrar ambas (solo cuando plot_exact=TRUE y n_paths=1)
  if(show_both) {
    # Trayectoria simulada
    lines(simulated_paths$time, simulated_paths$paths[,1], 
          col = "red", lwd = 2, lty = 1)
    
    # Trayectoria exacta
    lines(simulated_paths$time, simulated_paths$paths_exact,
          col = "black", lwd = 3, lty = 2)
    
  }
  
  # Caso 3: Mostrar múltiples paths simulados (n_paths > 1)
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
alpha_true <- -1 #mean reversion
beta_true <- 1
sigma_true <- 0.5
slope_fun <- function(t, alpha,...) rep(alpha, length(t))
trend_fun <- function(t,beta, sigma,...) beta + sigma*sqrt(t)
vol_fun <- function(t,...) rep(1, length(t))

# functions with true parameters
slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
vol_fun_t <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)

## derivadas ,cuando tiene t esque solo depende de t y no del parametro
# derivadas de slope:
dslope_dtheta1 <- function(t,...) rep(1, length(t))
dslope_dtheta2 <- function(t,...) rep(0, length(t))
dslope_dtheta3 <- function(t,...) rep(0, length(t))

# derivadas de trend
dtrend_dtheta1 <- function(t,...) rep(0, length(t))
dtrend_dtheta2 <- function(t,...) rep(1, length(t))
dtrend_dtheta3 <- function(t,...) sqrt(t)

# derivadas de vol
dvol_dtheta1 <- function(t,...) rep(0, length(t))
dvol_dtheta2 <- function(t,...) rep(0, length(t))
dvol_dtheta3 <- function(t,...) rep(0, length(t))

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

# params <- list(alpha = alpha_est, beta = beta_est, sigma = sigma_est)

partition_length <- 25
delta <- 1
expiration <- delta*partition_length
time_line <- seq(0, expiration, l = partition_length)
strike <- 9
discount <- 0 

bnd <- boundary(tol = 1e-3 ,strike = strike, time_line = time_line,
                discount = discount, slope = slope_fun_t, trend = trend_fun_t,
                vol = vol_fun_t, errors = TRUE, trapz_subdivs = 100, cumtrapz_subdivs_slope = 1000)

err <- bnd$errors
cat("error =", err, "\n")
bnd <- bnd$boundary
cat("boundary =", bnd, "\n")

# Save as a RDS file
folder_name <- "rds_files"
boundary_file <- file.path(folder_name, paste0("boundary_irrational.rds"))
saveRDS(bnd, boundary_file)

