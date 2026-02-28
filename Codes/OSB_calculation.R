
M  = 2 # number of paths and their boundary inference
for(i in 2:M){
  N = 100 # number of time steps
  delta = 1 # time step
  X <- sde.sim(model = "OU", theta = c(3,1,2), N = N, delta = delta) # (beta, -alpha, sigma)
  subdivision <- 70 # depends on T and n; for OU it does not matter since the functions are constants
  der_tol <- 1e-4
  
  alpha_test <- -0.5
  beta_test <- 1
  sigma_test <- 1
  
  alpha <- -1
  beta <- 3
  sigma <- 2
  
  # OU stationary
  slope_fun <- function(t, alpha, ...) rep(alpha, length(t))
  trend_fun <- function(t, beta,...) rep(beta, length(t))
  vol_fun <- function(t, sigma,...) rep(sigma, length(t))
  
  ## derivatives
  # slope derivatives:
  dslope_dtheta1 <- function(t,...) rep(1, length(t)) 
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
  
  inference <- tryCatch({
    infer_boundary(X = X, delta = delta, z_alpha = 0.1, partition_length = N, strike = 0, discount = 0,
                   theta_test = c(alpha_test, beta_test, sigma_test), slope_fun = slope_fun,
                   trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1,
                   dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3, 
                   dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, 
                   dtrend_dtheta3 = dtrend_dtheta3, dvol_dtheta1 = dvol_dtheta1,
                   dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3,
                   der_tol = der_tol, subdivision = subdivision)
  }, error = function(e){
    message("Error in boundary inference (iteration ", i, "): ", e$message)
  })
  
  if (is.null(inference)) next  # If Cholesky failed, move to the next iteration
  
  # save the specified file index
  file_index <- i  # Replace with the desired file index
  folder_name <- "rds_files"
  
  # Define file paths
  x_path_file <- file.path(folder_name, paste0("path_OU_", file_index, ".rds"))
  inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_OU_", file_index, ".rds"))
  
  # save the data to the .rds files
  saveRDS(X, x_path_file)
  saveRDS(inference, inferred_boundary_file)