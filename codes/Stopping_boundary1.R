####
## - The function "boundary" computes the optimal stopping boundary of the (exponentially) 
## discounted optimal stopping problem with finite horizon and gain function 
## G(x) = (strike - x)^+. The underlying process with parameters is a 
## time-dependent Ornstein-Uhlenbeck given by the SDE
##            dX_t = (slope(t)X_t + trend(t)) + vol(t)dW_t.
## - The boundary is computed by running a Picard iteration algorithm,
## that stops when the L2 distance between consecutive boundaries is less than
## tol, and solves the the free-boundary equation.
## - The boundary is computed at the time points provided in time_line
## - If errors == TRUE, a vector of the errors of the Picard scheme is provided
## alongside the boundary
####
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0, 
                      boundary_kernel, slope, trend, vol, errors = FALSE) {
  
  # Escribir encabezado en el archivo
  write("c1\tc2\ti1\tx1\ti2\tI3[i1]\tI3[i2]\tI1[i2]\tx2\tslope(i1)\ttrend(i1)\tvol(i1)\tmarginal_mean\tmarginal_var", file = "stopping_boundary_Abel1.txt")
  write("j\ti\tK1\tK2\ti2nbd[i]", file = "stopping_boundary_Abel2.txt")
  N <- length(time_line)       # Partition length
  expiration <- time_line[N]        # Expiration date
  delta <- time_line[2:N] - time_line[1:(N-1)]  # Length step vector
  
  # Creating errors vector if required
  if (errors) er <- c()
  
  # Pre-allocating boundary
  # bnd <- rep(min((slope(expiration) * pull(expiration) + discount * strike) /
  #                  (slope(expiration) + discount), strike), N)
  bnd <- rep(min((trend(expiration) + discount * strike)/(discount - slope(expiration)), strike), N)
  
  # Auxiliary pre-computations for marginal_mean and marginal_var
  f1 <- function(s) Vectorize(function(t) {
    exp(integrate(slope, lower = t, upper = s, 
                   subdivisions = 10000, rel.tol = 1e-10)$value)
  })
  I1 <- f1(expiration)(time_line)
  f2 <- function(t) trend(t) * f1(expiration)(t)
  I2 <- sapply(time_line, function(s) {
    integrate(f2, lower = s, upper = expiration, 
              subdivisions = 10000, rel.tol = 1e-10)$value
  })
  f3 <- function(t) (f1(expiration)(t) * vol(t))^2
  I3 <- sapply(time_line, function(s) {
    integrate(f3, lower = s, upper = expiration, 
              subdivisions = 10000, rel.tol = 1e-10)$value
  })
  
  # Kernel definition
  boundary_kernel <- function(c1, c2, i1, x1, i2, x2) {
    
    # Compute the marginal mean
    marginal_mean <- (x1 * I1[i1] + (I2[i1] -  I2[i2])) / I1[i2]
    # Compute the marginal standard deviation
    #cat("(I3[i1] - I3[i2]) / I1[i2]^2 =",(I3[i1] - I3[i2]) / I1[i2]^2, "\n")
    marginal_var <- (I3[i1] - I3[i2]) / I1[i2]^2
    marginal_sd <- sqrt((I3[i1] - I3[i2]) / I1[i2]^2)
    
    # Crear una línea con los valores separados por tabulación
    linea1 <- paste(c1, c2, i1, x1, i2,I3[i1],I3[i2],I1[i2], x2,slope(i1), trend(i1), vol(i1), marginal_mean, marginal_var, sep = "\t")
    # Escribir en el archivo, agregando nuevas líneas en cada llamada
    write(linea1, file = "stopping_boundary_Abel1.txt", append = TRUE, ncolumns = 1)
    
    # Compute standardized values
    x2 <- (x2 - marginal_mean) / marginal_sd
    # Compute normal distribution and density
    normal_dist <- pnorm(x2, mean = 0, sd = 1, lower.tail = T)
    normal_dens <- dnorm(x2, mean = 0, sd = 1)
    # Evaluate Kernel
    K <- exp(-discount * (time_line[i2] - time_line[i1])) * 
      ((c1 - c2 * marginal_mean) * normal_dist + c2 * marginal_sd * normal_dens) 
    
    return(K)
    
  }
  
  # Boundary computation
  e <- 1  # error in the while loop
  # Fixed point algorithm
  
  j <- 0
  while (e > tol) {
    j <- j + 1
    
    bnd_old <- bnd
    
    # 
    for (i in (N - 1):1) {  
      
      # print(paste0("updating boundary at t_", i, " = ", time_line[i]))
      # Evaluate the kernel
      #cat("Ahora K1", "\n")
      K1 <- boundary_kernel(c1 = strike, c2 = 1, 
                            i1 = i, x1 = bnd_old[i], 
                            i2 = N, x2 = strike)
      #cat("Ahora K2", "\n")
      K2 <- boundary_kernel(c1 = discount * strike + trend(time_line[(i+1):N]),
                            c2 = discount - slope(time_line[(i + 1):N]),
                            i1 = i, x1 = bnd_old[i], 
                            i2 = (i + 1):N, x2 = bnd_old[(i + 1):N])
      
      # Update the boundary at t_present
      bnd[i] <- strike - K1 - sum(K2 * delta[i:(N - 1)])
      # # Crear una línea con los valores separados por tabulación
      linea2 <- paste(j, i, K1, K2, bnd[i], sep = "\t")
      # # Escribir en el archivo, agregando nuevas líneas en cada llamada
      write(linea2, file = "stopping_boundary_Abel2.txt", append = TRUE, ncolumns = 1)
      
      if(any(bnd[i] > strike)){
        print("Imposible: Boundary above the strike price")
      }
      
    }
    
    # absolute L2 error
    e <- sum((bnd - bnd_old)^2)
    print(paste0("error: ", e))
    if (errors) er <- c(er, e)
    
  }
  
  if (errors) return(list(boundary = bnd, errors = er))
  
  return(bnd)
  
}


expiration <- 1
partition_length <- 100
delta <- 3
expiration <- partition_length*delta
time_line <- seq(0, expiration, l = partition_length+1)
strike <- 0
discount <- 0 

alpha_est <- -1
beta_est <- 3 
sigma_est <- 2

slope_fun_t <- function(t) rep(alpha_est, length(t))
trend_fun_t <- function(t) rep(beta_est, length(t))
vol_fun_t <- function(t) rep(sigma_est, length(t))

bnd <- boundary(tol = 1e-3 ,strike = strike, time_line = time_line,
                     discount = discount, boundary_kernel = boundary_kernel,
                     slope = slope_fun_t, trend = trend_fun_t, vol = vol_fun_t, errors = TRUE)
err <- bnd$errors
cat("error =", err, "\n")
bnd <- bnd$boundary
cat("boundary =", bnd, "\n")
