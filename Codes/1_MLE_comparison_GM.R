# OUB

# alpha_true <- 1
# beta_true <- 2
# sigma_true <- 2
# slope_fun <- function(t, alpha,...)  -alpha/tanh(alpha*(T-t))
# trend_fun <- function(t, alpha, beta, ...) xf*alpha/sinh(alpha*(T-t)) - beta*tanh(alpha*(T-t)/2)
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true)
# vol_fun_true   <- function(t) vol_fun(t, sigma = sigma_true)
#
# ## derivatives; when it has t it means it only depends on t and not on the parameter
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

# Hull-White model with periodic mean (difficult)
# alpha_true <- 3
# beta_true <- 1
# sigma_true <- 2
# slope_fun <- function(t,...) rep(-1, length(t)) # mean reverting
# trend_fun <- function(t, alpha, beta, sigma) alpha + beta*sin(sigma*t)
# vol_fun <- function(t,...) rep(1, length(t)) # unit volatility
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
#
# ## derivatives; when it has t it means it only depends on t and not on the parameter
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

# OU with time-dependent volatility (complex)
# alpha_true <- 1
# beta_true <- 0.5
# sigma_true <- 0.25
# slope_fun <- function(t,...) rep(-1, length(t)) # constant, mean reversion
# # trend_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t, alpha,...) rep(alpha, length(t))
# vol_fun <- function(t, beta, sigma,...) beta*exp(-sigma*t)
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
#
# # derivatives; when it has t it means it only depends on t and not on the parameter
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

# Hull-White model with square-root drift term
# alpha_true <- 1 # mean reversion
# beta_true <- 1
# sigma_true <- 0.5
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta, sigma,...) beta + sigma*sqrt(t)
# vol_fun <- function(t,...) rep(1, length(t))
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
#
# ## derivatives; when it has t it means it only depends on t and not on the parameter
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
# alpha_true <- -1 # mean reversion
# beta_true <- 3
# sigma_true <- 2
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t, beta,...) rep(beta, length(t))
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
#
# ## derivatives
# # slope derivatives:
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) # when it has t it means it only depends on t and not on the parameter
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
# dvol_dtheta3 <- function(t,...) rep(1, length(t))

# Hull-White periodic mean version 1
true_alpha <- -1 # mean reversion
true_beta  <- 3
true_sigma <- 0.5
slope_fun <- function(t, alpha,...) rep(alpha, length(t))
trend_fun <- function(t, beta, sigma,...) beta + sin(sigma*t)
vol_fun <- function(t,...) rep(1, length(t))

# functions with true parameters
slope_fun_true <- function(t) slope_fun(t, alpha = true_alpha, beta = true_beta, sigma = true_sigma)
trend_fun_true <- function(t) trend_fun(t, alpha = true_alpha, beta = true_beta, sigma = true_sigma)
vol_fun_true   <- function(t) vol_fun(t, alpha = true_alpha, beta = true_beta, sigma = true_sigma)

## derivatives
# slope derivatives:
dslope_dtheta1 <- function(t,...) rep(1, length(t))
dslope_dtheta2 <- function(t,...) rep(0, length(t))
dslope_dtheta3 <- function(t,...) rep(0, length(t))

# trend derivatives
dtrend_dtheta1 <- function(t,...) rep(0, length(t))
dtrend_dtheta2 <- function(t,...) rep(1, length(t))
dtrend_dtheta3 <- function(t, sigma,...) t*cos(sigma*t)

# vol derivatives
dvol_dtheta1 <- function(t,...) rep(0, length(t))
dvol_dtheta2 <- function(t,...) rep(0, length(t))
dvol_dtheta3 <- function(t,...) rep(0, length(t))

# Hull-White periodic mean version 2
# alpha_true <- 3
# beta_true <- 3
# sigma_true <- 2
# slope_fun <- function(t, alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta, sigma,...) beta + sigma*sin(t)
# vol_fun <- function(t,...) rep(1, length(t))
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
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

# OU with time-dependent decaying volatility
# alpha_true <- -0.5
# beta_true <- 3
# sigma_true <- 0.5
# slope_fun <- function(t,alpha,...) rep(alpha, length(t)) # mean reversion
# trend_fun <- function(t,beta,...) rep(beta, length(t))
# # vol_fun <- function(t,sigma,...) 5*exp(-sigma*t)*sin(sigma*t) # 1st
# vol_fun <- function(t,sigma,...) 5*exp(-sigma*t) # 2nd
# # vol_fun <- function(t, beta, sigma,...) beta + 5*exp(-sigma*t)
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
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
# # dvol_dtheta3 <- function(t, sigma,...) -t*exp(-sigma*t)*sin(sigma*t) + t*exp(-sigma*t)*cos(sigma*t) # 1st
# dvol_dtheta3 <- function(t, sigma,...) -5*t*exp(-sigma*t) # 2nd

# Hull-White logarithmic mean
# alpha_true <- -1
# beta_true <- 2
# sigma_true <- 2
# # slope_fun <- function(t,...) rep(-1, length(t)) # mean reversion
# slope_fun <- function(t, alpha,...) rep(alpha, length(t)) # mean reversion # 3rd
# # trend_fun <- function(t,alpha, beta,...) alpha + log(1 + beta*t) # 1st
# # trend_fun <- function(t,alpha, beta,...) alpha + beta*log(1 + t) # 2nd
# trend_fun <- function(t,beta,...) beta*log(1 + t) # 3rd
# vol_fun <- function(t, sigma,...) rep(sigma, length(t))
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
#
# ## derivatives
# # slope derivatives:
# # dslope_dtheta1 <- function(t,...) rep(0, length(t)) # 1st, 2nd
# dslope_dtheta1 <- function(t,...) rep(1, length(t)) # 3rd
# dslope_dtheta2 <- function(t,...) rep(0, length(t))
# dslope_dtheta3 <- function(t,...) rep(0, length(t))
#
# # trend derivatives
# # dtrend_dtheta1 <- function(t,...) rep(1, length(t)) # 2nd
# dtrend_dtheta1 <- function(t,...) rep(0, length(t)) # 3rd
# # dtrend_dtheta2 <- function(t, beta,...) t/(1 + beta*t) # 1st
# dtrend_dtheta2 <- function(t,...) log(1 + t)  # 2nd
# dtrend_dtheta3 <- function(t,...) rep(0, length(t))
#
# # vol derivatives
# dvol_dtheta1 <- function(t,...) rep(0, length(t))
# dvol_dtheta2 <- function(t,...) rep(0, length(t))
# dvol_dtheta3 <- function(t,...) rep(1, length(t))

# OU with time-dependent rational-function volatility
# alpha_true <- -1
# beta_true <- 3
# sigma_true <- 0.5
# slope_fun <- function(t,alpha,...) rep(alpha, length(t))
# trend_fun <- function(t,beta,...) rep(beta, length(t))
# vol_fun <- function(t,sigma,...) 1/(1 + sigma*t)
# vol_fun <- function(t,sigma,...) t/(1 + sigma*t^2)
#
# # functions with true parameters
# slope_fun_true <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_true <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_true   <- function(t) vol_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
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

dt <- 1
num_steps <- 100
T <- dt * num_steps
num_paths <- 1
x0 <- 1
num_subdivisions <- 100

# Test parameters
alpha_test <- -0.5 # mean reversion
beta_test <- 0.5
sigma_test <- 1
theta_test <- c(alpha_test, beta_test, sigma_test)

# It is of little use to pass slope_fun, ... because inside the function I compute it using the exact transition probability
simulated_path <- simulate_gm_process(
  num_paths, num_steps, dt, x0,
  slope_fun_true, trend_fun_true, vol_fun_true,
  true_alpha, true_beta, true_sigma,
  num_subdivisions
)

# Show captured errors
print(simulated_path$errors)

plot_gm_paths(simulated_path, plot_exact = TRUE)

# Uncomment for Hull-White square-root, Hull-White with logarithmic mean and rational volatility
X <- simulated_path$paths
# X <- simulated_path$paths_exact
# X <- sde.sim(model = "OU", theta = c(beta_true, -alpha_true, sigma_true), N = num_steps, delta = dt, M = 1) # (beta, -alpha, sigma)

# Estimation with likelihoodentero4
opt_result <- tryCatch({
  lk_estimates(num_subdivisions, X, T, theta_test = theta_test,
               slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
}, error = function(e) {
  message("Error in estimation with likelihoodentero4: ", e$message)
  NULL
})
summary(opt_result@coef)

# Estimated values
alpha_est <- opt_result@coef[1]
beta_est <- opt_result@coef[2]
sigma_est <- opt_result@coef[3]

print("optimized theta1 value")
print(alpha_est)
print("optimized theta2 value")
print(beta_est)
print("optimized sigma value")
print(sigma_est)

# Estimates with MLE_grad_loglikehood
grad_result <- tryCatch({
  MLE_grad_loglikelihood(
    initial_guess = c(alpha_test, beta_test, sigma_test),
    T, X, num_subdivisions,
    slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun,
    dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2, dslope_dtheta3 = dslope_dtheta3,
    dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
    dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3,
    method = "Newton", global = "cline"
  )
}, error = function(e) {
  message("Error in estimation with MLE_grad_loglikelihood: ", e$message)
  NULL
})
print(grad_result$x)

# Estimated values
alpha_est <- grad_result$x[1]
beta_est <- grad_result$x[2]
sigma_est <- grad_result$x[3]

print("grad theta1 value")
print(alpha_est)
print("grad theta2 value")
print(beta_est)
print("grad sigma value")
print(sigma_est)