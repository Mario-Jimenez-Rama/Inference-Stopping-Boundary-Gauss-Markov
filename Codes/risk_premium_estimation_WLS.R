library(splines)

# read sofr BEFORE anchor date

# 3 years
sofr_obj <- read_sofr_xlsx(
  path       = "MarketData/SOFR_RATES_2023_2025.xlsx",
  sheet      = "Results",
  start_date = "2023-01-01",
  time_axis  = "ACT360"
)

# # 1 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_Feb_2025_2025.xlsx",
#   sheet      = "Results",
#   start_date = "2025-02-05",
#   time_axis  = "ACT360"
# )

## 15 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2010_2025.xlsx",
#   sheet      = "Results",
#   start_date = "2023-01-01",
#   time_axis  = "ACT360"
# )

sofr_df <- sofr_obj$sofr_table

curve_xlsx <- read_excel("MarketData/SWPN_Calibration_Template_30092025_USD.xlsx",
                         sheet = "Curve")

mkt_curve_df <- curve_xlsx %>%
  transmute(
    Year_Frac = as.numeric(Year_Frac),
    Discount_Factor = as.numeric(Discount_Factor)
  )

# Build the spliced curve and functions

# 15 years
# res <- build_spliced_sofr_to_market_curve(
#   sofr_df      = sofr_df,
#   mkt_curve_df = mkt_curve_df,
#   start_date   = "2010-01-01",
#   anchor_date  = "2025-09-30",
#   time_axis    = "ACT360",
#   smooth_spar  = NULL
# )

# 3 years
res <- build_spliced_sofr_to_market_curve(
  sofr_df      = sofr_df,
  mkt_curve_df = mkt_curve_df,
  start_date   = "2023-01-01",
  anchor_date  = "2025-09-30",
  time_axis    = "ACT360",
  smooth_spar  = NULL
)

# # 1 years
# res <- build_spliced_sofr_to_market_curve(
#   sofr_df      = sofr_df,
#   mkt_curve_df = mkt_curve_df,
#   start_date   = "2025-02-05",
#   anchor_date  = "2025-09-30",
#   time_axis    = "ACT360",
#   smooth_spar  = NULL
# )

disc_fun <- res$disc_fun
fwd_fun  <- res$fwd_fun
dfdt_fun <- res$dfdt_fun
full_tbl <- res$curve_table

# For 3 months
# res <- build_curves(mkt_curve_df$Year_Frac, mkt_curve_df$Discount_Factor)
# disc_fun <- res$discount_smooth
# fwd_fun  <- res$inst_forward_smooth
# dfdt_fun <- res$dfdt_smooth

# # 3 years
sofr_obj <- read_sofr_xlsx(
  path       = "MarketData/SOFR_RATES_2023_2026.xlsx",
  sheet      = "Results",
  start_date = "2023-01-01",
  time_axis  = "ACT360"
)

# # 1 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_Feb_2025_2026.xlsx",
#   sheet      = "Results",
#   start_date = "2025-02-05",
#   time_axis  = "ACT360"
# )

# 15 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2010_2026.xlsx",
#   sheet      = "Results",
#   start_date = "2010-01-01",
#   time_axis  = "ACT360"
# )

# 3 months
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2025_2026.xlsx",
#   sheet      = "Results",
#   start_date = "2025-09-30",
#   time_axis  = "ACT360"
# )

r <- sofr_obj$sofr_table$SOFR / 100  # SOFR rates      
t  <- sofr_obj$time_grid_years       # time grid

a <- 0.0111838072 # mean-reverting
sigma <- 0.0084657428 # volatility

dt <- diff(t)
dr <- diff(r)

ti <- t[-length(t)]
ri <- r[-length(r)]

theta_i <- a * fwd_fun(ti) +
  dfdt_fun(ti) +
  sigma^2 * (1 - exp(-2 * a * ti)) / (2 * a) 

k <- 3  # drop first k increments
dt <- dt[-(1:k)]
dr <- dr[-(1:k)]
ti <- ti[-(1:k)]
ri <- ri[-(1:k)]
theta_i <- theta_i[-(1:k)]

y <- dr - (theta_i - a * ri) * dt   # residual increment

# Spline definition 1

df_spline <- 5
# cubic B spline
B  <- bs(ti, df = df_spline, degree = 3, intercept = TRUE)
# natural spline
# B <- ns(ti, df=df_spline, intercept=TRUE)

X  <- sigma * (dt * B)     # y ≈ X b + error
w  <- 1/dt                 # since Var(error_i) = sigma^2 * dt_i

fit <- lm(y ~ X - 1, weights = w)

# linear WLS

# Using the "standardized" form:
# z  <- y / (sigma * sqrt(dt))
# X2 <- cbind(1*sqrt(dt), ti*sqrt(dt))  # alpha + beta t
# 
# fit_lin <- lm(z ~ X2 - 1)  # OLS
# 
# cat("coefficent: ",coef(fit_lin))     # should be comparable to MLE if assumptions match

b_hat <- coef(fit)

Vb <- vcov(fit)            # covariance matrix of b_hat (already accounts for weights)

# compute market price risk, Splines

# cubic B spline
B_full <- bs(t, df = df_spline, degree = 3, intercept = TRUE)
# natural spline
# B_full <- ns(t, df=df_spline, intercept=TRUE)

lambda_hat <- as.numeric(B_full %*% b_hat)

# pointwise variance: diag(B Vb B')
lambda_var <- rowSums((B_full %*% Vb) * B_full)
lambda_se  <- sqrt(pmax(lambda_var, 0))


alpha <- 0.1
zcrit <- qnorm(1 - alpha/2)

lower <- lambda_hat - zcrit * lambda_se
upper <- lambda_hat + zcrit * lambda_se

# Save splines
# B-spline
df_spline <- 5
B <- bs(ti, df = df_spline, degree = 3, intercept = TRUE)

# save basis "spec"
bs_spec <- attributes(B)
# bs_spec$knots
# bs_spec$Boundary.knots
# bs_spec$degree
# bs_spec$intercept
# print("spline coefficients if lambda ")
# print(lambda_hat)
# print("lower confidence pointwise spline coefficients if lambda ")
# print(lower)
# print("upper confidence pointwise spline coefficients if lambda ")
# print(upper)

# Plot lanbda and the confidence curves


# plot theta 

plot(ti, theta_i, type="l")

# Define colors
col_wls   <- "#1F77B4"   # blue
col_wls_f <- adjustcolor(col_wls, alpha.f = 0.25)  # lighter fill
col_lin   <- "#D62728"   # red

plot(t, lambda_hat, type="n",
     xlab="t (years)", ylab=expression(lambda(t)),
     ylim = c(-5, 5))

# Confidence band (lighter shade of WLS color)
polygon(c(t, rev(t)),
        c(lower, rev(upper)),
        border = NA, col = col_wls_f)

# WLS spline estimate
lines(t, lambda_hat, col = col_wls, lwd = 2)

# Linear MLE
lambda_lin <- -0.2236956 + 0.1788985 * t
lines(t, lambda_lin, col = col_lin, lwd = 2)