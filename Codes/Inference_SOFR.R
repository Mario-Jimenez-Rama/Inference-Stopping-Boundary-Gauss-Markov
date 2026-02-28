library(readxl)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)

# ------------------------------------------------------------
# 1) Helper: build Tau schedule from Expiry/Tenor/fixed frequency
# ------------------------------------------------------------
build_tau <- function(expiry_years, tenor_years, fixed_freq_per_year) {
  # expiry_years: option expiry T
  # tenor_years: swap tenor (length)
  # fixed_freq_per_year: e.g. 1, 2, 4
  step <- 1 / fixed_freq_per_year
  T0 <- expiry_years
  Tn <- expiry_years + tenor_years
  # include T0 and Tn
  seq(T0, Tn, by = step)
}

build_curves <- function(time, discount, smooth_spar = NULL) {
  stopifnot(length(time) == length(discount))
  
  # Sort by time just in case
  o <- order(time)
  t  <- as.numeric(time[o])
  P  <- as.numeric(discount[o])
  lnP <- log(P)
  
  # ============================
  # SAFETY CHECKS 
  # ============================
  if (any(!is.finite(t)) || any(!is.finite(P)) || any(P <= 0)) {
    stop("Invalid inputs to build_curves(): non-finite or non-positive discount factors.")
  }
  
  if (any(duplicated(t))) {
    stop("Duplicate time points detected in build_curves().")
  }
  
  # ============================
  # 1) METHOD A: splinefun on ln(P)
  # ============================
  lnP_splinefun <- splinefun(t, lnP, method = "fmm")
  
  discount_splinefun <- function(x) {
    exp(lnP_splinefun(x))
  }
  
  inst_forward_splinefun <- function(x) {
    # f(t) = - d/dt ln P(t)
    -lnP_splinefun(x, deriv = 1)
  }
  
  dfdt_splinefun <- function(x) {
    # f'(t) = - d^2/dt^2 ln P(t)
    -lnP_splinefun(x, deriv = 2)
  }
  
  # ============================
  # 2) METHOD B: smooth.spline on ln(P)
  # ============================
  fit_lnP_smooth <- smooth.spline(t, lnP, spar = smooth_spar)
  
  discount_smooth <- function(x) {
    pred <- predict(fit_lnP_smooth, x)
    exp(pred$y)
  }
  
  inst_forward_smooth <- function(x) {
    dlnP_dt <- predict(fit_lnP_smooth, x, deriv = 1)$y
    -dlnP_dt
  }
  
  dfdt_smooth <- function(x) {
    d2lnP_dt2 <- predict(fit_lnP_smooth, x, deriv = 2)$y
    -d2lnP_dt2
  }
  
  # ============================
  # 3) METHOD C: stats::spline on ln(P)
  # ============================
  lnP_spline <- function(x) {
    stats::spline(t, lnP, xout = x, method = "fmm")$y
  }
  
  discount_spline <- function(x) {
    exp(lnP_spline(x))
  }
  
  inst_forward_spline <- function(x) {
    # central difference on ln P
    x <- as.numeric(x)
    h <- 1e-4 * diff(range(t))
    x_clipped <- pmin(max(t) - h, pmax(min(t) + h, x))
    
    dlnP_dt <- (lnP_spline(x_clipped + h) - lnP_spline(x_clipped - h)) / (2 * h)
    -dlnP_dt
  }
  
  dfdt_spline <- function(x) {
    # second derivative of ln P by central difference
    x <- as.numeric(x)
    h <- 1e-4 * diff(range(t))
    x_clipped <- pmin(max(t) - h, pmax(min(t) + h, x))
    
    lnP_plus  <- lnP_spline(x_clipped + h)
    lnP_0     <- lnP_spline(x_clipped)
    lnP_minus <- lnP_spline(x_clipped - h)
    
    d2lnP_dt2 <- (lnP_plus - 2 * lnP_0 + lnP_minus) / (h^2)
    -d2lnP_dt2
  }
  
  if (any(!is.finite(t)) || any(!is.finite(P)) || any(P <= 0)) {
    stop("Invalid inputs to build_curves(): check time grid and discount factors.")
  }
  
  if (any(duplicated(t))) {
    stop("Duplicate time points detected in build_curves().")
  }
  
  
  # Return all functions
  list(
    # METHOD A: splinefun
    discount_splinefun      = discount_splinefun,
    inst_forward_splinefun  = inst_forward_splinefun,
    dfdt_splinefun          = dfdt_splinefun,
    
    # METHOD B: smooth.spline
    discount_smooth         = discount_smooth,
    inst_forward_smooth     = inst_forward_smooth,
    dfdt_smooth             = dfdt_smooth,
    
    # METHOD C: stats::spline
    discount_spline         = discount_spline,
    inst_forward_spline     = inst_forward_spline,
    dfdt_spline             = dfdt_spline
  )
}


# Read SOFR from your Excel (descending dates) and get the time grid

read_sofr_xlsx <- function(path, sheet,
                           start_date,
                           date_col = "Effective Date",
                           rate_col = "Rate (%)",
                           time_axis = c("ACT365", "ACT360")) {
  
  time_axis <- match.arg(time_axis)
  denom <- if (time_axis == "ACT365") 365.0 else 360.0
  
  start_date <- as.Date(start_date)
  stopifnot(!is.na(start_date))
  
  # Read Excel safely
  raw <- read_excel(path, sheet = sheet, col_types = "text")
  
  sofr <- raw %>%
    dplyr::transmute(
      Date = lubridate::mdy(.data[[date_col]]),
      SOFR = as.numeric(.data[[rate_col]])
    ) %>%
    dplyr::filter(!is.na(Date), !is.na(SOFR)) %>%
    dplyr::distinct(Date, .keep_all = TRUE) %>%
    dplyr::arrange(Date)
  
  if (nrow(sofr) < 2) {
    stop("SOFR file has fewer than 2 valid observations.")
  }
  
  # ---- TIME GRID (ONLY FROM REAL SOFR DATES) ----
  sofr <- sofr %>%
    dplyr::filter(Date >= start_date) %>%
    dplyr::mutate(
      t_years = as.numeric(Date - start_date) / denom
    )
  
  if (nrow(sofr) == 0) {
    stop("No SOFR observations on or after start_date.")
  }
  
  list(
    sofr_table      = sofr,              # Date, SOFR, t_years
    time_grid_years = sofr$t_years        # EXACTLY what you asked for
  )
}


# Build realized discount factors from SOFR between start and anchor (fast, vectorized)

build_realized_df_curve_from_sofr <- function(sofr_df, start_date, anchor_date) {
  start_date  <- as.Date(start_date)
  anchor_date <- as.Date(anchor_date)
  stopifnot(start_date < anchor_date)
  
  x <- sofr_df %>%
    filter(Date >= start_date, Date <= anchor_date) %>%
    arrange(Date)
  
  if (nrow(x) < 2) {
    stop("Not enough SOFR observations between start_date and anchor_date.")
  }
  
  # rate r_i applies from Date_i up to Date_{i+1}
  x <- x %>%
    mutate(
      Date_next = lead(Date),
      dt_days   = as.numeric(Date_next - Date),
      tau       = dt_days / 360.0,        # ACT/360 for SOFR accrual
      r         = SOFR / 100.0
    )
  
  # keep accrual steps that end on/before anchor_date
  steps <- x %>%
    filter(!is.na(Date_next), Date_next <= anchor_date) %>%
    mutate(growth = 1 + r * tau)
  
  if (nrow(steps) == 0) stop("No valid accrual steps. Check dates / missing data.")
  
  # cumulative accrual from start -> each step end date
  acc <- cumprod(steps$growth)
  
  # Discount factor P(start, step_end) = 1 / acc
  df_steps <- tibble(
    Date = steps$Date_next,
    DF   = 1 / acc
  )
  
  # Include start point
  df_full <- bind_rows(
    tibble(Date = start_date, DF = 1.0),
    df_steps
  ) %>%
    arrange(Date) %>%
    distinct(Date, .keep_all = TRUE)
  
  df_full
}

# Build curves just from discount factor

build_curves <- function(time, discount, smooth_spar = NULL) {
  stopifnot(length(time) == length(discount))
  
  # Sort by time just in case
  o <- order(time)
  t  <- as.numeric(time[o])
  P  <- as.numeric(discount[o])
  lnP <- log(P)
  
  # ============================
  # 1) METHOD A: splinefun on ln(P)
  # ============================
  lnP_splinefun <- splinefun(t, lnP, method = "fmm")
  
  discount_splinefun <- function(x) {
    exp(lnP_splinefun(x))
  }
  
  inst_forward_splinefun <- function(x) {
    # f(t) = - d/dt ln P(t)
    -lnP_splinefun(x, deriv = 1)
  }
  
  dfdt_splinefun <- function(x) {
    # f'(t) = - d^2/dt^2 ln P(t)
    -lnP_splinefun(x, deriv = 2)
  }
  
  # ============================
  # 2) METHOD B: smooth.spline on ln(P)
  # ============================
  fit_lnP_smooth <- smooth.spline(t, lnP, spar = smooth_spar)
  
  discount_smooth <- function(x) {
    pred <- predict(fit_lnP_smooth, x)
    exp(pred$y)
  }
  
  inst_forward_smooth <- function(x) {
    dlnP_dt <- predict(fit_lnP_smooth, x, deriv = 1)$y
    -dlnP_dt
  }
  
  dfdt_smooth <- function(x) {
    d2lnP_dt2 <- predict(fit_lnP_smooth, x, deriv = 2)$y
    -d2lnP_dt2
  }
  
  # ============================
  # 3) METHOD C: stats::spline on ln(P)
  # ============================
  lnP_spline <- function(x) {
    stats::spline(t, lnP, xout = x, method = "fmm")$y
  }
  
  discount_spline <- function(x) {
    exp(lnP_spline(x))
  }
  
  inst_forward_spline <- function(x) {
    # central difference on ln P
    x <- as.numeric(x)
    h <- 1e-4 * diff(range(t))
    x_clipped <- pmin(max(t) - h, pmax(min(t) + h, x))
    
    dlnP_dt <- (lnP_spline(x_clipped + h) - lnP_spline(x_clipped - h)) / (2 * h)
    -dlnP_dt
  }
  
  dfdt_spline <- function(x) {
    # second derivative of ln P by central difference
    x <- as.numeric(x)
    h <- 1e-4 * diff(range(t))
    x_clipped <- pmin(max(t) - h, pmax(min(t) + h, x))
    
    lnP_plus  <- lnP_spline(x_clipped + h)
    lnP_0     <- lnP_spline(x_clipped)
    lnP_minus <- lnP_spline(x_clipped - h)
    
    d2lnP_dt2 <- (lnP_plus - 2 * lnP_0 + lnP_minus) / (h^2)
    -d2lnP_dt2
  }
  
  # Return all functions
  list(
    # METHOD A: splinefun
    discount_splinefun      = discount_splinefun,
    inst_forward_splinefun  = inst_forward_splinefun,
    dfdt_splinefun          = dfdt_splinefun,
    
    # METHOD B: smooth.spline
    discount_smooth         = discount_smooth,
    inst_forward_smooth     = inst_forward_smooth,
    dfdt_smooth             = dfdt_smooth,
    
    # METHOD C: stats::spline
    discount_spline         = discount_spline,
    inst_forward_spline     = inst_forward_spline,
    dfdt_spline             = dfdt_spline
  )
}


# Splice realized SOFR curve to the market curve at anchor, then build forward function

build_spliced_sofr_to_market_curve <- function(
    sofr_df,
    mkt_curve_df,     # must have Year_Frac (from anchor) and Discount_Factor (P(anchor,T))
    start_date,
    anchor_date,
    time_axis = c("ACT365", "ACT360"),
    smooth_spar = NULL
) {
  time_axis <- match.arg(time_axis)
  
  start_date  <- as.Date(start_date)
  anchor_date <- as.Date(anchor_date)
  stopifnot(start_date < anchor_date)
  
  # choose time axis for t (only for curve parametrization; SOFR accrual uses ACT/360 inside realized DF)
  denom <- ifelse(time_axis == "ACT360", 365.0, 360.0)
  
  # -----------------------------
  # A) realized P(start, d) from SOFR up to anchor
  # -----------------------------
  realized <- build_realized_df_curve_from_sofr(sofr_df, start_date, anchor_date) %>%
    dplyr::mutate(t = as.numeric(Date - start_date) / denom)
  
  # get P(start, anchor)
  DF_start_anchor <- realized %>%
    dplyr::filter(Date == anchor_date) %>%
    dplyr::pull(DF)
  
  if (length(DF_start_anchor) == 0) {
    # If anchor_date isn't exactly a fixing date, approximate using last available < anchor
    DF_start_anchor <- realized %>%
      dplyr::filter(Date < anchor_date) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::pull(DF)
    
    warning("anchor_date not present in SOFR dates; using last available date < anchor_date for DF_start_anchor.")
  }
  
  if (length(DF_start_anchor) != 1 || !is.finite(DF_start_anchor) || DF_start_anchor <= 0) {
    stop("DF_start_anchor is invalid (NA/Inf/non-positive). Check SOFR data and date coverage.")
  }
  
  # time from start to anchor in chosen axis
  t_anchor <- as.numeric(anchor_date - start_date) / denom
  
  # -----------------------------
  # B) market P(anchor, T) beyond anchor; scale to P(start, T)
  #     CRITICAL FIX: do NOT round to calendar dates.
  # -----------------------------
  mkt <- mkt_curve_df %>%
    dplyr::transmute(
      Year_Frac = as.numeric(Year_Frac),
      DF_anchor = as.numeric(Discount_Factor)
    ) %>%
    dplyr::filter(is.finite(Year_Frac), is.finite(DF_anchor)) %>%
    dplyr::arrange(Year_Frac)
  
  # ensure (0,1) exists
  if (nrow(mkt) == 0) stop("mkt_curve_df is empty after cleaning.")
  if (min(mkt$Year_Frac, na.rm = TRUE) > 0) {
    mkt <- dplyr::bind_rows(dplyr::tibble(Year_Frac = 0, DF_anchor = 1), mkt) %>%
      dplyr::arrange(Year_Frac)
  }
  
  # market points expressed on the SAME time axis: t = t_anchor + Year_Frac
  mkt_scaled <- mkt %>%
    dplyr::mutate(
      t  = t_anchor + Year_Frac,
      DF = DF_start_anchor * DF_anchor,
      Date = as.Date(NA)  # keep schema compatible; dates not needed for curve fit
    ) %>%
    dplyr::select(Date, t, DF)
  
  # -----------------------------
  # C) combine, sanitize, and build forward curve
  # -----------------------------
  full_curve <- dplyr::bind_rows(
    realized %>% dplyr::select(Date, t, DF),
    mkt_scaled
  ) %>%
    dplyr::filter(is.finite(t), is.finite(DF), DF > 0) %>%
    dplyr::arrange(t) %>%
    # ensure unique time points (keep the smallest DF if duplicates)
    dplyr::group_by(t) %>%
    dplyr::summarise(
      Date = suppressWarnings(min(Date, na.rm = TRUE)),
      DF   = min(DF),
      .groups = "drop"
    ) %>%
    dplyr::arrange(t)
  
  # If Date was all NA in a group, min(Date, na.rm=TRUE) becomes Inf; fix to NA.
  full_curve <- full_curve %>%
    dplyr::mutate(Date = ifelse(is.infinite(Date), as.Date(NA), as.Date(Date, origin = "1970-01-01")))
  
  # final sanity checks to avoid smooth.spline degeneracy
  if (nrow(full_curve) < 4) stop("Not enough curve points after cleaning (need >= 4).")
  if (diff(range(full_curve$t)) <= 0) stop("Time grid has zero range after cleaning.")
  if (any(duplicated(full_curve$t))) stop("Duplicate t remained after cleaning (unexpected).")
  
  # Build spline functions (your machinery)
  curves <- build_curves(full_curve$t, full_curve$DF, smooth_spar = smooth_spar)
  disc_fun <- curves$discount_smooth
  fwd_fun  <- curves$inst_forward_smooth
  dfdt_fun <- curves$dfdt_spline
  
  # optional: continuous-compounded zero rate z(t) = -log(P)/t
  full_curve <- full_curve %>%
    dplyr::mutate(zero_rate_cc = ifelse(t > 0, -log(DF) / t, NA_real_))
  
  list(
    curve_table      = full_curve,
    time_grid_years  = full_curve$t,    
    DF_start_anchor  = DF_start_anchor,
    disc_fun         = disc_fun,
    fwd_fun          = fwd_fun,
    dfdt_fun         = dfdt_fun,
    curves           = curves
  )
}

# Function to price many financial products using the Hull-White model

hw_pricer <- function(a, sigma, disc_fun, fwd_fun) {
  
  B <- function(t, T) {
    (1 - exp(-a * (T - t))) / a
  }
  
  A <- function(t, T) {
    BtT <- B(t, T)
    disc_fun(T) / disc_fun(t) *
      exp(
        BtT * fwd_fun(t) -
          (sigma^2 / (4 * a)) * (1 - exp(-2 * a * t)) * BtT^2
      )
  }
  
  zcb_price <- function(t, T, r_t) {
    A(t, T) * exp(-B(t, T) * r_t)
  }
  
  sigma_p <- function(T, S) {
    sigma * sqrt((1 - exp(-2 * a * T)) / (2 * a)) * B(T, S)
  }
  
  zcb_put <- function(T, S, K) {
    sp <- sigma_p(T, S)
    PT <- disc_fun(T)
    PS <- disc_fun(S)
    h  <- log(PS / (K * PT)) / sp + 0.5 * sp
    K * PT * pnorm(-h + sp) - PS * pnorm(-h)
  }
  
  zcb_call <- function(T, S, K) {
    sp <- sigma_p(T, S)
    PT <- disc_fun(T)
    PS <- disc_fun(S)
    h  <- log(PS / (K * PT)) / sp + 0.5 * sp
    PS * pnorm(h) - K * PT * pnorm(h - sp)
  }
  
  jamshidian_root <- function(r, T, Tau, K) {
    lhs <- 0
    for (i in 2:length(Tau)) {
      Delta <- Tau[i] - Tau[i-1]
      lhs <- lhs + Delta * K * zcb_price(T, Tau[i], r)
    }
    lhs - (1 - zcb_price(T, Tau[length(Tau)], r))
  }
  
  find_rstar <- function(T, Tau, K) {
    uniroot(
      jamshidian_root,
      interval = c(-50, 50),
      T = T,
      Tau = Tau,
      K = K,
      tol = 1e-12
    )$root
  }
  
  swaption <- function(Tau, N, K, payer = TRUE) {
    T <- Tau[1]
    S <- Tau[length(Tau)]
    r_star <- find_rstar(T, Tau, K)
    
    fixed_leg <- 0
    for (i in 2:length(Tau)) {
      Delta <- Tau[i] - Tau[i-1]
      Ki <- A(T, Tau[i]) * exp(-B(T, Tau[i]) * r_star)
      opt <- if (payer) zcb_put(T, Tau[i], Ki)
      else zcb_call(T, Tau[i], Ki)
      fixed_leg <- fixed_leg + Delta * K * opt
    }
    
    KN <- A(T, S) * exp(-B(T, S) * r_star)
    floating_leg <- if (payer) zcb_put(T, S, KN)
    else zcb_call(T, S, KN)
    
    N * (floating_leg + fixed_leg)
  }
  
  list(
    B = B,
    A = A,
    zcb_put = zcb_put,
    zcb_call = zcb_call,
    swaption = swaption
  )
}



# read sofr BEFORE anchor date

# 3 years
sofr_obj <- read_sofr_xlsx(
  path       = "MarketData/SOFR_RATES_2023_2025.xlsx",
  sheet      = "Results",
  start_date = "2023-01-01",
  time_axis  = "ACT360"
)
sofr_df <- sofr_obj$sofr_table

## 15 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2010_2025.xlsx",
#   sheet      = "Results",
#   start_date = "2023-01-01",
#   time_axis  = "ACT360"
# )
# sofr_df <- sofr_obj$sofr_table

# Read the market curve (your 30/09/2025 discount factors)

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

disc_fun <- res$disc_fun
fwd_fun  <- res$fwd_fun
dfdt_fun <- res$dfdt_fun
full_tbl <- res$curve_table

# For 3 months
# res <- build_curves(mkt_curve_df$Year_Frac, mkt_curve_df$Discount_Factor)
# disc_fun <- res$discount_smooth
# fwd_fun  <- res$inst_forward_smooth
# dfdt_fun <- res$dfdt_smooth

## Plot

t_eval <- seq(min(full_tbl$t), 45, length.out = 400)
plot(t_eval, disc_fun(t_eval), type = "l",
     xlab = "Years from start", ylab = "P(start,t)")

plot(t_eval, fwd_fun(t_eval), type = "l",
     xlab = "Years from start", ylab = "f(t)")
abline(h = 0, lty = 3)

# read sofr AFTER anchor date

# # 3 months
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2025_2026.xlsx",
#   sheet      = "Results",
#   start_date = "2025-09-30",
#   time_axis  = "ACT360"
# )

# 3 years
sofr_obj <- read_sofr_xlsx(
  path       = "MarketData/SOFR_RATES_2023_2026.xlsx",
  sheet      = "Results",
  start_date = "2023-01-01",
  time_axis  = "ACT360"
)

# 15 years
# sofr_obj <- read_sofr_xlsx(
#   path       = "MarketData/SOFR_RATES_2010_2026.xlsx",
#   sheet      = "Results",
#   start_date = "2010-01-01",
#   time_axis  = "ACT360"
# )

sofr_df <- sofr_obj$sofr_table      
time_grid  <- sofr_obj$time_grid_years 

head(sofr_df)
head(time_grid)

## Risk premium estimates

lambda_fun <- function(t_new) {
  B_new <- bs(
    t_new,
    knots = bs_spec$knots,
    Boundary.knots = bs_spec$Boundary.knots,
    degree = bs_spec$degree,
    intercept = bs_spec$intercept
  )
  as.numeric(B_new %*% b_hat)
}


# Hull-White one factor model

# alpha_true <- 0.01000 # mean-reverting
# beta_true  <- 0.5 # risk premium
# sigma_true <- 0.00858 # volatility
# 
# # HW in the physical measure
# 
# slope_fun <- function(t, alpha, ...) rep(-alpha, length(t))
# 
# trend_fun <- function(t, alpha, beta, sigma, ...) {
#   alpha * fwd_fun(t) +
#     dfdt_fun(t) +
#     sigma^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha) + sigma*beta
# }
# 
# vol_fun <- function(t, sigma, ...) rep(sigma, length(t))
# 
# # versions with true parameters (for simulation)
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t   <- function(t) vol_fun(t,   alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivative functions (if you need them later)
# dslope_dtheta1 <- function(t, ...) rep(-1, length(t))
# dslope_dtheta2 <- function(t, ...) rep(0, length(t))
# dslope_dtheta3 <- function(t, ...) rep(0, length(t))
# 
# dtrend_dtheta1 <- function(t, alpha, beta, sigma, ...) {
#   fwd_fun(t) +
#     sigma^2 * exp(-2 * alpha * t) * t / alpha -
#     sigma^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha^2)
# }
# 
# dtrend_dtheta2 <- function(t, sigma, ...) rep(sigma, length(t))
# 
# dtrend_dtheta3 <- function(t, alpha, beta, sigma, ...) {
#   sigma * (1 - exp(-2 * alpha * t)) / alpha + beta
# }
# 
# dvol_dtheta1 <- function(t, ...) rep(0, length(t))
# dvol_dtheta2 <- function(t, ...) rep(0, length(t))
# dvol_dtheta3 <- function(t, ...) rep(1, length(t))
# 
# # Test parameters
# alpha_test <- 0.01 # mean reversion
# beta_test <- 0.005
# sigma_test <- 1
# theta_test <- c(alpha_test, beta_test, sigma_test)

# HW in the risk-free measure

alpha_true <- 0.0111838072 # mean-reverting
beta_true <- 0.0084657428 # volatility
# 
# slope_fun <- function(t, alpha, ...) rep(-alpha, length(t))
# 
# trend_fun <- function(t, alpha, beta, ...) {
#   alpha * fwd_fun(t) +
#     dfdt_fun(t) +
#     beta^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha)
# }
# 
# vol_fun <- function(t, beta, ...) rep(beta, length(t))
# 
# # versions with true parameters (for simulation)
# slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# vol_fun_t   <- function(t) vol_fun(t,   alpha = alpha_true, beta = beta_true, sigma = sigma_true)
# 
# ## derivative functions (if you need them later)
# dslope_dtheta1 <- function(t, ...) rep(-1, length(t))
# dslope_dtheta2 <- function(t, ...) rep(0, length(t))
# dslope_dtheta3 <- function(t, ...) rep(0, length(t))
# 
# dtrend_dtheta1 <- function(t, alpha, beta, ...) {
#   fwd_fun(t) +
#     beta^2 * exp(-2 * alpha * t) * t / alpha -
#     beta^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha^2)
# }
# dtrend_dtheta2 <- function(t, alpha, beta, ...) {
#   beta * (1 - exp(-2 * alpha * t)) / alpha
# }
# 
# dtrend_dtheta3 <- function(t, ...) rep(0, length(t))
# 
# dvol_dtheta1 <- function(t, ...) rep(0, length(t))
# dvol_dtheta2 <- function(t, ...) rep(1, length(t))
# dvol_dtheta3 <- function(t, ...) rep(0, length(t))

# HW in the physical measure, risk premium estimated

slope_fun <- function(t, alpha, ...) rep(-alpha, length(t))

trend_fun <- function(t, alpha, beta, ...) {
  alpha * fwd_fun(t) +
    dfdt_fun(t) +
    beta^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha) + beta*lambda_fun(t)
}

vol_fun <- function(t, beta, ...) rep(beta, length(t))

# versions with true parameters (for simulation)
slope_fun_t <- function(t) slope_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
trend_fun_t <- function(t) trend_fun(t, alpha = alpha_true, beta = beta_true, sigma = sigma_true)
vol_fun_t   <- function(t) vol_fun(t,   alpha = alpha_true, beta = beta_true, sigma = sigma_true)

## derivative functions (if you need them later)
dslope_dtheta1 <- function(t, ...) rep(-1, length(t))
dslope_dtheta2 <- function(t, ...) rep(0, length(t))
dslope_dtheta3 <- function(t, ...) rep(0, length(t))

dtrend_dtheta1 <- function(t, alpha, beta, ...) {
  fwd_fun(t) +
    beta^2 * exp(-2 * alpha * t) * t / alpha -
    beta^2 * (1 - exp(-2 * alpha * t)) / (2 * alpha^2)
}
dtrend_dtheta2 <- function(t, alpha, beta, ...) {
  beta * (1 - exp(-2 * alpha * t)) / alpha
}

dtrend_dtheta3 <- lambda_fun(t)

dvol_dtheta1 <- function(t, ...) rep(0, length(t))
dvol_dtheta2 <- function(t, ...) rep(1, length(t))
dvol_dtheta3 <- function(t, ...) rep(0, length(t))

# Test parameters
alpha_test <- 0.1 # mean reversion
beta_test <- 1
theta_test <- c(alpha_test, beta_test)

## Simulation

# delta <- 1
# n_steps <- 50 
# T <- delta*n_steps
# n_paths <- 1  
# x0 <- 1
# subdivision <- 100
# 
# sim <- simulate_gm_process_nonuni(
#   n_paths = 1,
#   time_grid = time_grid,
#   x0 = 0.05,
#   slope_fun_t = slope_fun_t,
#   trend_fun_t = trend_fun_t,
#   vol_fun_t   = vol_fun_t,
#   alpha_gm = alpha_true,
#   beta_gm  = beta_true,
#   sigma_gm = sigma_true,
#   subdiv = subdivision
# )

# X <- sim$paths[, 1]

## sumulation to see if with longer grids the estimation works better
# build long grid
# n_years = 3
# time_grid_long <- replicate_time_grid_years(time_grid, n_years)
# 
# # 4) simulate on the long grid 
# set.seed(1)
# sim_long <- simulate_gm_process_nonuni(
#   n_paths = 1,
#   time_grid = time_grid_long,
#   x0 = 0.05,
#   slope_fun_t = slope_fun_t,
#   trend_fun_t = trend_fun_t,
#   vol_fun_t   = vol_fun_t,
#   alpha_gm = alpha_true,
#   beta_gm  = beta_true,
#   sigma_gm = sigma_true,
#   subdiv = subdivision
# )

# X_long <- sim_long$paths[, 1]

# Real market data

X <- sofr_obj$sofr_table$SOFR/100.0
subdivision <- 100

# Estimacion con likelihoodentero4, physical measure

# Estimates with MLE_grad_loglikehood
# result <- tryCatch({
#   MLE_grad_loglikelihood(initial_guess = c(alpha_test, beta_test, sigma_test),time_grid, X, subdivision, slope_fun = slope_fun,
#                          trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2,
#                          dslope_dtheta3 = dslope_dtheta3, dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2, dtrend_dtheta3 = dtrend_dtheta3,
#                          dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, dvol_dtheta3 = dvol_dtheta3, method = "Newton", global = "cline")
# }, error = function(e) {
#   message("Error en estimación con MLE_grad_loglikelihood: ", e$message)
#   NULL
# })
# print(result$x)

# opt_result <- tryCatch({
#   lk_estimates_nonuni(subdivision, X, time_grid, theta_test = theta_test,
#                slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
# }, error = function(e) {
#   message("Error en estimación con likelihoodentero4: ", e$message)
#   NULL
# })
# summary(opt_result@coef)

# Estimacion con likelihoodentero4, risk-free measure

opt_result <- tryCatch({
  lk_estimates_2d_nonuni(subdivision, X, time_grid, theta_test = theta_test,
                      slope_fun = slope_fun, trend_fun = trend_fun, vol_fun = vol_fun)
}, error = function(e) {
  message("Error en estimación con likelihoodentero4: ", e$message)
  NULL
})
summary(opt_result@coef)

# result <- tryCatch({
#   MLE_grad_loglikelihood_2d_nonuni(initial_guess = c(alpha_test, beta_test), time_grid, X, subdivision, slope_fun = slope_fun,
#                                    trend_fun = trend_fun, vol_fun = vol_fun, dslope_dtheta1 = dslope_dtheta1, dslope_dtheta2 = dslope_dtheta2,
#                                    dtrend_dtheta1 = dtrend_dtheta1, dtrend_dtheta2 = dtrend_dtheta2,
#                                    dvol_dtheta1 = dvol_dtheta1, dvol_dtheta2 = dvol_dtheta2, method = "Newton", global = "cline")
# }, error = function(e) {
#   message("Error en estimación con MLE_grad_loglikelihood_2d_nonuni: ", e$message)
#   NULL
# })
# print(result$x)
# 
# #Valores estimados
# alpha_est <- result$x[1]
# beta_est <- result$x[2]

#Valores estimados
alpha_est <- opt_result@coef[1]
beta_est <- opt_result@coef[2]
sigma_est <- opt_result@coef[3]

print("valor theta1 optim")
print(alpha_est)
print("valor theta2 optim")
print(beta_est)
print("Valor sigma optim")
print(sigma_est)

#Valores estimados
# alpha_est <- result$x[1]
# beta_est <- result$x[2]
# sigma_est <- result$x[3]
# 
# print("valor theta1 grad")
# print(alpha_est)
# print("valor theta2 grad")
# print(beta_est)
# print("Valor sigma grad")
# print(sigma_est)
  