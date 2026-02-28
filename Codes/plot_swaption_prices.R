library(purrr)
library(tidyr)

parse_payment_dates <- function(x) {
  if (is.na(x)) return(numeric(0))
  x <- trimws(x)
  x <- gsub("^\\[|\\]$", "", x)            # remove leading [ and trailing ]
  x <- trimws(x)
  if (nchar(x) == 0) return(numeric(0))
  as.numeric(strsplit(x, "\\s*,\\s*")[[1]])
}

# ------------------------------------------------------------
# 3) MAIN: read Excel and price
# ------------------------------------------------------------

# ---- YOU SET THESE ----

# Hull–White params from calibration (Q) and from MLE
a_calib <- 0.0111838072
sigma_calib <- 0.0084657428

a_mle <- 0.0111838072
#sigma_mle <- beta_est
sigma_mle <- 0.01206489
# load swaption sheet
xlsx_path <- "MarketData/SWPN_Calibration_Template_30092025_USD.xlsx" 

df <- read_excel(xlsx_path, sheet = "Template")

# clean + compute model prices
df <- df %>%
  mutate(
    Payment_Dates = lapply(Payment_Dates, parse_payment_dates),
    Expiry = as.numeric(Expiry),
    Tenor  = as.numeric(Tenor),
    Strike = as.numeric(Strike),
    Price  = as.numeric(Price),
    Notional = as.numeric(Notional),
    payer = tolower(Opt_Type) == "payer",
    K_dec = ifelse(Strike > 1, Strike/100, Strike)   # strike in decimals for pricing
  )

# create again the curves of discount and forward
curve <- read_excel(xlsx_path, sheet = "Curve") %>%
  transmute(
    t = as.numeric(Year_Frac),
    DF = as.numeric(Discount_Factor)
  ) %>% arrange(t)

curves_mkt <- build_curves(curve$t, curve$DF, smooth_spar = NULL)

disc_fun <- curves_mkt$discount_smooth
fwd_fun  <- curves_mkt$inst_forward_smooth

max_curve <- max(curve$t, na.rm = TRUE)

price_one <- function(pricer, Tau, N, K, payer, Expiry, max_curve) {
  # domain checks
  if (any(!is.finite(Tau))) return(NA_real_)
  if (max(Tau) > max_curve) return(NA_real_)   # outside curve
  if (!is.finite(K) || K <= 0) return(NA_real_)
  if (!is.finite(Expiry) || Expiry < 0) return(NA_real_)
  
  pv <- tryCatch(pricer$swaption(Tau, N = N, K = K, payer = payer),
                 error = function(e) NA_real_)
  if (is.na(pv)) return(NA_real_)
  
  # Convert PV to "forward premium" like your Python: PV / DF(Expiry)
  dfT <- disc_fun(Expiry)
  if (!is.finite(dfT) || dfT <= 0) return(NA_real_)
  pv / dfT
}

pr_calib <- hw_pricer(a_calib, sigma_calib, disc_fun, fwd_fun)
pr_mle   <- hw_pricer(a_mle,   sigma_mle,   disc_fun, fwd_fun)

df <- df %>%
  mutate(
    Price_Calib = pmap_dbl(list(Payment_Dates, Notional, K_dec, payer, Expiry),
                           \(Tau, N, K, pay, T) price_one(pr_calib, Tau, N, K, pay, T, max_curve)
    ),
    Price_MLE = pmap_dbl(list(Payment_Dates, Notional, K_dec, payer, Expiry),
                         \(Tau, N, K, pay, T) price_one(pr_mle, Tau, N, K, pay, T, max_curve)
    )
  )

summary(df$Price_Calib)
summary(df$Price_MLE)
sum(is.na(df$Price_Calib))
sum(is.na(df$Price_MLE))

# ------------------------------------------------------------
# 4) Plot: ATM term structure by tenor
# ------------------------------------------------------------

# keep only ATM rows (adjust depending on your file conventions)
atm <- df %>%
  dplyr::filter(Moneyness == "ATM") %>%
  dplyr::filter(Tenor %in% c(5, 10, 20, 30)) %>%
  dplyr::select(Expiry, Tenor, Price, Price_Calib, Price_MLE) %>%
  tidyr::pivot_longer(
    cols = c(Price, Price_Calib, Price_MLE),
    names_to = "Source",
    values_to = "SwaptionPrice"
  ) %>%
  dplyr::mutate(Source = dplyr::recode(
    Source,
    Price = "Market",
    Price_Calib = "Model (Calib)",
    Price_MLE = "Model (MLE)"
  ))

ggplot(atm, aes(x = Expiry, y = SwaptionPrice, color = Source, linetype = Source)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Tenor, scales = "free_y", ncol = 2,
             labeller = labeller(Tenor = function(x) paste0("ATM, ", x, "Y Tenor"))) +
  labs(
    x = "Expiry (Years)",
    y = "Price"
  ) +
  scale_color_manual(
    name = "",
    values = c("Market" = "blue", "Model (Calib)" = "red", "Model (MLE)" = "darkgreen")
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Market" = "solid", "Model (Calib)" = "dashed", "Model (MLE)" = "dashed")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",  # You can change to "bottom", "left", "right"
    legend.key.size = unit(0.5, "cm"),  # Small symbol size
    legend.text = element_text(size = 8),  # Small text
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0),  # Reduce margin
    legend.box.spacing = unit(0.2, "cm")  # Small spacing
  )

