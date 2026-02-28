library(readxl)
library(dplyr)
library(lubridate)

# -----------------------------
# 1) ACM data (monthly, month-end)
# -----------------------------
acm_df <- read_excel("MarketData/ACMTermPremium.xls") %>%
  mutate(
    DATE = as.Date(DATE, format = "%d-%b-%Y"),
    DATE = ceiling_date(DATE, "month") - days(1)   # force month-end
  ) %>%
  filter(DATE >= as.Date("2023-01-01"),
         DATE <= as.Date("2026-01-31"))            # monthly -> allow full month

# columns to correlate (everything except DATE)
acm_cols <- setdiff(names(acm_df), "DATE")

# -----------------------------
# 2) lambda(t) -> monthly mean at month-end
# -----------------------------
lambda_monthly <- data.frame(
  DATE = as.Date(sofr_obj$sofr_table$Date),
  lambda = lambda_hat
) %>%
  mutate(DATE = ceiling_date(DATE, "month") - days(1)) %>%  # month-end stamp
  group_by(DATE) %>%
  summarise(lambda_m = mean(lambda, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 3) Join and compute correlations
# -----------------------------
compare_df <- inner_join(lambda_monthly, acm_df, by = "DATE")

# Correlation lambda_m with every ACM column
cors <- sapply(acm_cols, function(v) {
  cor(compare_df$lambda_m, compare_df[[v]], use = "complete.obs")
})

# Put into a nice table
cor_tbl <- tibble(
  variable = names(cors),
  correlation = as.numeric(cors),
  n_pairs = sapply(acm_cols, function(v) sum(complete.cases(compare_df$lambda_m, compare_df[[v]])))
) %>%
  arrange(desc(abs(correlation)))

print(cor_tbl, n = Inf)

pvals <- sapply(acm_cols, function(v) {
  x <- compare_df$lambda_m
  y <- compare_df[[v]]
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor.test(x[ok], y[ok])$p.value
})

cor_tbl <- cor_tbl %>%
  mutate(p_value = pvals[variable])

print(cor_tbl, n = Inf)

top_vars <- head(cor_tbl$variable, 2)

for (v in top_vars) {
  plot(compare_df$DATE, compare_df$lambda_m, type="l",
       xlab="Date", ylab="Value", main=paste("lambda_m vs", v))
  lines(compare_df$DATE, scale(compare_df[[v]]), col="red")
  legend("topleft", legend=c("lambda_m", paste0("scaled ", v)),
         col=c("black","red"), lwd=2, bty="n")
}
