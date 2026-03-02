import QuantLib as ql
import pandas as pd
import numpy as np
import math

# ========== Inputs ==========
xlsx_path = "MarketData/SWPN_Calibration_Template_30092025_USD.xlsx"
eval_date = ql.Date(30, 9, 2025)  # Sep 30, 2025

# ========== Read data ==========
df_swpn = pd.read_excel(xlsx_path, sheet_name="Template")
df_curve = pd.read_excel(xlsx_path, sheet_name="Curve")

# ===== filter to Expiry <= 10Y =====
# MAX_EXPIRY_Y = 10.0
# df_swpn = df_swpn[df_swpn["Expiry"] <= MAX_EXPIRY_Y].copy()

max_curve_years = df_curve["Year_Frac"].max()
df_swpn = df_swpn[(df_swpn["Expiry"] + df_swpn["Tenor"]) <= max_curve_years - 0.25]

print("Remaining quotes:", len(df_swpn))
print("Max Expiry:", df_swpn["Expiry"].max())
print("Max Tenor:", df_swpn["Tenor"].max())
print("Max Expiry+Tenor:", (df_swpn["Expiry"] + df_swpn["Tenor"]).max())

# ========== Evaluation date ==========
ql.Settings.instance().evaluationDate = eval_date
calendar = ql.UnitedStates(ql.UnitedStates.FederalReserve)

# ========== Build discount curve from year fractions + discount factors ==========
# Curve sheet has: Year_Frac, Discount_Factor
t = df_curve["Year_Frac"].astype(float).to_numpy()
dfs = df_curve["Discount_Factor"].astype(float).to_numpy()

# Convert times -> dates using Actual/365F mapping (consistent and common for building dates)
# curve_day_count = ql.Actual365Fixed()
# curve_dates = [eval_date]
# for ti in t[1:]:
#     days = int(round(ti * 360.0))
#     curve_dates.append(eval_date + days)

# discount_curve = ql.DiscountCurve(curve_dates, list(dfs), curve_day_count, calendar)
# yts = ql.RelinkableYieldTermStructureHandle(discount_curve)

curve_day_count = ql.Actual360()   # <-- match curve's Year_Frac basis (very likely)
curve_dates = [eval_date]
for ti in t[1:]:
    days = int(round(ti * 360.0))
    curve_dates.append(eval_date + days)

discount_curve = ql.DiscountCurve(curve_dates, list(dfs), curve_day_count, calendar)
yts = ql.RelinkableYieldTermStructureHandle(discount_curve)

# SOFR index linked to this curve
sofr = ql.Sofr(yts)

# ========== Swaption quotes preprocessing ==========
# Expiry is given as Act/360-style year fraction, e.g. 1.013889 ~ 1Y
expiries = df_swpn["Expiry"].astype(float).to_numpy()
tenors_y = df_swpn["Tenor"].astype(int).to_numpy()
vol_bps = df_swpn["Volatility (Bps)"].astype(float).to_numpy()

# Convert Expiry -> months tenor that QuantLib SwaptionHelper accepts
expiry_months = np.rint(expiries * (360.0/365.0) * 12.0).astype(int)
expiry_months = np.maximum(expiry_months, 1)

# Normal vol in absolute rate units
normal_vols = vol_bps / 10000.0

# ========== Model + engine ==========
model = ql.HullWhite(yts)  # params: a, sigma
engine = ql.JamshidianSwaptionEngine(model)

fixed_leg_tenor = ql.Period(1, ql.Years)
fixed_leg_daycounter = ql.Thirty360(ql.Thirty360.BondBasis)  # common USD fixed leg
float_leg_daycounter = ql.Actual360()

# ========== Build calibration helpers ==========
helpers = []
for mths, ten, vol in zip(expiry_months, tenors_y, normal_vols):
    vol_quote = ql.SimpleQuote(float(vol))
    vol_handle = ql.QuoteHandle(vol_quote)

    maturity = ql.Period(int(mths), ql.Months)
    length = ql.Period(int(ten), ql.Years)

    h = ql.SwaptionHelper(
        maturity, length, vol_handle, sofr,
        fixed_leg_tenor, fixed_leg_daycounter,
        float_leg_daycounter, yts,
        ql.BlackCalibrationHelper.RelativePriceError,
        ql.nullDouble(), 1.0, ql.Normal
    )
    h.setPricingEngine(engine)
    helpers.append(h)

# ========== Calibrate ==========
opt = ql.LevenbergMarquardt(1e-8, 1e-8, 1e-8)
end = ql.EndCriteria(10000, 100, 1e-8, 1e-8, 1e-8)

model.calibrate(helpers, opt, end)
a, sigma = model.params()
print(f"Calibrated Hull-White params on {eval_date}:")
print(f"  a     = {a:.10f}")
print(f"  sigma = {sigma:.10f}")

# ========== Diagnostics ==========
rows = []
for i, h in enumerate(helpers):
    rows.append({
        "ExpiryMonths": int(expiry_months[i]),
        "TenorYears": int(tenors_y[i]),
        "MarketNormalVol": float(normal_vols[i]),
        "MarketValue": float(h.marketValue()),
        "ModelValue": float(h.modelValue()),
        "RelError": float(h.calibrationError()),
    })

rep = pd.DataFrame(rows)
rmse = math.sqrt((rep["RelError"]**2).mean())
mae = rep["RelError"].abs().mean()

print(f"\nCalibration error (relative price error):")
print(f"  RMSE = {rmse:.6e}")
print(f"  MAE  = {mae:.6e}")


