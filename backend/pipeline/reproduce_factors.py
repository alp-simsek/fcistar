#!/usr/bin/env python3
"""
reproduce_factors.py — validation: reproduce the official FCI-G sub-indices from our public
proxies + the Fed's own multipliers, to confirm the construction and our data transforms.

Two outputs:
  data/weights/weights_official.csv   the Fed multipliers (refs/multipliers.csv) with named
                                      columns + lag; baseline 3yr index uses rows lag=1..12.
  validation/reproduction_*.{csv,json}  per-variable reproduced vs published, with fit stats.

Construction (from the Fed's R replication code, refs/public_fci_publicrelease.R):
  Factor_j(t) = - sum_{k=0..11} multiplier[k, j] * dX_j(t - 3k)
where dX_j is the QUARTERLY (3-month) change of driver j, and the 12 lags are quarter-spaced
(t, t-3, ..., t-33 months) — i.e. 12 quarterly changes sampled monthly. The trailing -1 gives
the published sign. Per the Technical Appendix, the driver level entering the 3-month change is:
  rates (FFR/10Y/mortgage/BBB): 3-month average of daily values, differenced  -> pp
  dollar:                       100*log of the 3-month average daily level, differenced
  stock/house:                  100*log of the end-of-period level, differenced

We validate the three variables that drive most of the index: STOCK, DOLLAR, BBB.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "data" / "raw"
WEIGHTS_DIR = HERE / "data" / "weights"
VALID_DIR = HERE / "validation"
OUTPUT_DIR = HERE / "data" / "output"
OFFICIAL = HERE / "data" / "fcig_official" / "fci_g_public_monthly_3yr.csv"
MULTIPLIERS = HERE / "refs" / "multipliers.csv"

# Column order of multipliers.csv / input_data.csv (Fed's code-compatible order). Names match the
# official published sub-index columns for readability.
VAR_ORDER = ["FFR", "10Yr Treasury", "Mortgage Rate", "BBB", "Stock Market", "House Prices", "Dollar"]

N_LAGS = 12  # baseline 3-year index uses the first 12 (quarter-spaced) multiplier rows


def build_weights_official() -> pd.DataFrame:
    """refs/multipliers.csv -> data/weights/weights_official.csv (named, with lag)."""
    m = pd.read_csv(MULTIPLIERS)            # columns "1".."7", 20 rows
    m.columns = VAR_ORDER
    m.insert(0, "lag", range(1, len(m) + 1))
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    m.to_csv(WEIGHTS_DIR / "weights_official.csv", index=False)
    (WEIGHTS_DIR / "weights_metadata.json").write_text(json.dumps({
        "source": "Fed FCI-G replication ZIP (public_fci_publicrelease.zip) -> multipliers.csv; "
                  "ZIP linked from the FEDS Note as 'FCI-G Code (ZIP)'.",
        "shape": list(m.shape),
        "var_order": VAR_ORDER,
        "convention": "Raw (beta_4 - beta_0) GDP-growth multipliers, exactly as the Fed ships "
                      "them. Used by the R code as: published factor = -1 * sum(multiplier * dX), "
                      "where stock & house drivers enter as the NEGATIVE 3-month log change.",
        "vs_table_1": "Magnitudes are identical to Table 1 of the note. Table 1 re-orients every "
                      "weight so tighter conditions (higher rates/dollar, lower stock/house) raise "
                      "FCI-G: Table1 = -1*multiplier for FFR/10Y/Mortgage/BBB/Dollar, and "
                      "Table1 = +1*multiplier for Stock/House.",
        "lags": "20 quarters (full IRF). Baseline 3yr index uses lag 1..12; 1yr variant lag 1..4. "
                "Beyond ~3 years all weights are 0 EXCEPT House (GDP effects persist ~5 years), but "
                "House's lag 13..20 are unused by the baseline anyway.",
    }, indent=2))
    return m


# --------------------------------------------------------------------------------------
# Monthly driver transforms
# --------------------------------------------------------------------------------------

def _load_daily(name: str) -> pd.Series:
    df = pd.read_csv(RAW / f"{name}.csv", parse_dates=["date"])
    return df.set_index("date")["value"].sort_index()


def _monthly_index(s: pd.Series) -> pd.PeriodIndex:
    return pd.period_range(s.index.to_period("M").min(), s.index.to_period("M").max(), freq="M")


def three_month_avg(s: pd.Series) -> pd.Series:
    """Exact daily average over the trailing 3 calendar months, indexed by month (Period)."""
    p = s.index.to_period("M")
    msum = s.groupby(p).sum()
    mcnt = s.groupby(p).count()
    full = pd.period_range(msum.index.min(), msum.index.max(), freq="M")
    msum, mcnt = msum.reindex(full), mcnt.reindex(full)
    return (msum.rolling(3).sum() / mcnt.rolling(3).sum()).rename("A")


def eom_level(s: pd.Series) -> pd.Series:
    """End-of-month level, indexed by month (Period), on a continuous monthly grid."""
    p = s.index.to_period("M")
    last = s.groupby(p).last()
    return last.reindex(pd.period_range(p.min(), p.max(), freq="M"))


def driver_level(name: str, transform: str) -> pd.Series:
    """Monthly transformed level X(t) per the appendix convention."""
    s = _load_daily(name)
    if transform == "avg_level":          # rates: 3mo daily average, in pp
        return three_month_avg(s)
    if transform == "avg_log":            # dollar: 100*log of 3mo daily average
        return 100.0 * np.log(three_month_avg(s))
    if transform == "eom_log":            # dollar-like price: +100*log of end-of-month level
        return 100.0 * np.log(eom_level(s))
    if transform == "eom_neglog":         # stock/house: oriented so increase = tightening
        return -100.0 * np.log(eom_level(s))
    raise ValueError(transform)


def reproduce_factor(level: pd.Series, mult_col: pd.Series) -> pd.Series:
    """Factor(t) = - sum_{k=0..11} mult[k] * dX(t-3k), dX = 3-month change of level."""
    dX = level - level.shift(3)
    out = pd.Series(0.0, index=level.index)
    for k in range(N_LAGS):
        out = out + float(mult_col.iloc[k]) * dX.shift(3 * k)
    return -out


def build_predicted_bbb(weights: pd.DataFrame, official: pd.DataFrame) -> dict:
    """Write data/raw/predicted_bbb.csv: a long-history BBB yield derived from Moody's Baa,
    `predicted_bbb = alpha + beta * Baa`, used as the BBB INPUT for backcasting.
      beta : scales Baa's variation to the ICE BBB index, so the BBB FACTOR built from
             predicted_bbb reproduces the published factor (beta ~ 1.34, from regressing the
             published BBB factor on the Baa-reconstructed factor over the full sample).
      alpha: aligns the level to ICE BBB over the 2023+ overlap (cosmetic — only changes matter).
    Stored in raw/ (per its role as a backcasting input); regenerated by this script, not fetched.
    """
    baa = _load_daily("dbaa")
    factor_baa = reproduce_factor(driver_level("dbaa", "avg_level"), weights["BBB"]).reindex(official.index)
    df = pd.concat([official["BBB"].rename("p"), factor_baa.rename("r")], axis=1).dropna()
    beta = float(np.polyfit(df["r"].to_numpy(), df["p"].to_numpy(), 1)[0])
    ice = _load_daily("bbb")
    ov = pd.concat([ice.rename("i"), baa.rename("b")], axis=1).dropna()
    alpha = float(ov["i"].mean() - beta * ov["b"].mean())
    pred = (alpha + beta * baa).rename("value").rename_axis("date").reset_index()
    pred["date"] = pred["date"].dt.strftime("%Y-%m-%d")
    pred.to_csv(RAW / "predicted_bbb.csv", index=False)
    return {"beta": round(beta, 4), "alpha": round(alpha, 4), "overlap_n": int(len(ov))}


# --------------------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------------------

VARS = [
    # name,      proxy file,      transform,    multiplier col,   published column
    ("ffr",      "dff",          "avg_level",  "FFR",           "FFR"),
    ("t10yr",    "svenpy10",     "avg_level",  "10Yr Treasury", "10Yr Treasury"),
    ("mortgage", "mortgage30us", "avg_level",  "Mortgage Rate", "Mortgage Rate"),
    # BBB backcast uses `predicted_bbb` (= alpha + beta*Baa, ICE-calibrated, built by
    # build_predicted_bbb), fed through the SAME construction. The live nowcast instead uses ICE
    # `bbb` (2023+, exact) anchored to the release.
    ("bbb",      "predicted_bbb","avg_level",  "BBB",           "BBB"),
    ("stock",    "wilshire5000", "eom_neglog", "Stock Market",  "Stock Market"),
    ("house",    "zillow_zhvi",  "eom_neglog", "House Prices",  "House Prices"),
    ("dollar",   "dollar_broad", "avg_log",    "Dollar",        "Dollar"),
]


def fit_stats(pub: pd.Series, rep: pd.Series) -> dict:
    df = pd.concat([pub.rename("pub"), rep.rename("rep")], axis=1).dropna()
    if len(df) < 6:
        return {"n": int(len(df)), "note": "insufficient overlap to validate"}
    pub_, rep_ = df["pub"].to_numpy(), df["rep"].to_numpy()
    corr = float(np.corrcoef(pub_, rep_)[0, 1])
    b, a = np.polyfit(rep_, pub_, 1)            # pub = a + b*rep
    resid = pub_ - (a + b * rep_)
    return {
        "n": int(len(df)),
        "sample": f"{df.index.min()}..{df.index.max()}",
        "corr": round(corr, 5),
        "r2": round(corr ** 2, 5),
        "rmse_raw": round(float(np.sqrt(np.mean((pub_ - rep_) ** 2))), 5),
        "slope": round(float(b), 4),
        "intercept": round(float(a), 5),
        "rmse_after_fit": round(float(np.sqrt(np.mean(resid ** 2))), 5),
    }


def main() -> None:
    weights = build_weights_official()
    print(f"weights_official.csv written ({weights.shape[0]}x{weights.shape[1]-1}).\n")

    official = pd.read_csv(OFFICIAL, parse_dates=["date"])
    official["m"] = official["date"].dt.to_period("M")
    official = official.set_index("m")

    bbb_cal = build_predicted_bbb(weights, official)   # writes data/raw/predicted_bbb.csv
    print(f"predicted_bbb.csv written: beta={bbb_cal['beta']} alpha={bbb_cal['alpha']} "
          f"(ICE/Baa overlap n={bbb_cal['overlap_n']}).\n")

    out = pd.DataFrame(index=official.index)
    stats = {}
    for name, proxy, transform, mcol, pubcol in VARS:
        level = driver_level(proxy, transform)
        rep = reproduce_factor(level, weights[mcol]).reindex(official.index)
        pub = official[pubcol]
        out[f"{name}_published"] = pub
        out[f"{name}_reproduced"] = rep
        stats[name] = {"proxy": proxy, "transform": transform, **fit_stats(pub, rep)}

    # Capstone: total index = sum of reproduced factors (BBB via the Baa regression), falling back
    # to the published factor only for early months before a variable's 3yr window fills. This is
    # the backcast: the whole index rebuilt from public proxies alone.
    total_pub = official["FCI-G Index (baseline)"]
    total_hat = pd.Series(0.0, index=official.index)
    for name, *_ in VARS:
        total_hat = total_hat + out[f"{name}_reproduced"].fillna(out[f"{name}_published"])
    out["index_published"] = total_pub
    out["index_reproduced"] = total_hat
    stats["index"] = {"proxy": "7 proxies (BBB via Baa regression)", **fit_stats(total_pub, total_hat)}

    VALID_DIR.mkdir(parents=True, exist_ok=True)
    out.index = out.index.to_timestamp("M").strftime("%Y-%m-%d")
    out.to_csv(VALID_DIR / "reproduction_factors.csv")
    (VALID_DIR / "reproduction_stats.json").write_text(json.dumps(stats, indent=2))

    # Standalone BACKCAST: the index + 7 factor contributions rebuilt from public proxies (BBB via
    # the Baa regression). The derived counterpart to the raw proxies in data/raw/. Where a proxy
    # doesn't reach back far enough (house pre-2003, dollar pre-2009), that factor falls back to the
    # official published contribution, so fci_g is exactly the row sum of the 7 factor columns.
    backcast = pd.DataFrame({"date": out.index})
    backcast["fci_g"] = out["index_reproduced"].values
    for name, *_ in VARS:
        backcast[name] = out[f"{name}_reproduced"].fillna(out[f"{name}_published"]).values
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    backcast.to_csv(OUTPUT_DIR / "fcig_backcast.csv", index=False)

    print(f"{'var':9s} {'proxy':13s} {'n':>4s} {'corr':>7s} {'r2':>7s} "
          f"{'rmse_raw':>9s} {'slope':>6s} {'rmse_fit':>9s}   sample")
    for name, st in stats.items():
        if "note" in st:
            print(f"{name:9s} {st['proxy']:13s} {st['n']:4d}   -- {st['note']} --")
            continue
        print(f"{name:9s} {st['proxy']:13s} {st['n']:4d} {st['corr']:7.4f} {st['r2']:7.4f} "
              f"{st['rmse_raw']:9.4f} {st['slope']:6.3f} {st['rmse_after_fit']:9.4f}   {st['sample']}")


if __name__ == "__main__":
    main()
