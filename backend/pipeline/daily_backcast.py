#!/usr/bin/env python3
"""
daily_backcast.py — DAILY pseudo-real-time backcast of FCI-G (= the nowcast engine run at every
historical date), using TRAILING 3-month windows measured from each date d (not snapped to the
calendar month). This rolls the lag structure continuously, so there is no month-boundary jump,
and at month-ends the trailing windows coincide with the official calendar windows (so it still
matches the published monthly value).

At date d, for each factor j:
  Factor_j(d) = - sum_{k=0..11} mult[k] * ( A_j(d - 3k mo) - A_j(d - (3k+3) mo) )
where A_j is the driver's transformed level evaluated with a trailing window ending at the date:
  - rates  (avg_level): A_j(t) = mean of daily values over (t-3mo, t]                  [pp]
  - dollar (avg_log)  : A_j(t) = 100*log( mean daily level over (t-3mo, t] )
  - stock/house (eom) : A_j(t) = -100*log( level at t )   (end-of-period, spot; signed)

Everything uses only data through d (causal). Two factors splice the EXACT recent Fed source onto a
long-history backcast proxy: BBB = ICE `bbb` (2023+) over Baa-derived `predicted_bbb`; mortgage =
Optimal Blue `obmmi30` (2017+) over Freddie `mortgage30us`. 10Y = `svenpy10` (GSW), with its ~1-week
lagging tail extended by `dgs10` changes. Purpose: (i) fill daily gaps where official (monthly)
FCI-G doesn't exist; (ii) test the nowcast.

Output: data/output/fcig_daily_backcast.csv  (date, fci_g, + 7 factor contributions).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from reproduce_factors import (
    VARS, OUTPUT_DIR, OFFICIAL,
    build_weights_official, build_predicted_bbb, _load_daily,
)

N_LAGS = 12  # baseline 3-year lookback: 12 quarterly changes


def trailing_avg_daily(s: pd.Series, full: pd.DatetimeIndex) -> pd.Series:
    """Trailing 3-calendar-month average of daily s, evaluated for every day in `full`."""
    sf = s.reindex(full)
    csum = sf.fillna(0.0).cumsum()
    ccnt = sf.notna().astype(int).cumsum()
    sh = full - pd.DateOffset(months=3)                       # window start (exclusive)
    S = csum.to_numpy() - np.nan_to_num(csum.reindex(sh).to_numpy())
    C = ccnt.to_numpy() - np.nan_to_num(ccnt.reindex(sh).to_numpy())
    return pd.Series(np.where(C > 0, S / np.where(C > 0, C, 1), np.nan), index=full)


# Two factors splice a long-history backcast proxy with the EXACT recent Fed source:
#   BBB      : Baa-derived predicted_bbb (history)  +  ICE bbb (2023+, exact)
#   mortgage : Freddie mortgage30us (history)       +  Optimal Blue obmmi30 (2017+, exact)
SPLICE = {"predicted_bbb": ("predicted_bbb", "bbb"), "mortgage30us": ("mortgage30us", "obmmi30")}
# 10Y reference svenpy10 (GSW) lags ~1 week; extend its TAIL with dgs10 changes for the nowcast.
BRIDGE = {"svenpy10": "dgs10"}


def _splice(history: pd.Series, recent: pd.Series) -> pd.Series:
    """Join a long-history proxy to the exact recent source, anchored at the handoff.

    Steps: (1) handoff date sp = first date `recent` exists; (2) level gap
    offset = recent(sp) - history(sp); (3) return history(before sp) + offset, then recent(from sp).
    The two sources sit at different levels (Freddie survey vs Optimal Blue lock; Baa vs ICE BBB),
    so the constant `offset` makes them meet exactly at sp. Since the index uses only 3-MONTH
    CHANGES (never the level), this shift cancels in every difference within the history portion —
    it only affects the single 3-month change that straddles sp, removing the spurious inter-source
    level jump there. The spliced series' absolute level is otherwise irrelevant.
    """
    sp = recent.index.min()
    offset = float(recent.loc[sp] - history.asof(sp))
    return pd.concat([history[history.index < sp] + offset, recent]).sort_index()


def _bridge_tail(ref: pd.Series, fill: pd.Series) -> pd.Series:
    """Extend `ref` past its last date using `fill`'s daily CHANGES (same anchor-on-level,
    extend-by-changes idea as _splice, applied at the recent edge): for dates after ref's last
    print, ref_last + (fill(d) - fill(at ref's last date)). Used for the 10Y nowcast tail
    (svenpy10 lags ~1wk; dgs10 is current). Only the few trailing days are affected."""
    last = ref.index.max()
    tail = fill.index[fill.index > last]
    if len(tail) == 0:
        return ref
    ext = float(ref.loc[last]) + (fill.loc[tail] - float(fill.asof(last)))
    return pd.concat([ref, ext]).sort_index()


def daily_input(proxy: str) -> pd.Series:
    """The daily driver series for a factor, applying splices/bridges where defined."""
    if proxy in SPLICE:
        hist, rec = SPLICE[proxy]
        return _splice(_load_daily(hist), _load_daily(rec))
    if proxy in BRIDGE:
        return _bridge_tail(_load_daily(proxy), _load_daily(BRIDGE[proxy]))
    return _load_daily(proxy)


def trailing_level(proxy: str, transform: str, full: pd.DatetimeIndex) -> pd.Series:
    """Driver's transformed level A_j(t) on a daily calendar, per the appendix convention."""
    s = daily_input(proxy)
    if transform == "avg_level":
        return trailing_avg_daily(s, full)
    if transform == "avg_log":
        return 100.0 * np.log(trailing_avg_daily(s, full))
    if transform == "eom_neglog":
        return -100.0 * np.log(s.sort_index().reindex(full, method="ffill"))
    raise ValueError(transform)


def main() -> None:
    weights = build_weights_official()
    official = pd.read_csv(OFFICIAL, parse_dates=["date"])
    official["m"] = official["date"].dt.to_period("M")
    official = official.set_index("m")
    build_predicted_bbb(weights, official)                   # ensure data/raw/predicted_bbb.csv

    days = _load_daily("wilshire5000").index                 # daily grid = equity trading days
    proxies = {v[1] for v in VARS}
    start = min(_load_daily(p).index.min() for p in proxies)
    full = pd.date_range(start, days.max(), freq="D")        # daily calendar for A_j lookups

    out = pd.DataFrame(index=days)
    for name, proxy, transform, mcol, _pub in VARS:
        A = trailing_level(proxy, transform, full)
        mult = weights[mcol].to_numpy()
        # A_j at d, d-3mo, ..., d-36mo  (13 points -> 12 quarterly changes)
        Av = [A.reindex(days - pd.DateOffset(months=3 * k)).to_numpy() for k in range(N_LAGS + 1)]
        fac = np.zeros(len(days))
        for k in range(N_LAGS):
            fac -= mult[k] * (Av[k] - Av[k + 1])
        out[name] = fac

    names = [v[0] for v in VARS]
    out["fci_g"] = out[names].sum(axis=1, min_count=len(names))   # defined only where all 7 exist
    out = out.dropna(subset=["fci_g"])

    # ---- test: daily value at each month-end vs the official monthly number ----
    me = out.groupby(out.index.to_period("M")).tail(1)
    cmp = pd.DataFrame({"daily": me["fci_g"].to_numpy(),
                        "official": official["FCI-G Index (baseline)"].reindex(me.index.to_period("M")).to_numpy()},
                       index=me.index).dropna()
    corr = float(np.corrcoef(cmp["daily"], cmp["official"])[0, 1])
    rmse = float(np.sqrt(np.mean((cmp["daily"] - cmp["official"]) ** 2)))

    save = out[["fci_g"] + names].copy()
    save.index = save.index.strftime("%Y-%m-%d")
    save.index.name = "date"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    save.to_csv(OUTPUT_DIR / "fcig_daily_backcast.csv")

    print(f"daily backcast (trailing windows): {len(save)} days, {save.index[0]}..{save.index[-1]}")
    print(f"month-end daily vs official monthly: n={len(cmp)}  corr={corr:.4f}  r2={corr**2:.4f}  rmse={rmse:.4f}")
    print(f"wrote {OUTPUT_DIR / 'fcig_daily_backcast.csv'}")


if __name__ == "__main__":
    main()
