#!/usr/bin/env python3
"""
nowcast_fcig.py — ANCHORED daily FCI-G nowcast.

For each sub-index j and date d:
    anchored_j(d) = official_j(m) + ( recon_j(d) - recon_j(end of m) )
where m = the latest official month whose month-end <= d, and recon is the unanchored daily
reconstruction (data/output/fcig_daily_backcast.csv). The total is the SUM of the anchored
sub-indices (equivalent to anchoring the total, since both official and recon sum across factors).

This passes exactly through every official monthly point (so per-factor proxy gaps like mortgage
disappear from the level) and uses the reconstruction only for the within-period / forward
DIFFERENCES. The stretch after the last official release is the live nowcast.

Output: data/output/fcig_daily_nowcast.csv  (date, fci_g, + 7 factor contributions).
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "data" / "output"
RECON = OUT / "fcig_daily_backcast.csv"
OFFICIAL = HERE / "data" / "fcig_official" / "fci_g_public_monthly_3yr.csv"

OFF_MAP = {"FFR": "ffr", "10Yr Treasury": "t10yr", "Mortgage Rate": "mortgage", "BBB": "bbb",
           "Stock Market": "stock", "House Prices": "house", "Dollar": "dollar"}
FACS = ["ffr", "t10yr", "mortgage", "bbb", "stock", "house", "dollar"]


def main() -> None:
    recon = pd.read_csv(RECON, parse_dates=["date"], index_col="date")
    off = pd.read_csv(OFFICIAL, parse_dates=["date"]).rename(columns=OFF_MAP)
    off = off[["date"] + FACS].sort_values("date")
    off["period"] = off["date"].dt.to_period("M")          # official date is the month-end

    # recon at the last trading day of each month (the FCI^now_m reference), keyed by month
    recon_end = recon.groupby(recon.index.to_period("M")).tail(1)
    recon_end = recon_end.set_index(recon_end.index.to_period("M"))[FACS]

    # for each daily date, the latest official month-end <= date  ->  (official values, anchor month)
    left = pd.DataFrame({"date": recon.index}).sort_values("date")
    asof = pd.merge_asof(left, off.rename(columns={"date": "mend"}),
                         left_on="date", right_on="mend", direction="backward").set_index("date")
    ref = recon_end.reindex(asof["period"]).set_axis(asof.index)   # recon at end of the anchor month

    anc = pd.DataFrame(index=recon.index)
    for c in FACS:
        anc[c] = asof[c].to_numpy() + recon[c].to_numpy() - ref[c].to_numpy()
    anc["fci_g"] = anc[FACS].sum(axis=1)                    # total = sum of anchored sub-indices
    anc = anc.dropna(subset=["fci_g"])

    out = anc[["fci_g"] + FACS].copy()
    out.index = out.index.strftime("%Y-%m-%d")
    out.index.name = "date"
    out.to_csv(OUT / "fcig_daily_nowcast.csv")

    # ---- sanity: anchored hits official exactly at trading-day month-ends ----
    om = off.set_index("period")
    chk = []
    for p, g in recon.groupby(recon.index.to_period("M")):
        d_end = g.index.max()
        if p in om.index and pd.Timestamp(d_end) == om.loc[p, "date"]:   # month-end is a trading day
            chk.append(abs(anc["fci_g"].get(d_end, np.nan) - om.loc[p, FACS].sum()))
    last_off = off["date"].max()
    tail = anc.index[anc.index > last_off]
    print(f"anchored nowcast: {len(out)} days, {out.index[0]}..{out.index[-1]}")
    print(f"passes through official at month-ends: max |anchored - official| = {np.nanmax(chk):.2e} (n={len(chk)})")
    print(f"live nowcast tail (after last official {last_off.date()}): {len(tail)} days; "
          f"latest fci_g = {anc['fci_g'].iloc[-1]:.3f}")
    print(f"wrote {OUT / 'fcig_daily_nowcast.csv'}")


if __name__ == "__main__":
    main()
