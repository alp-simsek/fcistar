#!/usr/bin/env python3
"""
build_fci_nowcast.py — WEBSITE-ONLY glue: turn the daily FCI-G nowcast into the frontend contract.

Runs after the ported nowcast chain (get_fcig -> daily_backcast -> nowcast_fcig). Reads:
  - this pipeline's daily anchored nowcast  : data/output/fcig_daily_nowcast.csv
  - the official monthly FCI-G (+ 7 factors): data/fcig_official/fci_g_public_monthly_3yr.csv
  - the FCI* model output (quarterly)       : ../data/output/fcistar.csv
Writes:
  ../data/output/fci_nowcast.csv     : date, fci, fci_gap, kind   (kind in {official, nowcast})
  ../data/output/fci_components.csv  : date, kind, <7 component contributions>   (kind in
                                       {quarter, official, nowcast}) — Figure 1 decomposition
and adds nowcast fields to ../data/output/metadata.json.

The FCI line is quarterly through the last estimated quarter (fcistar.csv), then MONTHLY official
points, then DAILY nowcast. The components mirror that cadence and sum to FCI by construction. The
gap holds FCI* at its latest estimate: fci_gap = fci - fcistar_last.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PIPE_OUT = HERE / "data" / "output"
OFFICIAL = HERE / "data" / "fcig_official" / "fci_g_public_monthly_3yr.csv"
SITE_OUT = HERE.parent / "data" / "output"
FCISTAR = SITE_OUT / "fcistar.csv"
META = SITE_OUT / "metadata.json"

OFFICIAL_INDEX_COL = "FCI-G Index (baseline)"
# official component column -> standard name (matches the daily nowcast columns)
OFF_COMP = {"FFR": "ffr", "10Yr Treasury": "t10yr", "Mortgage Rate": "mortgage", "BBB": "bbb",
            "Stock Market": "stock", "House Prices": "house", "Dollar": "dollar"}
COMPS = list(OFF_COMP.values())


def main() -> None:
    fcistar = pd.read_csv(FCISTAR, parse_dates=["date"]).sort_values("date")
    last_quarter = fcistar["date"].max()
    fcistar_last = float(fcistar.loc[fcistar["date"] == last_quarter, "fcistar"].iloc[0])

    off = pd.read_csv(OFFICIAL, parse_dates=["date"]).sort_values("date")
    off_by_month = off.set_index(off["date"].dt.to_period("M"))
    daily = pd.read_csv(PIPE_OUT / "fcig_daily_nowcast.csv", parse_dates=["date"]).sort_values("date")

    # ---- fci_nowcast.csv: recent monthly official + daily nowcast (total + gap) ----
    off_recent = off[off["date"] > last_quarter]
    official_rows = pd.DataFrame({"date": off_recent["date"], "fci": off_recent[OFFICIAL_INDEX_COL],
                                  "fci_gap": off_recent[OFFICIAL_INDEX_COL] - fcistar_last, "kind": "official"})
    last_official = official_rows["date"].max() if not official_rows.empty else last_quarter
    tail = daily[daily["date"] > last_official]
    nowcast_rows = pd.DataFrame({"date": tail["date"], "fci": tail["fci_g"],
                                 "fci_gap": tail["fci_g"] - fcistar_last, "kind": "nowcast"})
    nowcast_through = nowcast_rows["date"].max() if not nowcast_rows.empty else last_official

    out = pd.concat([official_rows, nowcast_rows]).sort_values("date")
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    out.to_csv(SITE_OUT / "fci_nowcast.csv", index=False)

    # ---- fci_components.csv: 7-component decomposition for Figure 1 ----
    #   quarter rows : official components at each FCI* quarter-end (full history, matched by month)
    #   official rows: official components for released months past the quarter
    #   nowcast rows : daily components for the not-yet-released tail
    q = []
    for d in fcistar["date"]:
        p = d.to_period("M")
        if p in off_by_month.index:
            r = off_by_month.loc[p]
            q.append({"date": d, "kind": "quarter", **{c: r[oc] for oc, c in OFF_COMP.items()}})
    qdf = pd.DataFrame(q)
    odf = pd.DataFrame({"date": off_recent["date"], "kind": "official",
                        **{c: off_recent[oc] for oc, c in OFF_COMP.items()}})
    ndf = pd.DataFrame({"date": tail["date"], "kind": "nowcast", **{c: tail[c] for c in COMPS}})
    comp = pd.concat([qdf, odf, ndf]).sort_values("date")
    comp["date"] = comp["date"].dt.strftime("%Y-%m-%d")
    comp[["date", "kind"] + COMPS].to_csv(SITE_OUT / "fci_components.csv", index=False)

    meta = json.loads(META.read_text()) if META.exists() else {}
    meta.update({
        "last_quarter": last_quarter.strftime("%Y-%m-%d"),
        "last_official_date": pd.Timestamp(last_official).strftime("%Y-%m-%d"),
        "nowcast_through": pd.Timestamp(nowcast_through).strftime("%Y-%m-%d"),
        "fcistar_last": round(fcistar_last, 6),
    })
    META.write_text(json.dumps(meta, indent=2))

    print(f"fci_nowcast.csv: {len(official_rows)} official month(s) + {len(nowcast_rows)} nowcast day(s)")
    print(f"fci_components.csv: {len(qdf)} quarter + {len(odf)} official + {len(ndf)} nowcast rows")
    print(f"  last_quarter={meta['last_quarter']}  last_official={meta['last_official_date']}  "
          f"nowcast_through={meta['nowcast_through']}")


if __name__ == "__main__":
    main()
