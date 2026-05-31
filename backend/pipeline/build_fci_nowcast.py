#!/usr/bin/env python3
"""
build_fci_nowcast.py — WEBSITE-ONLY glue: turn the daily FCI-G nowcast into the frontend contract.

Runs after the ported nowcast chain (get_fcig -> daily_backcast -> nowcast_fcig). Reads:
  - this pipeline's daily anchored nowcast  : data/output/fcig_daily_nowcast.csv
  - the official monthly FCI-G              : data/fcig_official/fci_g_public_monthly_3yr.csv
  - the FCI* model output (quarterly)       : ../data/output/fcistar.csv
Writes the recent-extension series the frontend appends to fcistar.csv:
  ../data/output/fci_nowcast.csv  : date, fci, fci_gap, kind   (kind in {official, nowcast})
and adds nowcast fields to ../data/output/metadata.json.

Plotting intent (see website/CLAUDE.md): the FCI line is quarterly through the last estimated
quarter (from fcistar.csv, where FCI* exists), then MONTHLY official points for released months past
the quarter (kind=official), then DAILY nowcast for the not-yet-released tail (kind=nowcast). No
daily points between the last quarter and the monthly official (avoids crowding). The gap holds FCI*
at its latest estimate: fci_gap = fci - fcistar_last.
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


def main() -> None:
    fcistar = pd.read_csv(FCISTAR, parse_dates=["date"]).sort_values("date")
    last_quarter = fcistar["date"].max()
    fcistar_last = float(fcistar.loc[fcistar["date"] == last_quarter, "fcistar"].iloc[0])

    # official MONTHLY points released after the last estimated quarter
    off = pd.read_csv(OFFICIAL, parse_dates=["date"]).sort_values("date")
    off = off[off["date"] > last_quarter]
    official_rows = pd.DataFrame({
        "date": off["date"],
        "fci": off[OFFICIAL_INDEX_COL],
        "fci_gap": off[OFFICIAL_INDEX_COL] - fcistar_last,
        "kind": "official",
    })
    last_official = official_rows["date"].max() if not official_rows.empty else last_quarter

    # DAILY anchored nowcast for the not-yet-released tail (after the last official month-end)
    daily = pd.read_csv(PIPE_OUT / "fcig_daily_nowcast.csv", parse_dates=["date"]).sort_values("date")
    tail = daily[daily["date"] > last_official]
    nowcast_rows = pd.DataFrame({
        "date": tail["date"],
        "fci": tail["fci_g"],
        "fci_gap": tail["fci_g"] - fcistar_last,
        "kind": "nowcast",
    })
    nowcast_through = nowcast_rows["date"].max() if not nowcast_rows.empty else last_official

    out = pd.concat([official_rows, nowcast_rows]).sort_values("date")
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    out.to_csv(SITE_OUT / "fci_nowcast.csv", index=False)

    meta = json.loads(META.read_text()) if META.exists() else {}
    meta.update({
        "last_quarter": last_quarter.strftime("%Y-%m-%d"),
        "last_official_date": pd.Timestamp(last_official).strftime("%Y-%m-%d"),
        "nowcast_through": pd.Timestamp(nowcast_through).strftime("%Y-%m-%d"),
        "fcistar_last": round(fcistar_last, 6),
    })
    META.write_text(json.dumps(meta, indent=2))

    print(f"fci_nowcast.csv: {len(official_rows)} official month(s) + {len(nowcast_rows)} nowcast day(s)")
    print(f"  last_quarter={meta['last_quarter']}  last_official={meta['last_official_date']}  "
          f"nowcast_through={meta['nowcast_through']}  fcistar_last={meta['fcistar_last']}")


if __name__ == "__main__":
    main()
