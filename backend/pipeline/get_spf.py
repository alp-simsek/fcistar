#!/usr/bin/env python3
"""
*******************************************************************************
**** Download the Survey of Professional Forecasters (SPF) forecasts for real
**** GDP growth and core PCE inflation, and reshape them into a filter-ready
**** panel. Pulls both the MEDIAN forecast and the cross-sectional 25th/75th
**** PERCENTILES (for the FCI* sensitivity analysis).
****
**** PORTED FROM THE RESEARCH COPY at code/empirics/SPF/get_spf.py (FCIstar
**** Dropbox). That copy is canonical: change it there first, then port here.
**** The two differ only in I/O paths.
****
**** Source: Federal Reserve Bank of Philadelphia, SPF data files (public).
****   Median forecasts:
****     - Median_RGDP_Growth.xlsx  (sheet 'Median_Growth'):  DRGDP2..DRGDP6,
****         annualized % growth of real GDP.
****     - Median_COREPCE_Level.xlsx (sheet 'Median_Level'):  COREPCE1..COREPCE6,
****         annualized % core PCE inflation rate.
****   Cross-sectional dispersion (25th / 75th percentiles across forecasters):
****     - Dispersion_RGDP.xlsx (sheet 'D2'):  P25/P75 of Q/Q growth (ann. %).
****     - Dispersion_COREPCE.xlsx (sheet 'D1'): P25/P75 of the level (ann. %).
****
**** Each survey row is (YEAR, QUARTER) = the quarter the survey was taken in.
**** Horizon columns are forecasts for that quarter and the next four (h=0..4).
**** Median files: column suffix 2+h (suffix 2 = horizon 0 = survey quarter).
**** Dispersion files: per horizon a (P25, P75, IQR) triple after the date, so
**** P25 = col 1+3h, P75 = col 2+3h (h=0 = survey quarter).
****
**** Output:
****   data/raw/*.csv                  parsed wide files (one row per survey)
****   data/output/spf_forecasts.csv   tidy long panel, one row per
****       (survey_quarter, target_quarter) with drgdp / corepce at the 25th,
****       50th (median), 75th percentiles.
****
**** Units are SPF native annualized percent. The GDP growth becomes a quarterly
**** log change (drgdp/400) only later, in the forecast step.
****
**** Selection: for a target quarter Q, use the row with the LATEST survey_quarter
**** S <= Q (horizon Q - S in 0..4) -- always available (at worst the survey from
**** one quarter ago, horizon 1).
*******************************************************************************
"""
from __future__ import annotations

import io
import logging
import re
import time
import zipfile
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("get_spf")

HERE = Path(__file__).resolve().parent
RAW_DIR = HERE / "data" / "spf" / "raw"
OUT_DIR = HERE / "data" / "spf" / "output"

# Phil Fed SPF data files. The `hash` query parameter is part of the published
# download URL; if the Phil Fed re-publishes a file the hash changes, so a 404
# here is the signal to refresh these URLs.
SERIES = {
    "drgdp": {
        "name": "Real GDP growth (annualized %, median)",
        "url": (
            "https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/"
            "survey-of-professional-forecasters/data-files/files/Median_RGDP_Growth.xlsx"
            "?sc_lang=en&hash=084D7EAE781F8DBAB17975CDB199BFBB"
        ),
        "sheet": "Median_Growth",
        "prefix": "DRGDP",
        "raw_csv": "median_rgdp_growth.csv",
    },
    "corepce": {
        "name": "Core PCE inflation (annualized %, median)",
        "url": (
            "https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/"
            "survey-of-professional-forecasters/data-files/files/Median_COREPCE_Level.xlsx"
            "?sc_lang=en&hash=FA36DE2D0348728E2CA5C172FAFB1545"
        ),
        "sheet": "Median_Level",
        "prefix": "COREPCE",
        "raw_csv": "median_corepce_level.csv",
    },
}

# Cross-sectional dispersion files (25th / 75th percentiles across forecasters).
# Same horizon convention; per horizon the sheet has a (P25, P75, IQR) triple.
DISPERSION = {
    "drgdp": {
        "name": "Real GDP growth dispersion (P25/P75 of Q/Q growth)",
        "url": (
            "https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/"
            "survey-of-professional-forecasters/data-files/files/Dispersion_RGDP.xlsx"
            "?sc_lang=en&hash=346EE80BCC5FAC8F41F39289A86189DD"
        ),
        "sheet": "D2",  # 75th/25th percentiles of Q/Q growth, annualized %
        "raw_csv": "dispersion_rgdp_growth.csv",
    },
    "corepce": {
        "name": "Core PCE dispersion (P25/P75 of the level)",
        "url": (
            "https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/"
            "survey-of-professional-forecasters/data-files/files/Dispersion_COREPCE.xlsx"
            "?sc_lang=en&hash=D398DE7B9C9586744A0F91A5639EFA81"
        ),
        "sheet": "D1",  # 75th/25th percentiles of the annualized-rate level
        "raw_csv": "dispersion_corepce.csv",
    },
}

MAX_HORIZON = 4  # survey covers the current quarter (h=0) through +4
_QPAT = re.compile(r"^\s*(\d{4})Q([1-4])\s*$")


def _download(url: str) -> bytes:
    """GET with retry/backoff (the Phil Fed CDN occasionally rate-limits/5xx)."""
    for attempt in range(5):
        r = requests.get(url, timeout=120)
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
        return r.content
    r.raise_for_status()
    return r.content


def _sanitize_xlsx(raw: bytes) -> io.BytesIO:
    """
    Some SPF .xlsx carry a `modified` date in docProps/core.xml that openpyxl
    rejects (it expects a datetime, the file stores a date), which aborts the
    whole read. Rewrite that one metadata part with minimal valid properties.
    Everything else in the workbook is copied through untouched.
    """
    MIN_CORE = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<cp:coreProperties'
        ' xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties"'
        ' xmlns:dc="http://purl.org/dc/elements/1.1/"></cp:coreProperties>'
    )
    out = io.BytesIO()
    with zipfile.ZipFile(io.BytesIO(raw)) as zin, zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.namelist():
            data = zin.read(item)
            if item == "docProps/core.xml":
                data = MIN_CORE.encode("utf-8")
            zout.writestr(item, data)
    out.seek(0)
    return out


def fetch_wide(key: str) -> pd.DataFrame:
    """Download one MEDIAN file and return its parsed wide sheet (one row per survey)."""
    spec = SERIES[key]
    logger.info("Downloading %s", spec["name"])
    clean = _sanitize_xlsx(_download(spec["url"]))
    df = pd.read_excel(clean, sheet_name=spec["sheet"])
    df = df.dropna(subset=["YEAR", "QUARTER"]).copy()
    df["YEAR"] = df["YEAR"].astype(int)
    df["QUARTER"] = df["QUARTER"].astype(int)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(RAW_DIR / spec["raw_csv"], index=False)
    return df


def fetch_dispersion(key: str) -> pd.DataFrame:
    """Download one DISPERSION file and return the raw sheet (header=None; title rows kept)."""
    spec = DISPERSION[key]
    logger.info("Downloading %s", spec["name"])
    clean = _sanitize_xlsx(_download(spec["url"]))
    raw = pd.read_excel(clean, sheet_name=spec["sheet"], header=None)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    raw.to_csv(RAW_DIR / spec["raw_csv"], index=False, header=False)
    return raw


def _emit(records, sy, sq, h, payload):
    """Append one (survey, target, horizon) record with the given value payload."""
    t_index = (sy * 4 + (sq - 1)) + h
    ty, tq = divmod(t_index, 4)
    records.append({"survey_year": sy, "survey_q": sq, "target_year": ty,
                    "target_q": tq + 1, "horizon": h, **payload})


def melt_to_long(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """Melt a wide MEDIAN sheet to one row per (survey, target). Median lives in column 2+h."""
    prefix = SERIES[key]["prefix"]
    records = []
    for _, row in df.iterrows():
        sy, sq = int(row["YEAR"]), int(row["QUARTER"])
        for h in range(MAX_HORIZON + 1):
            col = f"{prefix}{2 + h}"
            if col in df.columns and not pd.isna(row[col]):
                _emit(records, sy, sq, h, {key: float(row[col])})
    return pd.DataFrame.from_records(records)


def melt_dispersion(raw: pd.DataFrame, key: str) -> pd.DataFrame:
    """
    Melt a DISPERSION sheet to one row per (survey, target) with P25/P75.
    Rows are found by a 'YYYYQq' survey date in column 0; per horizon h the
    columns are P25 = 1+3h, P75 = 2+3h (the third of each triple is the IQR).
    """
    records = []
    for _, row in raw.iterrows():
        sd = row[0]
        if not isinstance(sd, str):
            continue
        m = _QPAT.match(sd)
        if not m:
            continue
        sy, sq = int(m.group(1)), int(m.group(2))
        for h in range(MAX_HORIZON + 1):
            p25, p75 = row[1 + 3 * h], row[2 + 3 * h]
            if pd.isna(p25) and pd.isna(p75):
                continue
            _emit(records, sy, sq, h, {
                f"{key}_p25": None if pd.isna(p25) else float(p25),
                f"{key}_p75": None if pd.isna(p75) else float(p75),
            })
    return pd.DataFrame.from_records(records)


def _qlabel(year: int, q: int) -> str:
    return f"{year}Q{q}"


def build_panel() -> pd.DataFrame:
    """Download median + dispersion for both series and merge into one tidy panel."""
    keys = ["survey_year", "survey_q", "target_year", "target_q", "horizon"]
    panel = None
    for k in ("drgdp", "corepce"):
        med = melt_to_long(fetch_wide(k), k)                 # k = median (p50)
        disp = melt_dispersion(fetch_dispersion(k), k)        # k_p25, k_p75
        merged = pd.merge(med, disp, on=keys, how="outer")
        panel = merged if panel is None else pd.merge(panel, merged, on=keys, how="outer")

    panel["survey_quarter"] = [_qlabel(y, q) for y, q in zip(panel["survey_year"], panel["survey_q"])]
    panel["target_quarter"] = [_qlabel(y, q) for y, q in zip(panel["target_year"], panel["target_q"])]
    panel = panel.sort_values(["target_year", "target_q", "survey_year", "survey_q"]).reset_index(drop=True)

    cols = [
        "target_quarter", "survey_quarter", "horizon",
        "drgdp", "drgdp_p25", "drgdp_p75",
        "corepce", "corepce_p25", "corepce_p75",
        "target_year", "target_q", "survey_year", "survey_q",
    ]
    for c in cols:
        if c not in panel.columns:
            panel[c] = None
    return panel[cols]


def latest_forecast_for(panel: pd.DataFrame, target_year: int, target_q: int) -> dict | None:
    """
    Most recent survey's forecast for the target quarter (largest survey quarter,
    horizon 0..MAX_HORIZON). Returns medians and 25th/75th percentiles, or None.
    """
    sub = panel[(panel["target_year"] == target_year) & (panel["target_q"] == target_q)]
    sub = sub[(sub["horizon"] >= 0) & (sub["horizon"] <= MAX_HORIZON)]
    if sub.empty:
        return None
    row = sub.sort_values("horizon").iloc[0]   # smallest horizon = most recent survey

    def g(col):
        v = row.get(col)
        return None if v is None or pd.isna(v) else float(v)

    return {
        "target_quarter": row["target_quarter"],
        "survey_quarter": row["survey_quarter"],
        "horizon": int(row["horizon"]),
        "drgdp": g("drgdp"), "drgdp_p25": g("drgdp_p25"), "drgdp_p75": g("drgdp_p75"),
        "corepce": g("corepce"), "corepce_p25": g("corepce_p25"), "corepce_p75": g("corepce_p75"),
    }


def main() -> None:
    panel = build_panel()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "spf_forecasts.csv"
    panel.to_csv(out_path, index=False)
    logger.info("Wrote %d rows -> %s", len(panel), out_path)

    last = panel.sort_values(["survey_year", "survey_q"]).iloc[-1]
    sy, sq = int(last["survey_year"]), int(last["survey_q"])
    logger.info("Latest survey in file: %s", _qlabel(sy, sq))
    fc = latest_forecast_for(panel, sy, sq)
    if fc:
        logger.info("  %s: drgdp[25/50/75]=%.2f/%.2f/%.2f  corepce[25/50/75]=%.2f/%.2f/%.2f",
                    fc["target_quarter"], fc["drgdp_p25"], fc["drgdp"], fc["drgdp_p75"],
                    fc["corepce_p25"], fc["corepce"], fc["corepce_p75"])


if __name__ == "__main__":
    main()
