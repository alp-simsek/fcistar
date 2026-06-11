#!/usr/bin/env python3
"""
*******************************************************************************
**** Download the Survey of Professional Forecasters (SPF) median forecasts
**** for real GDP growth and core PCE inflation, and reshape them into a
**** filter-ready panel.
****
**** PORTED FROM THE RESEARCH COPY at code/empirics/SPF/get_spf.py (FCIstar
**** Dropbox). That copy is canonical: make changes there first, then port
**** here. The two differ only in I/O paths.
****
**** Source: Federal Reserve Bank of Philadelphia, SPF data files (public).
****   - Median_RGDP_Growth.xlsx  (sheet 'Median_Growth'):  DRGDP2..DRGDP6,
****       annualized % growth of real GDP.
****   - Median_COREPCE_Level.xlsx (sheet 'Median_Level'):  COREPCE1..COREPCE6,
****       annualized % core PCE inflation rate.
****
**** Each survey row is (YEAR, QUARTER) = the quarter the survey was taken in.
**** The horizon columns are forecasts for that quarter and the next four:
****       suffix 2 -> survey quarter itself (horizon 0, a "nowcast"),
****       suffix 3 -> +1 quarter, ... suffix 6 -> +4 quarters.
****
**** Output (under backend/pipeline/data/spf/, gitignored):
****   raw/median_rgdp_growth.csv     parsed wide file (one row per survey)
****   raw/median_corepce_level.csv   parsed wide file (one row per survey)
****   output/spf_forecasts.csv       tidy long panel, one row per
****       (survey_quarter, target_quarter) with median drgdp and corepce.
****
**** Selecting a forecast: for a target quarter Q, use the row with the LATEST
**** survey_quarter S <= Q (horizon Q - S in 0..4). Always available -- at worst
**** the survey from one quarter ago (horizon 1).
*******************************************************************************
"""
from __future__ import annotations

import io
import logging
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

# Phil Fed SPF median data files. The `hash` query parameter is part of the
# published download URL; if the Phil Fed re-publishes a file the hash changes,
# so a 404 here is the signal to refresh these URLs.
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

MAX_HORIZON = 4  # survey covers the current quarter (h=0) through +4


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
    """Download one SPF file and return its parsed wide sheet (one row per survey)."""
    spec = SERIES[key]
    logger.info("Downloading %s", spec["name"])
    raw = _download(spec["url"])
    clean = _sanitize_xlsx(raw)
    df = pd.read_excel(clean, sheet_name=spec["sheet"])
    df = df.dropna(subset=["YEAR", "QUARTER"]).copy()
    df["YEAR"] = df["YEAR"].astype(int)
    df["QUARTER"] = df["QUARTER"].astype(int)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(RAW_DIR / spec["raw_csv"], index=False)
    return df


def melt_to_long(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """
    Melt a wide survey sheet into one row per (survey_quarter, target_quarter).

    For horizon h in 0..MAX_HORIZON the forecast lives in column f"{prefix}{2+h}"
    and refers to the quarter h steps after the survey quarter.
    """
    spec = SERIES[key]
    prefix = spec["prefix"]
    records = []
    for _, row in df.iterrows():
        sy, sq = int(row["YEAR"]), int(row["QUARTER"])
        for h in range(0, MAX_HORIZON + 1):
            col = f"{prefix}{2 + h}"
            if col not in df.columns:
                continue
            val = row[col]
            if pd.isna(val):
                continue
            t_index = (sy * 4 + (sq - 1)) + h
            ty, tq = divmod(t_index, 4)
            records.append(
                {
                    "survey_year": sy,
                    "survey_q": sq,
                    "target_year": ty,
                    "target_q": tq + 1,
                    "horizon": h,
                    key: float(val),
                }
            )
    return pd.DataFrame.from_records(records)


def _qlabel(year: int, q: int) -> str:
    return f"{year}Q{q}"


def build_panel() -> pd.DataFrame:
    """Download both series and merge into the tidy long forecast panel."""
    drgdp_long = melt_to_long(fetch_wide("drgdp"), "drgdp")
    corepce_long = melt_to_long(fetch_wide("corepce"), "corepce")

    keys = ["survey_year", "survey_q", "target_year", "target_q", "horizon"]
    panel = pd.merge(drgdp_long, corepce_long, on=keys, how="outer")

    panel["survey_quarter"] = [
        _qlabel(y, q) for y, q in zip(panel["survey_year"], panel["survey_q"])
    ]
    panel["target_quarter"] = [
        _qlabel(y, q) for y, q in zip(panel["target_year"], panel["target_q"])
    ]
    panel = panel.sort_values(
        ["target_year", "target_q", "survey_year", "survey_q"]
    ).reset_index(drop=True)

    cols = [
        "target_quarter", "survey_quarter", "horizon",
        "drgdp", "corepce",
        "target_year", "target_q", "survey_year", "survey_q",
    ]
    return panel[cols]


def latest_forecast_for(panel: pd.DataFrame, target_year: int, target_q: int) -> dict | None:
    """
    Return the most recent survey's forecast for the given target quarter
    (largest survey quarter with horizon in 0..MAX_HORIZON), or None.
    """
    sub = panel[(panel["target_year"] == target_year) & (panel["target_q"] == target_q)]
    sub = sub[(sub["horizon"] >= 0) & (sub["horizon"] <= MAX_HORIZON)]
    if sub.empty:
        return None
    row = sub.sort_values("horizon").iloc[0]
    return {
        "target_quarter": row["target_quarter"],
        "survey_quarter": row["survey_quarter"],
        "horizon": int(row["horizon"]),
        "drgdp": None if pd.isna(row["drgdp"]) else float(row["drgdp"]),
        "corepce": None if pd.isna(row["corepce"]) else float(row["corepce"]),
    }


def main() -> None:
    panel = build_panel()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "spf_forecasts.csv"
    panel.to_csv(out_path, index=False)
    logger.info("Wrote %d rows -> %s", len(panel), out_path)

    last_survey = panel.sort_values(["survey_year", "survey_q"]).iloc[-1]
    sy, sq = int(last_survey["survey_year"]), int(last_survey["survey_q"])
    logger.info("Latest survey in file: %s", _qlabel(sy, sq))


if __name__ == "__main__":
    main()
