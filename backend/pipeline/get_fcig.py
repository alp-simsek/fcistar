#!/usr/bin/env python3
"""
get_fcig.py — incremental data fetcher for the daily FCI-G nowcast.

See CLAUDE.md (same folder) for the full plan. This script ONLY fetches/updates data;
it does not estimate weights or build the nowcast.

What it writes (under ./data/):
  fcig_official/fci_g_public_monthly_3yr.csv   official monthly index + 7 sub-indices
  raw/<series>.csv                             one high-frequency proxy per file [date,value]
  raw/metadata.json                            per-series source, last date, last fetch

Incremental update logic: deep history is frozen; on each run we re-fetch only a trailing
REVISION_WINDOW of days and overwrite that tail (FRED, Zillow, and the official FCI-G all get
revised). If a raw CSV does not exist yet, we fetch the full history from HISTORY_START.

Requires: FRED_API_KEY in the environment (same key the website pipeline uses).
  set -a && source ../../../website/.env && set +a    # or export FRED_API_KEY=...

Usage:
  python get_fcig.py                 # update everything
  python get_fcig.py --only dff bbb  # update a subset of raw series
  python get_fcig.py --full          # ignore existing files; refetch full history
"""
from __future__ import annotations

import argparse
import io
import json
import logging
import os
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# --------------------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OFFICIAL_DIR = DATA / "fcig_official"
RAW_DIR = DATA / "raw"

HISTORY_START = "1990-01-01"   # deep history start when a series is fetched for the first time
REVISION_WINDOW = 180          # days of tail to re-fetch on each incremental update

FRED_API_URL = "https://api.stlouisfed.org/fred/series/observations"

# Raw proxy series sourced from FRED: local_name -> FRED series id.
# Per CLAUDE.md Section 4. svenpy10 (GSW, reference) is fetched separately; dgs10 is the daily
# bridge for the live tail. wilshire5000 proxies the Dow Jones U.S. Total Stock Market index.
FRED_SERIES = {
    "dff": "DFF",                    # effective federal funds rate (daily)
    "dgs10": "DGS10",                # 10y CMT par yield (daily) — bridge for svenpy10's tail
    "bbb": "BAMLC0A4CBBBEY",         # ICE BofA BBB effective yield — EXACT Fed series, but FRED
                                     #   serves only a trailing ~3yr (ICE licensing). Used for
                                     #   NOWCASTING the recent BBB factor (anchored to the release).
    "dbaa": "DBAA",                  # Moody's Baa yield (daily, 1986+) — long-history BBB proxy for
                                     #   BACKCASTING & validation. Mapped to the published BBB factor
                                     #   by regression (slope ~1.34); see CLAUDE.md / reproduce_factors.py.
    "dollar_broad": "DTWEXBGS",      # nominal broad USD index (daily, 2006+)
    "mortgage30us": "MORTGAGE30US",  # Freddie 30y fixed mortgage (weekly) — long history (1990+),
                                     #   used for BACKCASTING the mortgage factor.
    "obmmi30": "OBMMIC30YF",         # Optimal Blue OBMMI 30yr conforming (daily, 2017+) — the EXACT
                                     #   Fed source since 2016; used for NOWCASTING (spliced on Freddie).
}

# Official monthly FCI-G file (index + 7 factor contributions). Same host/path the website
# pipeline uses for the quarterly file.
OFFICIAL_FCIG = {
    "filename": "fci_g_public_monthly_3yr.csv",
    "url": (
        "https://www.federalreserve.gov/econres/notes/feds-notes/"
        "fci_g_public_monthly_3yr.csv"
    ),
}

# Gurkaynak-Sack-Wright nominal yield curve (contains SVENPY10). URL VERIFY-AT-RUNTIME:
# the Board has moved this file before. Header has comment rows before the CSV header.
GSW = {
    "local_name": "svenpy10",
    "url": "https://www.federalreserve.gov/data/yield-curve-tables/feds200628.csv",
    "column": "SVENPY10",
    "skiprows_search": "Date",   # find the header row that starts with this token
}

# Equity proxy: the Fed's total-market stock index. ^DWCF (Dow Jones U.S. Total Stock Market) is
# the EXACT index the FCI-G uses and Yahoo keeps it current, so it is the primary/live source.
# ^DWCF only starts in 1995, so we backfill the pre-1995 history with ^W5000 (FT Wilshire 5000,
# the same index back to 1989). The two are numerically the same index (daily log-return corr
# 0.998, level ratio ~1.00 on their 1995-2026 overlap) and reproduce the published Stock factor at
# slope ~1.00, R^2 ~0.996 — far closer than the S&P 500 (slope 0.96). FRED dropped its Wilshire
# series, so Yahoo is the source for both.
#
# History note: ^W5000 was originally the primary (for its 1989 start), but Yahoo froze ^W5000 on
# 2026-06-26, which stalled the whole nowcast because ^W5000 set the daily grid (daily_backcast.py).
# ^DWCF is the maintained series, so it is now the base and ^W5000 supplies only the frozen pre-1995
# tail — a dead ^W5000 no longer affects the live edge. The stored file keeps the name
# wilshire5000.csv (downstream refers to it by that local name). See fetch_yahoo_equity.
YAHOO_EQUITY = {
    "local_name": "wilshire5000",
    "url": "https://query1.finance.yahoo.com/v8/finance/chart/%5EDWCF",           # primary + live tail
    "history_url": "https://query1.finance.yahoo.com/v8/finance/chart/%5EW5000",  # pre-1995 backfill only
}

# Zillow national ZHVI (monthly). URL VERIFY-AT-RUNTIME: Zillow rotates these public CSV paths.
# We filter to RegionName == "United States" and melt wide->long. House prices carry a small
# weight and are held between prints, so this is the least time-critical series.
ZILLOW = {
    "local_name": "zillow_zhvi",
    "url": (
        "https://files.zillowstatic.com/research/public_csvs/zhvi/"
        "Metro_zhvi_uc_sfrcondo_tier_0.33_0.67_sm_sa_month.csv"
    ),
    "region": "United States",
}


def get_logger() -> logging.Logger:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    return logging.getLogger("get_fcig")


# --------------------------------------------------------------------------------------
# Incremental CSV helpers (one file per series, columns: date,value)
# --------------------------------------------------------------------------------------

def _read_existing(path: Path) -> pd.Series | None:
    if not path.exists():
        return None
    df = pd.read_csv(path, parse_dates=["date"])
    return df.set_index("date")["value"].sort_index()


def _fetch_start(existing: pd.Series | None, full: bool) -> str:
    """First date to fetch: full history, or the revision tail before the last stored date."""
    if full or existing is None or existing.empty:
        return HISTORY_START
    last = existing.index.max()
    return (last - timedelta(days=REVISION_WINDOW)).strftime("%Y-%m-%d")


def _merge_save(path: Path, existing: pd.Series | None, fresh: pd.Series, full: bool) -> pd.Series:
    """Overwrite the refetched tail onto frozen history and persist."""
    fresh = fresh.dropna().sort_index()
    if existing is None or full:
        merged = fresh
    else:
        cutoff = fresh.index.min()
        merged = pd.concat([existing[existing.index < cutoff], fresh]).sort_index()
    merged = merged[~merged.index.duplicated(keep="last")]
    out = merged.rename("value").rename_axis("date").reset_index()
    out["date"] = out["date"].dt.strftime("%Y-%m-%d")
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return merged


# --------------------------------------------------------------------------------------
# Source-specific fetchers
# --------------------------------------------------------------------------------------

def fetch_fred(series_id: str, start: str, api_key: str, session: requests.Session) -> pd.Series:
    params = {
        "series_id": series_id,
        "api_key": api_key,
        "file_type": "json",
        "observation_start": start,
        "sort_order": "asc",
    }
    # FRED rate-limits bursts; retry 429/5xx with exponential backoff.
    for attempt in range(5):
        r = session.get(FRED_API_URL, params=params, timeout=60)
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
        break
    else:
        r.raise_for_status()
    obs = r.json().get("observations", [])
    if not obs:
        raise ValueError(f"No observations for {series_id} from {start}")
    df = pd.DataFrame(obs)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"].replace(".", np.nan), errors="coerce")
    return df.dropna(subset=["date"]).set_index("date")["value"].sort_index()


def fetch_gsw_svenpy10(session: requests.Session) -> pd.Series:
    """GSW yield curve CSV -> SVENPY10 daily series. Always full history (small file)."""
    r = session.get(GSW["url"], timeout=120)
    r.raise_for_status()
    text = r.text
    lines = text.splitlines()
    header_idx = next(i for i, ln in enumerate(lines) if ln.startswith(GSW["skiprows_search"]))
    df = pd.read_csv(io.StringIO("\n".join(lines[header_idx:])))
    df = df.rename(columns={df.columns[0]: "date"})
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    s = pd.to_numeric(df[GSW["column"]], errors="coerce")
    return pd.Series(s.values, index=df["date"]).dropna().sort_index()


def _yahoo_chart(session: requests.Session, url: str) -> pd.Series:
    """Daily close series from Yahoo's keyless chart API (full history)."""
    params = {"period1": 0, "period2": 9999999999, "interval": "1d"}
    r = session.get(url, params=params, headers={"User-Agent": "Mozilla/5.0"}, timeout=60)
    r.raise_for_status()
    res = r.json()["chart"]["result"][0]
    ts = pd.to_datetime(res["timestamp"], unit="s").normalize()
    close = res["indicators"]["quote"][0]["close"]
    return pd.Series(close, index=ts).dropna().sort_index()


def fetch_yahoo_equity(session: requests.Session) -> pd.Series:
    """Fed total-market equity index: ^DWCF, backfilled before its 1995 start with ^W5000.

    ^DWCF (Dow Jones U.S. Total Stock Market — the exact index the FCI-G uses) is the live source
    of truth. It starts in 1995, so for dates before ^DWCF's first print we splice in ^W5000 (the
    same index back to 1989), anchored to ^DWCF at ^DWCF's start by ratio:
        wilshire(d) = dwcf(start) * w5(d) / w5(start)   for d < start.
    The FCI uses only 3-month log-CHANGES, so the ~1% level offset cancels and the join adds no
    jump. ^W5000 is frozen at Yahoo (last print 2026-06-26) but that no longer matters — it feeds
    only the unchanging pre-1995 history. If ^W5000 is unavailable, ^DWCF alone (1995+) is used.
    """
    dw = _yahoo_chart(session, YAHOO_EQUITY["url"])
    history_url = YAHOO_EQUITY.get("history_url")
    if history_url:
        try:
            w5 = _yahoo_chart(session, history_url)
            start = dw.index.min()
            anchor = w5.asof(start)
            pre = w5.index[w5.index < start]
            if len(pre) and pd.notna(anchor) and anchor > 0:
                hist = float(dw.loc[start]) * w5.loc[pre] / float(anchor)
                dw = pd.concat([hist, dw]).sort_index()
        except Exception as exc:  # noqa: BLE001 — ^DWCF (1995+) alone if the backfill fails
            logging.getLogger("get_fcig").warning(
                "^W5000 history backfill failed (%s); using ^DWCF from %s only",
                exc, dw.index.min().date())
    return dw


def fetch_zillow_zhvi(session: requests.Session) -> pd.Series:
    """Zillow ZHVI national monthly series (wide -> long)."""
    r = session.get(ZILLOW["url"], timeout=120)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text))
    row = df[df["RegionName"] == ZILLOW["region"]]
    if row.empty:
        raise ValueError(f"Region {ZILLOW['region']!r} not found in Zillow file")
    date_cols = [c for c in df.columns if c[:4].isdigit() and "-" in c]
    s = row[date_cols].iloc[0]
    s.index = pd.to_datetime(s.index, errors="coerce")
    return pd.to_numeric(s, errors="coerce").dropna().sort_index()


def download_official(logger: logging.Logger) -> Path:
    out = OFFICIAL_DIR / OFFICIAL_FCIG["filename"]
    out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading official monthly FCI-G -> %s", out.name)
    r = requests.get(OFFICIAL_FCIG["url"], timeout=60)
    r.raise_for_status()
    out.write_bytes(r.content)
    return out


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------

def update_raw(only: list[str] | None, full: bool, logger: logging.Logger) -> dict:
    session = requests.Session()
    meta: dict[str, dict] = {}
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    def record(name: str, series: pd.Series, source: str) -> None:
        meta[name] = {
            "source": source,
            "last_date": series.index.max().strftime("%Y-%m-%d"),
            "n_obs": int(series.shape[0]),
            "last_fetch": now,
        }

    # FRED series (incremental)
    api_key = os.getenv("FRED_API_KEY")
    for name, sid in FRED_SERIES.items():
        if only and name not in only:
            continue
        if not api_key:
            raise EnvironmentError("FRED_API_KEY is not set (source website/.env).")
        path = RAW_DIR / f"{name}.csv"
        existing = None if full else _read_existing(path)
        start = _fetch_start(existing, full)
        logger.info("FRED %s (%s) from %s", name, sid, start)
        fresh = fetch_fred(sid, start, api_key, session)
        merged = _merge_save(path, existing, fresh, full)
        record(name, merged, f"FRED:{sid}")
        time.sleep(0.5)  # be gentle on the FRED rate limit

    # GSW SVENPY10, Yahoo equity, Zillow ZHVI (full refetch; small files, defensive on URL drift)
    for cfg, fetcher in ((GSW, fetch_gsw_svenpy10),
                         (YAHOO_EQUITY, fetch_yahoo_equity),
                         (ZILLOW, fetch_zillow_zhvi)):
        name = cfg["local_name"]
        if only and name not in only:
            continue
        path = RAW_DIR / f"{name}.csv"
        try:
            logger.info("Fetching %s from %s", name, cfg["url"])
            fresh = fetcher(session)
            merged = _merge_save(path, None, fresh, True)
            record(name, merged, cfg["url"])
        except Exception as exc:  # noqa: BLE001 — keep the rest of the run alive
            logger.warning("FAILED %s (%s): %s — VERIFY URL in get_fcig.py", name, cfg["url"], exc)

    (RAW_DIR / "metadata.json").write_text(json.dumps(meta, indent=2))
    return meta


def main() -> None:
    ap = argparse.ArgumentParser(description="Incremental data fetch for the daily FCI-G nowcast.")
    ap.add_argument("--only", nargs="*", default=None, help="subset of raw series local-names")
    ap.add_argument("--full", action="store_true", help="ignore existing files; refetch full history")
    ap.add_argument("--skip-official", action="store_true", help="do not redownload the official FCI-G file")
    args = ap.parse_args()

    logger = get_logger()
    if not args.skip_official and not args.only:
        download_official(logger)
    meta = update_raw(args.only, args.full, logger)
    logger.info("Done. Updated %d raw series.", len(meta))
    for name, m in sorted(meta.items()):
        logger.info("  %-14s last=%s  n=%d", name, m["last_date"], m["n_obs"])


if __name__ == "__main__":
    main()
