from __future__ import annotations

"""
Assemble the quarterly dataset used by the FCI* estimation pipeline.

This script plays the role of the MATLAB assemble_data_covid.m code:
- download / read the raw source series
- convert them to a common quarterly panel
- construct data_hlm and covid_dummies
- save the intermediate files used by estimate_CCS_struct
"""

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import requests

FRED_API_URL = "https://api.stlouisfed.org/fred/series/observations"

# FRED series used in the assembled estimation panel.
FRED_SERIES = {
    "GDPC1": "real_gdp",
    "PCEPILFE": "core_pce_index",
    "FEDFUNDS": "fed_funds_rate",
    "INTDSRUSM193N": "ny_discount_rate",
    "NFCI": "chicago_fci",
    "GDPPOT": "cbo_potential_gdp",
}

# Public Fed FCI files used in place of the static MATLAB inputs.
FED_FCI_FILES = {
    "fci_g_public_quarterly_3yr.csv": {
        "url": (
            "https://www.federalreserve.gov/econres/notes/feds-notes/"
            "fci_g_public_quarterly_3yr.csv"
        ),
        "expected_value_col": "FCI-G Index (baseline)",
    },
    "fci_g_public_quarterly_1yr.csv": {
        "url": (
            "https://www.federalreserve.gov/econres/notes/feds-notes/"
            "fci_g_public_quarterly_1yr.csv"
        ),
        "expected_value_col": "FCI-G Index (one-year lookback)",
    },
}


def get_logger() -> logging.Logger:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    return logging.getLogger("assemble_data")


def repo_root() -> Path:
    """Assumes this file lives at <repo>/backend/pipeline/assemble_data.py."""
    return Path(__file__).resolve().parents[2]


def backend_root() -> Path:
    return repo_root() / "backend"


def raw_root() -> Path:
    return backend_root() / "data" / "raw"


def quarter_index(start_date: str, end_date: str) -> pd.PeriodIndex:
    return pd.period_range(
        start=pd.Timestamp(start_date).to_period("Q"),
        end=pd.Timestamp(end_date).to_period("Q"),
        freq="Q",
    )


def quarter_end_dates(periods: pd.PeriodIndex) -> pd.DatetimeIndex:
    return periods.to_timestamp(how="end").normalize()


def quarter_decimal(periods: pd.PeriodIndex) -> np.ndarray:
    return np.array([p.year + 0.25 * (p.quarter - 1) for p in periods], dtype=float)


def fetch_fred_series(
    series_id: str,
    start_date: str,
    end_date: str,
    api_key: str,
    session: requests.Session,
) -> pd.Series:
    """Download one FRED series at its native frequency."""
    params = {
        "series_id": series_id,
        "api_key": api_key,
        "file_type": "json",
        "observation_start": start_date,
        "observation_end": end_date,
        "sort_order": "asc",
    }
    r = session.get(FRED_API_URL, params=params, timeout=60)
    r.raise_for_status()
    payload = r.json()

    obs = payload.get("observations", [])
    if not obs:
        raise ValueError(f"No observations returned for {series_id}")

    df = pd.DataFrame(obs)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"].replace(".", np.nan), errors="coerce")
    df = df.dropna(subset=["date"])

    s = df.set_index("date")["value"].sort_index()
    s.name = series_id
    return s


def to_quarterly_mean(series: pd.Series) -> pd.Series:
    """Convert a dated series to quarterly means."""
    out = series.groupby(series.index.to_period("Q")).mean()
    out.index = pd.PeriodIndex(out.index, freq="Q")
    out.name = series.name
    return out


def build_fred_panel(
    series_ids: Iterable[str],
    start_date: str,
    end_date: str,
    api_key: str,
    periods: pd.PeriodIndex,
    logger: logging.Logger,
) -> pd.DataFrame:
    """Download the requested FRED series and align them to the common quarterly index."""
    session = requests.Session()
    panel = pd.DataFrame(index=periods)

    for sid in series_ids:
        logger.info("Downloading %s", sid)
        raw = fetch_fred_series(sid, start_date, end_date, api_key, session)
        q = to_quarterly_mean(raw)
        panel[sid] = q.reindex(periods)

    return panel


def download_file(url: str, out_path: Path, logger: logging.Logger) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s -> %s", url, out_path)
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    out_path.write_bytes(r.content)


def download_new_fci_files(fci_dir: Path, logger: logging.Logger) -> dict[str, Path]:
    """Download the public quarterly FCI-G files from the Fed website."""
    downloaded = {}
    for filename, meta in FED_FCI_FILES.items():
        out_path = fci_dir / filename
        download_file(meta["url"], out_path, logger)
        downloaded[filename] = out_path
    return downloaded


def infer_date_col(df: pd.DataFrame) -> str:
    preferred = {"date", "time", "datetime", "observation_date"}
    for col in df.columns:
        if str(col).strip().lower() in preferred:
            return col

    raise ValueError(
        "Could not infer date column. Expected one of: "
        "date, time, datetime, observation_date."
    )


def validate_fci_csv_schema(
    path: Path,
    expected_value_col: str,
) -> tuple[pd.DataFrame, str]:
    """
    Validate that the downloaded FCI CSV has the expected date column and
    value column.
    """
    df = pd.read_csv(path)

    if df.empty:
        raise ValueError(f"{path} is empty.")

    date_col = infer_date_col(df)

    if expected_value_col not in df.columns:
        raise ValueError(
            f"{path} is missing expected column '{expected_value_col}'. "
            f"Available columns: {list(df.columns)}"
        )

    parsed_dates = pd.to_datetime(df[date_col], errors="coerce")
    if parsed_dates.notna().sum() == 0:
        raise ValueError(
            f"{path} has date column '{date_col}', but no parsable dates were found."
        )

    parsed_values = pd.to_numeric(df[expected_value_col], errors="coerce")
    if parsed_values.notna().sum() == 0:
        raise ValueError(
            f"{path} has expected column '{expected_value_col}', "
            "but no numeric values were found."
        )

    return df, date_col


def read_quarterly_series(
    path: Path,
    series_name: str,
    expected_value_col: str,
) -> pd.Series:
    """Read one downloaded FCI CSV as a quarterly series."""
    df, date_col = validate_fci_csv_schema(path, expected_value_col)

    s = pd.Series(
        pd.to_numeric(df[expected_value_col], errors="coerce").to_numpy(),
        index=pd.to_datetime(df[date_col], errors="coerce"),
        name=series_name,
    ).dropna()

    out = s.groupby(s.index.to_period("Q")).mean()
    out.index = pd.PeriodIndex(out.index, freq="Q")
    return out


def build_covid_index(hlw_path: Path, periods: pd.PeriodIndex) -> pd.Series:
    """
    Build the quarterly COVID index used in estimation.

    Live / paper-baseline rule:
    - read sheet 'US input data' from the HLW Excel file
    - use column 6
    - keep the values through 2022Q4
    - set 2023Q1 onward to zero
    - preserve NaN values from the source (do not fill with zero)
    """
    hlw = pd.read_excel(hlw_path, sheet_name="US input data")
    if hlw.shape[1] < 6:
        raise ValueError(
            f"{hlw_path} does not have at least 6 columns in sheet 'US input data'."
        )

    values = pd.to_numeric(hlw.iloc[:, 5], errors="coerce").to_numpy()

    out = pd.Series(0.0, index=periods, name="covid_index")
    target = periods[
        (periods >= pd.Period("1960Q1", freq="Q")) &
        (periods <= pd.Period("2022Q4", freq="Q"))
    ]
    n = min(len(values), len(target))
    out.loc[target[:n]] = values[:n]
    out.loc[out.index >= pd.Period("2023Q1", freq="Q")] = 0.0
    return out


def annualize_discount_rate(s: pd.Series) -> pd.Series:
    return 100.0 * ((1.0 + s / 36000.0) ** 365.0 - 1.0)


def build_covid_dummies(periods: pd.PeriodIndex) -> pd.DataFrame:
    """Construct the 4-column quarterly COVID dummy matrix used in estimation."""
    df = pd.DataFrame(
        {
            "covid_dummy_0": 1.0,
            "covid_dummy_2020": 0.0,
            "covid_dummy_2021": 0.0,
            "covid_dummy_2022": 0.0,
        },
        index=periods,
    )

    df.loc[pd.period_range("2020Q2", "2020Q4", freq="Q"), :] = [0.0, 1.0, 0.0, 0.0]
    df.loc[pd.period_range("2021Q1", "2021Q4", freq="Q"), :] = [0.0, 0.0, 1.0, 0.0]
    df.loc[pd.period_range("2022Q1", "2022Q4", freq="Q"), :] = [0.0, 0.0, 0.0, 1.0]

    return df


def build_data_hlm(
    fred_panel: pd.DataFrame,
    covid_index: pd.Series,
    fci_3yr: pd.Series,
    fci_1yr: pd.Series,
    periods: pd.PeriodIndex,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the main quarterly panel used in estimation.

    This is the Python counterpart of the MATLAB data_hlm construction.
    The output columns match the variables later used by estimate_CCS_struct.
    """
    gdp_tot = np.log(fred_panel["GDPC1"])
    ffr_base = fred_panel["FEDFUNDS"]
    ffr_ny = fred_panel["INTDSRUSM193N"]
    fci = fred_panel["NFCI"]
    gdp_pot_cbo = np.log(fred_panel["GDPPOT"])

    ffr_base_d = annualize_discount_rate(ffr_base)
    ffr_ny_d = annualize_discount_rate(ffr_ny)

    # Use the FRBNY discount rate in the early part of the sample, matching MATLAB.
    ffr = ffr_base_d.copy()
    early_mask = periods <= pd.Period("1964Q4", freq="Q")
    ffr.loc[early_mask] = ffr_ny_d.loc[early_mask]

    pce_core = 400.0 * np.log(fred_panel["PCEPILFE"] / fred_panel["PCEPILFE"].shift(1))
    pce_core.iloc[0] = np.nan

    # MATLAB: movmean(pce_core,[3 0],'includenan')
    expected_pce = pce_core.rolling(window=4, min_periods=4).mean()

    data_hlm = pd.DataFrame(
        {
            "date": quarter_end_dates(periods),
            "date_decimal": quarter_decimal(periods),
            "gdp_tot": gdp_tot.to_numpy(),
            "pce_core": pce_core.to_numpy(),
            "expected_pce": expected_pce.to_numpy(),
            "ffr": ffr.to_numpy(),
            "covid_index": covid_index.reindex(periods).to_numpy(),
            "fci_chicago_fed": fci.to_numpy(),
            "fci_new_3yr": fci_3yr.reindex(periods).to_numpy(),
            "fci_new_1yr": fci_1yr.reindex(periods).to_numpy(),
            "gdp_pot_cbo": gdp_pot_cbo.to_numpy(),
        }
    )

    covid_dummies = build_covid_dummies(periods).copy()
    covid_dummies.insert(0, "date", quarter_end_dates(periods))
    covid_dummies.insert(1, "date_decimal", quarter_decimal(periods))

    return data_hlm, covid_dummies


def trim_to_latest_common_sample(
    fred_panel: pd.DataFrame,
    data_hlm: pd.DataFrame,
    covid_dummies: pd.DataFrame,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Trim the sample to the latest quarter for which all required inputs are available."""
    required_cols = [
        "gdp_tot",
        "pce_core",
        "expected_pce",
        "ffr",
        "covid_index",
        "fci_chicago_fed",
        "fci_new_3yr",
        "fci_new_1yr",
        "gdp_pot_cbo",
    ]

    valid_mask = data_hlm[required_cols].notna().all(axis=1)
    valid_positions = np.flatnonzero(valid_mask.to_numpy())

    if len(valid_positions) == 0:
        raise ValueError("No quarter has all required inputs available.")

    last_valid_pos = int(valid_positions[-1])
    sample_end = data_hlm.loc[last_valid_pos, "date"]
    logger.info("Trimming sample through latest common quarter: %s", sample_end.date())

    data_hlm_trimmed = data_hlm.iloc[: last_valid_pos + 1].reset_index(drop=True)
    covid_dummies_trimmed = covid_dummies.iloc[: last_valid_pos + 1].reset_index(drop=True)
    fred_panel_trimmed = fred_panel.iloc[: last_valid_pos + 1].copy()

    return fred_panel_trimmed, data_hlm_trimmed, covid_dummies_trimmed


def save_intermediate_outputs(
    fred_panel: pd.DataFrame,
    data_hlm: pd.DataFrame,
    covid_dummies: pd.DataFrame,
    out_dir: Path,
) -> None:
    """Save the assembled intermediate files used by the estimation pipeline."""
    out_dir.mkdir(parents=True, exist_ok=True)

    fred_raw = fred_panel.copy()
    fred_raw.insert(0, "date", quarter_end_dates(fred_panel.index))
    fred_raw.insert(1, "date_decimal", quarter_decimal(fred_panel.index))
    fred_raw.to_csv(out_dir / "data_fred_raw.csv", index=False)

    data_hlm.to_csv(out_dir / "data_hlm.csv", index=False)
    covid_dummies.to_csv(out_dir / "covid_dummies.csv", index=False)

    metadata = {
        "series_name_var": [
            "Date",
            "Date Decimal",
            "Real GDP",
            "PCE Core Inflation",
            "Expected Inflation",
            "FFR",
            "Covid Index",
            "FCI - Chicago Fed",
            "New FCI - 3 year lookback",
            "New FCI - 1 year lookback",
            "CBO Potential GDP",
        ]
    }
    with open(out_dir / "assemble_data_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--start-date", default="1959-01-01")
    p.add_argument(
        "--end-date",
        default=None,
        help="Observation end date. Defaults to today's date if omitted.",
    )
    p.add_argument(
        "--assembled-dir",
        default=str(raw_root() / "assembled"),
        help="Directory for intermediate outputs from assemble_data",
    )
    p.add_argument(
        "--use-cached-fred",
        action="store_true",
        help="Reuse existing assembled/data_fred_raw.csv",
    )
    return p.parse_args()


def main() -> None:
    logger = get_logger()
    args = parse_args()

    if args.end_date is None:
        args.end_date = pd.Timestamp.today().strftime("%Y-%m-%d")

    periods = quarter_index(args.start_date, args.end_date)
    assembled_dir = Path(args.assembled_dir)

    hlw_path = (
        raw_root()
        / "HLW_data"
        / "Holston_Laubach_Williams_current_estimates.xlsx"
    )

    fci_dir = raw_root() / "fci_index"
    downloaded_fci = download_new_fci_files(fci_dir, logger)
    fci_3yr_path = downloaded_fci["fci_g_public_quarterly_3yr.csv"]
    fci_1yr_path = downloaded_fci["fci_g_public_quarterly_1yr.csv"]

    fred_cache = assembled_dir / "data_fred_raw.csv"
    if args.use_cached_fred and fred_cache.exists():
        logger.info("Using cached FRED data from %s", fred_cache)
        cached = pd.read_csv(fred_cache, parse_dates=["date"])
        cached["quarter"] = cached["date"].dt.to_period("Q")
        fred_panel = cached.set_index("quarter")[list(FRED_SERIES.keys())]
        fred_panel.index = pd.PeriodIndex(fred_panel.index, freq="Q")
        fred_panel = fred_panel.reindex(periods)
    else:
        api_key = os.getenv("FRED_API_KEY")
        if not api_key:
            raise EnvironmentError("FRED_API_KEY is not set.")
        fred_panel = build_fred_panel(
            FRED_SERIES.keys(),
            args.start_date,
            args.end_date,
            api_key,
            periods,
            logger,
        )

    fci_3yr = read_quarterly_series(
        fci_3yr_path,
        "fci_new_3yr",
        FED_FCI_FILES["fci_g_public_quarterly_3yr.csv"]["expected_value_col"],
    )
    fci_1yr = read_quarterly_series(
        fci_1yr_path,
        "fci_new_1yr",
        FED_FCI_FILES["fci_g_public_quarterly_1yr.csv"]["expected_value_col"],
    )
    covid_index = build_covid_index(hlw_path, periods)

    data_hlm, covid_dummies = build_data_hlm(
        fred_panel=fred_panel,
        covid_index=covid_index,
        fci_3yr=fci_3yr,
        fci_1yr=fci_1yr,
        periods=periods,
    )

    fred_panel, data_hlm, covid_dummies = trim_to_latest_common_sample(
        fred_panel=fred_panel,
        data_hlm=data_hlm,
        covid_dummies=covid_dummies,
        logger=logger,
    )

    save_intermediate_outputs(
        fred_panel=fred_panel,
        data_hlm=data_hlm,
        covid_dummies=covid_dummies,
        out_dir=assembled_dir,
    )

    logger.info("assemble_data completed.")


if __name__ == "__main__":
    main()
