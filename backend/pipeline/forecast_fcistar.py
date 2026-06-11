#!/usr/bin/env python3
"""
*******************************************************************************
**** One-quarter-ahead FCI* forecast.
****
**** Extends FCI* one quarter past the last estimated quarter by running the
**** Kalman filter ONE more step with the SAME fixed parameters -- we do NOT
**** re-estimate (that happens on the monthly cron when actual data prints).
****
**** Inputs:
****   ../data/output/theta_opt3.json                 fixed parameter vector
****   ../data/output/forecast_inputs/                 committed by the monthly
****       data_hlm.csv, covid_dummies.csv,            estimation run (Phase 2):
****       filter_settings.json (xi_0, P_0, spec)      the panel + exact filter
****                                                   initial conditions.
****   data/spf/output/spf_forecasts.csv               SPF medians (get_spf.py).
****   ../data/output/fcistar.csv                      committed in-sample FCI*,
****                                                   used as a correctness gate.
****
**** Method:
****   target quarter Q = last estimated quarter + 1.
****   Append one quarter to data_hlm using SPF, treated as if realized:
****       gdp_tot_Q = gdp_tot_{Q-1} + drgdp_Q / 400   (chain SPF growth onto
****                                                     the actual prior level)
****       pce_core_Q = corepce_Q                       (annualized %, as-is)
****   FCI enters the measurement lagged, so the filter step needs no FCI_Q --
****   the lagged regressor FCI_{Q-1} is the last actual quarter's FCI.
****   Rebuild the step-3 matrices over the extended sample and run the filter
****   with the committed xi_0 / P_0. The filtered state at Q gives FCI*_Q.
****
**** Because filtering is causal and xi_0 / P_0 are anchored at the sample
**** start, the in-sample filtered states are reproduced EXACTLY; only the
**** quarter Q is new. We assert that reproduction as a correctness gate.
****
**** Output:
****   ../data/output/fcistar_forecast.csv  one row: the forecast quarter with
****       target_quarter, date, fcistar, y_gap, the SPF inputs, and the survey
****       vintage. (The forecast FCI gap, which needs the daily FCI nowcast for
****       FCI_Q, is added in the next step.)
*******************************************************************************
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from _aux_fun.evaluate_kalman_covid import evaluate_kalman_covid
from _aux_fun.map_theta_to_mats_CCS_struct_step3_backward_inf import (
    map_theta_to_mats_CCS_struct_step3_backward_inf,
)
from estimate_CCS_struct import build_output_series, build_step3_matrices, output_root
from get_spf import latest_forecast_for

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger("forecast_fcistar")

HERE = Path(__file__).resolve().parent
OUT_DIR = output_root()
FINPUTS = OUT_DIR / "forecast_inputs"
SPF_CSV = HERE / "data" / "spf" / "output" / "spf_forecasts.csv"

REPRO_TOL = 1e-6  # max allowed |forecast in-sample FCI* - committed FCI*|


def next_quarter(year: int, q: int) -> tuple[int, int]:
    idx = year * 4 + (q - 1) + 1
    y, r = divmod(idx, 4)
    return y, r + 1


def quarter_end_timestamp(year: int, q: int) -> pd.Timestamp:
    return pd.Period(f"{year}Q{q}", freq="Q").to_timestamp(how="end").normalize()


def quarter_decimal(year: int, q: int) -> float:
    return year + 0.25 * (q - 1)


def load_settings() -> dict:
    return json.loads((FINPUTS / "filter_settings.json").read_text(encoding="utf-8"))


def load_theta() -> np.ndarray:
    payload = json.loads((OUT_DIR / "theta_opt3.json").read_text(encoding="utf-8"))
    return np.asarray(payload["theta_opt3"], dtype=float)


def append_forecast_quarter(
    data_hlm: pd.DataFrame,
    covid_dummies: pd.DataFrame,
    target_year: int,
    target_q: int,
    drgdp: float,
    corepce: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Append one quarter (the SPF forecast) to the panel and covid dummies."""
    ts = quarter_end_timestamp(target_year, target_q)
    dec = quarter_decimal(target_year, target_q)

    gdp_prev = float(data_hlm["gdp_tot"].iloc[-1])
    gdp_new = gdp_prev + drgdp / 400.0  # SPF annualized % growth -> quarterly log change

    # Only the columns the filter uses need values; the rest stay NaN. The
    # contemporaneous FCI (fci_new_3yr) is NOT used by the filter (FCI enters
    # lagged) -- it is left NaN here and supplied later for the gap.
    new_row = {col: np.nan for col in data_hlm.columns}
    new_row.update(
        {"date": ts, "date_decimal": dec, "gdp_tot": gdp_new, "pce_core": corepce, "covid_index": 0.0}
    )
    data_ext = pd.concat([data_hlm, pd.DataFrame([new_row])], ignore_index=True)

    # Copy the last (normal-period) covid-dummy pattern forward.
    cov_row = covid_dummies.iloc[-1].copy()
    cov_row["date"] = ts
    cov_row["date_decimal"] = dec
    covid_ext = pd.concat([covid_dummies, pd.DataFrame([cov_row])], ignore_index=True)

    return data_ext, covid_ext


def run() -> pd.DataFrame:
    settings = load_settings()
    theta = load_theta()

    data_hlm = pd.read_csv(FINPUTS / "data_hlm.csv", parse_dates=["date"])
    covid_dummies = pd.read_csv(FINPUTS / "covid_dummies.csv", parse_dates=["date"])

    last_q = pd.Timestamp(settings["last_quarter"])
    ly, lq = last_q.year, (last_q.month - 1) // 3 + 1
    ty, tq = next_quarter(ly, lq)
    target_label = f"{ty}Q{tq}"
    logger.info("Last estimated quarter %dQ%d -> forecasting %s", ly, lq, target_label)

    panel = pd.read_csv(SPF_CSV)
    fc = latest_forecast_for(panel, ty, tq)
    if fc is None or fc["drgdp"] is None or fc["corepce"] is None:
        raise RuntimeError(f"No SPF forecast available for {target_label}.")
    logger.info(
        "SPF for %s: drgdp=%.4f corepce=%.4f (from %s survey, horizon %d)",
        target_label, fc["drgdp"], fc["corepce"], fc["survey_quarter"], fc["horizon"],
    )

    data_ext, covid_ext = append_forecast_quarter(
        data_hlm, covid_dummies, ty, tq, fc["drgdp"], fc["corepce"]
    )

    # Rebuild the step-3 matrices. recent_data==3 auto-extends the sample end to
    # the latest available quarter, which is now the appended forecast quarter.
    inputs = build_step3_matrices(
        data_hlm=data_ext,
        covid_dummies_df=covid_ext,
        recent_data=int(settings["recent_data"]),
        set_covid_2022=bool(settings["set_covid_2022"]),
        sample_end_decimal=None,
    )

    filter_settings = {
        "detrend_y": int(settings["detrend_y"]),
        "smoother": 0,  # the forecast only needs the filtered (one-sided) state
        "states_smoother": settings["states_smoother"],
        "xi_0": np.asarray(settings["xi_0"], dtype=float),
        "P_0": np.asarray(settings["P_0"], dtype=float),
        "fixed_param": np.asarray(settings["fixed_param"], dtype=float),
        "covid_dummies": inputs["covid_dummies"],
    }

    # The filter's covariance recursion emits benign floating-point warnings on
    # the COVID quarters (huge measurement variance) -- they occur identically in
    # the in-sample estimation and don't affect the (finite) state estimates.
    # Silence them here so the daily forecast log stays readable.
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        estimates = evaluate_kalman_covid(
            inputs["y_mat"],
            inputs["x_mat_step3"],
            theta,
            map_theta_to_mats_CCS_struct_step3_backward_inf,
            filter_settings,
        )

    df = build_output_series(theta_opt3=theta, estimates_step3=estimates, inputs=inputs)

    # Correctness gate: the in-sample part must reproduce the committed FCI*.
    committed = pd.read_csv(OUT_DIR / "fcistar.csv")
    insample = df.iloc[:-1].reset_index(drop=True)
    if len(insample) != len(committed):
        raise RuntimeError(
            f"In-sample length {len(insample)} != committed {len(committed)}; "
            "the appended quarter changed the sample window unexpectedly."
        )
    max_diff = float(np.max(np.abs(insample["fcistar"].to_numpy() - committed["fcistar"].to_numpy())))
    logger.info("In-sample FCI* reproduction: max |diff| = %.3e", max_diff)
    if max_diff > REPRO_TOL:
        raise RuntimeError(
            f"In-sample FCI* does not reproduce the committed series "
            f"(max diff {max_diff:.3e} > tol {REPRO_TOL:.0e})."
        )

    frow = df.iloc[-1]
    forecast = pd.DataFrame(
        [
            {
                "target_quarter": target_label,
                "date": frow["date"],
                "fcistar": round(float(frow["fcistar"]), 6),
                "y_gap": round(float(frow["y_gap"]), 6),
                "drgdp": fc["drgdp"],
                "corepce": fc["corepce"],
                "survey_quarter": fc["survey_quarter"],
                "horizon": fc["horizon"],
            }
        ]
    )
    return forecast


def main() -> None:
    forecast = run()
    out_path = OUT_DIR / "fcistar_forecast.csv"
    forecast.to_csv(out_path, index=False)
    row = forecast.iloc[0]
    logger.info("Wrote %s", out_path)
    logger.info(
        "FORECAST %s (%s): FCI* = %.4f, y_gap = %.4f  [SPF %s, h=%d]",
        row["target_quarter"], row["date"], row["fcistar"], row["y_gap"],
        row["survey_quarter"], int(row["horizon"]),
    )


if __name__ == "__main__":
    main()
