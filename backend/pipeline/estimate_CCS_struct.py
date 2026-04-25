from __future__ import annotations

"""
Step-3 FCI* estimation script.

Purpose
-------
This script runs the step-3 state-space estimation used to construct the
website outputs:
    - FCI*
    - FCI gap
    - output gap

Relative to the full MATLAB script, this version keeps the estimation logic
that is needed for the backend pipeline and drops the paper-specific table
and figure generation code.

Main flow
---------
1. Load assembled inputs created by assemble_data.py.
2. Build the step-3 y_t and x_t matrices used by the Kalman likelihood.
3. Estimate the model by maximum likelihood.
   - pass 1: optimize from several deterministic starting values
   - update P_0 using the pass-1 filtered covariance
   - pass 2: re-optimize from the pass-1 winner
4. Rerun the filter/smoother at the optimum.
5. Write fcistar.csv, metadata.json, theta_opt3.json, and optional SE files.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from statsmodels.tsa.filters.hp_filter import hpfilter

from _aux_fun.likelihood_kalman_covid import likelihood_kalman_covid
from _aux_fun.evaluate_kalman_covid import evaluate_kalman_covid
from _aux_fun.map_theta_to_mats_CCS_struct_step3_backward_inf import (
    map_theta_to_mats_CCS_struct_step3_backward_inf,
)
from _aux_fun.unpack_params_struct_step3 import unpack_params_struct_step3
from _aux_fun.compute_se_kalman import compute_se_kalman


def get_logger() -> logging.Logger:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    return logging.getLogger("estimate_CCS_struct")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def backend_root() -> Path:
    return repo_root() / "backend"


def raw_root() -> Path:
    return backend_root() / "data" / "raw"


def output_root() -> Path:
    return backend_root() / "data" / "output"


def vintages_root() -> Path:
    return output_root() / "vintages"
    

def assembled_root() -> Path:
    return raw_root() / "assembled"


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def find_decimal_idx(date_decimal: np.ndarray, target: float) -> int:
    """Locate the row corresponding to a decimal-date target such as 2024.75."""
    idx = np.where(np.isclose(date_decimal, target))[0]
    if len(idx) == 0:
        raise ValueError(f"Could not find date_decimal == {target}")
    return int(idx[0])


def sample_window(
    recent_data: int,
    sample_end_decimal: float | None = None,
    latest_available_decimal: float | None = None,
) -> tuple[float, float]:
    """
    Replicate the MATLAB sample selection controlled by recent_data.

    For recent_data == 3 (the live-site specification), keep the MATLAB start
    date but let the sample end default to the latest quarter available in
    assembled data_hlm.csv unless the user explicitly passes
    --sample-end-decimal.
    """
    windows = {
        0: (1961.00, 2019.75),
        1: (1990.00, 2019.75),
        2: (1975.00, 2024.75),
        3: (1990.25, 2025.50),
        4: (1961.00, 2024.75),
    }
    if recent_data not in windows:
        raise ValueError(f"Unsupported recent_data={recent_data}")

    start_dec, end_dec = windows[recent_data]

    if sample_end_decimal is not None:
        end_dec = sample_end_decimal
    elif recent_data == 3:
        if latest_available_decimal is None:
            raise ValueError(
                "latest_available_decimal must be provided when recent_data == 3 "
                "and sample_end_decimal is not specified."
            )
        end_dec = latest_available_decimal

    return start_dec, end_dec


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the assembled data matrices produced by assemble_data.py."""
    data_hlm_path = assembled_root() / "data_hlm.csv"
    covid_dummies_path = assembled_root() / "covid_dummies.csv"
    if not data_hlm_path.exists():
        raise FileNotFoundError(f"Missing {data_hlm_path}. Run assemble_data.py first.")
    if not covid_dummies_path.exists():
        raise FileNotFoundError(f"Missing {covid_dummies_path}. Run assemble_data.py first.")
    return (
        pd.read_csv(data_hlm_path, parse_dates=["date"]),
        pd.read_csv(covid_dummies_path, parse_dates=["date"]),
    )


def build_step3_matrices(
    data_hlm: pd.DataFrame,
    covid_dummies_df: pd.DataFrame,
    recent_data: int,
    set_covid_2022: bool,
    sample_end_decimal: float | None = None,
) -> dict[str, Any]:
    """
    Build the step-3 observables and lagged regressors.

    MATLAB backward_inf == 1 uses:
        y_t = [GDP_t, core inflation_t]
        x_t = [y_{t-1}, FCI_{t-1}, pi_{t-1}, pi_{t-2:4}, d_t, d_{t-1}]

    Notes
    -----
    - GDP is scaled by 100 exactly as in MATLAB.
    - The HP filter is run on the full sample before trimming.
    - If set_covid_2022 == 1, the COVID index is forced to zero from 2023Q1 onward.
    """
    date = data_hlm["date"].copy()
    date_decimal = data_hlm["date_decimal"].to_numpy()

    gdp_tot = 100.0 * data_hlm["gdp_tot"].to_numpy()
    pce_core = data_hlm["pce_core"].to_numpy()
    covid_index = data_hlm["covid_index"].to_numpy().copy()
    if set_covid_2022:
        covid_index[date_decimal >= 2023.0] = 0.0

    fci_chicago = data_hlm["fci_chicago_fed"].to_numpy()
    fci_new_3yr = data_hlm["fci_new_3yr"].to_numpy()
    fci_use = fci_new_3yr if recent_data == 3 else fci_chicago

    latest_available_decimal = float(data_hlm["date_decimal"].dropna().max())

    start_dec, end_dec = sample_window(
    recent_data=recent_data,
    sample_end_decimal=sample_end_decimal,
    latest_available_decimal=latest_available_decimal,
    )
    start_idx = find_decimal_idx(date_decimal, start_dec)
    end_idx = find_decimal_idx(date_decimal, end_dec)

    _, trend = hpfilter(gdp_tot, lamb=36000)
    nat_output = np.asarray(trend)

    y_mat_aux = np.column_stack([gdp_tot, pce_core])

    # x_t = [y_{t-1}, FCI_{t-1}, pi_{t-1}, pi_{t-2:4}, d_t, d_{t-1}]
    x_mat_step3_aux = np.full((len(data_hlm), 6), np.nan)
    x_mat_step3_aux[1:, 0] = gdp_tot[:-1]
    x_mat_step3_aux[1:, 1] = fci_use[:-1]
    x_mat_step3_aux[1:, 2] = pce_core[:-1]

    # pi_{t-2:4} in the MATLAB code is a shrinking 3-quarter moving average
    # built from inflation up to t-2.
    past_inf_proxy = pd.Series(pce_core[:-2]).rolling(window=3, min_periods=1).mean().to_numpy()
    x_mat_step3_aux[2:, 3] = past_inf_proxy
    x_mat_step3_aux[:, 4] = covid_index
    x_mat_step3_aux[1:, 5] = covid_index[:-1]

    y_mat = y_mat_aux[start_idx : end_idx + 1, :]
    x_mat_step3 = x_mat_step3_aux[start_idx : end_idx + 1, :]

    covid_cols = ["covid_dummy_0", "covid_dummy_2020", "covid_dummy_2021", "covid_dummy_2022"]
    covid_dummies = covid_dummies_df.iloc[start_idx : end_idx + 1][covid_cols].to_numpy()

    # MATLAB initializes potential-output states using the 4 pre-sample HP-filtered values.
    output_nat_init = nat_output[start_idx - 1 : start_idx - 5 : -1]
    if output_nat_init.shape[0] != 4:
        raise ValueError("Failed to construct output_nat_init with 4 lags.")

    return {
        "date": date,
        "date_decimal": date_decimal,
        "gdp_tot": gdp_tot,
        "fci_use": fci_use,
        "nat_output": nat_output,
        "y_mat": y_mat,
        "x_mat_step3": x_mat_step3,
        "covid_dummies": covid_dummies,
        "output_nat_init": output_nat_init,
        "start_idx": start_idx,
        "end_idx": end_idx,
    }


def step3_initial_conditions(output_nat_init: np.ndarray, nat_output: np.ndarray, start_idx: int) -> np.ndarray:
    """
    Initial state vector for step 3.

    State order:
        [y*_t, y*_{t-1}, y*_{t-2},
         g_t, g_{t-1}, g_{t-2},
         FCI*_t, FCI*_{t-1},
         delta_t, delta_{t-1}, delta_{t-2}]
    """
    xi_0 = np.zeros(11)
    xi_0[0:3] = output_nat_init[0:3]
    xi_0[3] = nat_output[start_idx - 1] - nat_output[start_idx - 2]
    xi_0[4] = nat_output[start_idx - 2] - nat_output[start_idx - 3]
    xi_0[5] = nat_output[start_idx - 3] - nat_output[start_idx - 4]
    xi_0[6:8] = np.array([0.0, 0.0])
    xi_0[8:11] = np.array([0.0, 0.0, 0.0])
    return xi_0


def x_to_z(x: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> np.ndarray:
    """Map bounded parameters x to unconstrained coordinates z."""
    x = np.asarray(x, dtype=float)
    lb = np.asarray(lb, dtype=float)
    ub = np.asarray(ub, dtype=float)
    z = np.empty_like(x, dtype=float)
    for i in range(len(x)):
        lbi, ubi, xi = lb[i], ub[i], x[i]
        if np.isclose(lbi, ubi):
            z[i] = 0.0
        elif np.isfinite(lbi) and np.isfinite(ubi):
            s = 2.0 * (xi - lbi) / (ubi - lbi) - 1.0
            z[i] = np.arcsin(np.clip(s, -1.0, 1.0))
        elif np.isfinite(lbi) and not np.isfinite(ubi):
            z[i] = np.sqrt(max(xi - lbi, 0.0))
        elif not np.isfinite(lbi) and np.isfinite(ubi):
            z[i] = np.sqrt(max(ubi - xi, 0.0))
        else:
            z[i] = xi
    return z


def z_to_x(z: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> np.ndarray:
    """Inverse map from unconstrained coordinates z back to bounded x."""
    z = np.asarray(z, dtype=float)
    lb = np.asarray(lb, dtype=float)
    ub = np.asarray(ub, dtype=float)
    x = np.empty_like(z, dtype=float)
    for i in range(len(z)):
        lbi, ubi, zi = lb[i], ub[i], z[i]
        if np.isclose(lbi, ubi):
            x[i] = lbi
        elif np.isfinite(lbi) and np.isfinite(ubi):
            x[i] = lbi + (np.sin(zi) + 1.0) * (ubi - lbi) / 2.0
            x[i] = min(max(x[i], lbi), ubi)
        elif np.isfinite(lbi) and not np.isfinite(ubi):
            x[i] = lbi + zi**2
        elif not np.isfinite(lbi) and np.isfinite(ubi):
            x[i] = ubi - zi**2
        else:
            x[i] = zi
    return x


def fminsearchbnd_like(
    objective_x: Callable[[np.ndarray], float],
    x0: np.ndarray,
    lb: np.ndarray,
    ub: np.ndarray,
    maxiter: int = 20000,
    maxfev: int = 20000,
    logger: logging.Logger | None = None,
) -> tuple[np.ndarray, float, Any]:
    """
    Small Python analogue of MATLAB's fminsearchbnd:
    optimize in transformed coordinates, then map back to bounded parameters.
    """
    z0 = x_to_z(np.asarray(x0, float), np.asarray(lb, float), np.asarray(ub, float))

    def objective_z(z: np.ndarray) -> float:
        x = z_to_x(z, lb, ub)
        val = float(objective_x(x))
        if not np.isfinite(val):
            return 1.0e12
        return val

    if logger is not None:
        logger.info("Starting fminsearchbnd-like Nelder-Mead optimization")

    res = minimize(
        objective_z,
        z0,
        method="Nelder-Mead",
        options={"maxiter": maxiter, "maxfev": maxfev, "disp": True, "xatol": 1e-8, "fatol": 1e-8, "adaptive": False},
    )
    xopt = z_to_x(res.x, lb, ub)
    return xopt, float(objective_x(xopt)), res


def optimize_with_fixed_bounds(
    objective_full: Callable[[np.ndarray], float],
    x0_full: np.ndarray,
    lb_full: np.ndarray,
    ub_full: np.ndarray,
    logger: logging.Logger,
) -> tuple[np.ndarray, float, Any]:
    """
    Optimize only over the free parameters.

    Parameters with LB == UB are treated as fixed and are reinserted into
    the full theta vector before each objective evaluation.
    """
    free_mask = ~np.isclose(lb_full, ub_full)
    fixed_template = np.asarray(x0_full, float).copy()
    fixed_template[~free_mask] = lb_full[~free_mask]

    def objective_free(theta_free: np.ndarray) -> float:
        theta_full = fixed_template.copy()
        theta_full[free_mask] = theta_free
        return objective_full(theta_full)

    theta_free_opt, fval, res = fminsearchbnd_like(
        objective_x=objective_free,
        x0=np.asarray(x0_full, float)[free_mask],
        lb=np.asarray(lb_full, float)[free_mask],
        ub=np.asarray(ub_full, float)[free_mask],
        maxiter=20000,
        maxfev=20000,
        logger=logger,
    )
    theta_full_opt = fixed_template.copy()
    theta_full_opt[free_mask] = theta_free_opt
    return theta_full_opt, fval, res


def clip_to_bounds(x: np.ndarray, lb: np.ndarray, ub: np.ndarray) -> np.ndarray:
    """Clip a candidate starting value to the admissible box constraints."""
    x = np.asarray(x, dtype=float).copy()
    for i in range(len(x)):
        if np.isfinite(lb[i]):
            x[i] = max(x[i], lb[i])
        if np.isfinite(ub[i]):
            x[i] = min(x[i], ub[i])
    return x


def build_step3_spec(restrict_af: int, set_covid_2022: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Step-3 parameter specification.

    Parameter order:
        0  a_y1
        1  a_y2
        2  b_pi
        3  b_y
        4  sigma_yt
        5  sigma_pi
        6  sigma_ys
        7  sigma_pi_e
        8  sigma_delta
        9  phi
        10 kappa_2020
        11 kappa_2021
        12 kappa_2022
        13 rho_pi
        14 rho_delta
        15 eta
        16 theta17   (the parameter used to construct a_f)

    The arrays x0_step3, LB_step3, UB_step3, and fixed_value_step3 follow
    the MATLAB step-3 setup. Fixed parameters are imposed by setting LB = UB
    = fixed_value.
    """
    x0_step3 = np.array(
        [
            0.3,    # a_y1
            0.0,    # a_y2
            0.5,    # b_pi
            0.1,    # b_y
            0.4140, # sigma_yt
            0.5011, # sigma_pi
            0.2,    # sigma_ys
            0.2,    # sigma_pi_e
            0.5,    # sigma_delta
            -0.1,   # phi
            7.0,    # kappa_2020
            3.0,    # kappa_2021
            2.0,    # kappa_2022
            0.9,    # rho_pi
            0.7,    # rho_delta
            0.8,    # eta
            1.0,    # theta17
        ],
        dtype=float,
    )

    if set_covid_2022:
        fix_params_step3 = np.array([1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, restrict_af], dtype=int)
    else:
        fix_params_step3 = np.array([1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, restrict_af], dtype=int)

    lb_step3 = np.array([-2.0, -2.0, 0.0, 0.025, 0.01, 0.00001, 0.0001, 0.00001, 0.001, -np.inf, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=float)
    ub_step3 = np.array([2.0, 2.0, 1.0, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1.0, 1.0, 1.0, 5.0], dtype=float)
    fixed_value_step3 = np.array([0.0, 0.0, 0.689, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.7, 0.9, 1.0], dtype=float)

    for i in range(len(x0_step3)):
        if fix_params_step3[i] == 1:
            lb_step3[i] = fixed_value_step3[i]
            ub_step3[i] = fixed_value_step3[i]
            x0_step3[i] = fixed_value_step3[i]

    return x0_step3, lb_step3, ub_step3, fix_params_step3


# Multistart note:
# Step 3 uses four deterministic starting values to reduce sensitivity to local optima.
# Start 1 is the original MATLAB baseline x0_step3.
# Starts 2-4 keep the same fixed parameters and perturb only the free
# parameters around the MATLAB baseline, mainly along variance and
# persistence dimensions.
def build_multistart_starts(
    x0_full: np.ndarray,
    lb_full: np.ndarray,
    ub_full: np.ndarray,
    set_covid_2022: bool,
    n_starts: int = 4,
) -> list[np.ndarray]:
    base = np.asarray(x0_full, dtype=float)
    lb = np.asarray(lb_full, dtype=float)
    ub = np.asarray(ub_full, dtype=float)
    fixed_mask = np.isclose(lb, ub)

    starts: list[np.ndarray] = [base.copy()]
    sigma_idx = [4, 5, 6, 8]
    kappa_free_idx = [10] if set_covid_2022 else [10, 11, 12]

    # Start 2: lower variances and lower persistence.
    x2 = base.copy()
    for j in sigma_idx:
        x2[j] = 0.5 * base[j]
    x2[14] = 0.4
    x2[15] = 0.4
    for j in kappa_free_idx:
        x2[j] = max(1.5, 0.5 * base[j])
    x2 = clip_to_bounds(x2, lb, ub)
    x2[fixed_mask] = lb[fixed_mask]
    starts.append(x2)

    # Start 3: higher variances and medium persistence.
    x3 = base.copy()
    for j in sigma_idx:
        x3[j] = 1.5 * base[j]
    x3[14] = 0.7
    x3[15] = 0.7
    for j in kappa_free_idx:
        x3[j] = 1.5 * base[j]
    x3 = clip_to_bounds(x3, lb, ub)
    x3[fixed_mask] = lb[fixed_mask]
    starts.append(x3)

    # Start 4: high persistence with the other parameters near baseline.
    x4 = base.copy()
    x4[14] = 0.9
    x4[15] = 0.9
    x4 = clip_to_bounds(x4, lb, ub)
    x4[fixed_mask] = lb[fixed_mask]
    starts.append(x4)

    return starts[: max(1, n_starts)]


def build_se_constraints(lb_step3: np.ndarray, ub_step3: np.ndarray) -> dict[str, Any]:
    """Constraints used later by compute_se_kalman / compute_opg."""
    def const_f(theta: np.ndarray) -> np.ndarray:
        theta = np.asarray(theta, dtype=float).reshape(-1)
        val = (
            (theta[0] + theta[1] >= 1.0)
            + (theta[15] >= 1.0)
            + (theta[14] >= 1.0)
            + np.sum(theta[4:9] < 0.0)
        )
        return np.array([val], dtype=float)

    return {"LB": np.asarray(lb_step3, dtype=float), "UB": np.asarray(ub_step3, dtype=float), "const_f": const_f}


def estimate_step3(
    y_mat: np.ndarray,
    x_mat_step3: np.ndarray,
    covid_dummies: np.ndarray,
    output_nat_init: np.ndarray,
    nat_output: np.ndarray,
    start_idx: int,
    restrict_af: int,
    set_covid_2022: bool,
    logger: logging.Logger,
    n_starts: int = 4,
    compute_se: bool = False,
    se_draws: int = 5000,
) -> tuple[np.ndarray, Any, float, dict[str, Any] | None]:
    """
    Estimate the step-3 state-space model.

    MATLAB-style two-pass procedure:
        pass 1: optimize with the initial P_0
        update P_0 using the filtered covariance at the pass-1 optimum
        pass 2: re-optimize from the pass-1 winner

    The multistart extension affects only the pass-1 starting values.
    """
    lambda_g = 0.0667
    settings_step3 = {
        "detrend_y": 0,
        "smoother": 0,
        "states_smoother": list(range(1, 12)),
        "xi_0": step3_initial_conditions(output_nat_init, nat_output, start_idx),
        "P_0": 0.25 * np.eye(11),
        "fixed_param": [lambda_g],
        "covid_dummies": covid_dummies,
    }
    settings_step3["P_0"][0:3, 0:3] = np.eye(3)

    x0_step3, lb_step3, ub_step3, _ = build_step3_spec(restrict_af=restrict_af, set_covid_2022=set_covid_2022)

    def objective_full(theta_full: np.ndarray) -> float:
        # scipy minimizes, so we use the negative log-likelihood.
        return -likelihood_kalman_covid(
            y_mat,
            x_mat_step3,
            theta_full,
            map_theta_to_mats_CCS_struct_step3_backward_inf,
            settings_step3,
        )

    start_points = build_multistart_starts(x0_step3, lb_step3, ub_step3, set_covid_2022, n_starts=n_starts)
    logger.info("Starting step 3 optimization (pass 1, %d starts)", len(start_points))

    best_theta_aux = None
    best_obj_aux = np.inf
    for i, theta_start in enumerate(start_points, start=1):
        logger.info("Pass 1 start %d/%d | objective at start = %.6f", i, len(start_points), objective_full(theta_start))
        theta_opt_i, fval_i, _ = optimize_with_fixed_bounds(
            objective_full=objective_full, x0_full=theta_start, lb_full=lb_step3, ub_full=ub_step3, logger=logger
        )
        logger.info("Pass 1 result %d/%d | objective = %.6f", i, len(start_points), fval_i)
        if fval_i < best_obj_aux:
            best_obj_aux = fval_i
            best_theta_aux = theta_opt_i

    if best_theta_aux is None:
        raise RuntimeError("Pass 1 failed to produce a candidate optimum.")

    # Update P_0 from the pass-1 filtered covariance, exactly in the spirit
    # of the MATLAB step-3 rerun.
    theta_opt3_aux = best_theta_aux
    estimates_step3_aux = evaluate_kalman_covid(
        y_mat, x_mat_step3, theta_opt3_aux, map_theta_to_mats_CCS_struct_step3_backward_inf, settings_step3
    )
    settings_step3["P_0"] = np.asarray(_get(estimates_step3_aux, "state_contemporaneous_pred")[0])

    logger.info("Starting step 3 optimization (pass 2)")
    theta_opt3, fval, _ = optimize_with_fixed_bounds(
        objective_full=objective_full, x0_full=theta_opt3_aux, lb_full=lb_step3, ub_full=ub_step3, logger=logger
    )

    settings_final = dict(settings_step3)
    settings_final["smoother"] = 1
    estimates_step3 = evaluate_kalman_covid(
        y_mat, x_mat_step3, theta_opt3, map_theta_to_mats_CCS_struct_step3_backward_inf, settings_final
    )

    se_output = None
    if compute_se:
        logger.info("Computing parameter and state standard errors")
        constraints = build_se_constraints(lb_step3, ub_step3)

        def fun_se(theta: np.ndarray):
            return evaluate_kalman_covid(
                y_mat, x_mat_step3, theta, map_theta_to_mats_CCS_struct_step3_backward_inf, settings_final
            )

        se_output = compute_se_kalman(theta_opt=theta_opt3, constraints=constraints, evaluate_kalman=fun_se, N_draws=se_draws)

    return theta_opt3, estimates_step3, fval, se_output


def build_output_series(theta_opt3: np.ndarray, estimates_step3: Any, inputs: dict[str, Any]) -> pd.DataFrame:
    """
    Transform the estimated states into the backend output contract:
        date, fci, fcistar, fci_gap, y_gap

    State positions are written first in MATLAB's 1-based notation and then
    translated where needed to Python indexing.
    """
    param = unpack_params_struct_step3(theta_opt3)
    start_idx = inputs["start_idx"]
    end_idx = inputs["end_idx"]
    x_mat_step3 = inputs["x_mat_step3"]
    gdp_tot = inputs["gdp_tot"]
    fci_use = inputs["fci_use"]
    date = inputs["date"]

    y_n_pos = 1
    g_pos = 4
    fci_pos = 7
    delta_pos = 9
    covid_d_pos_step3 = (4, 5)

    state_contemp = np.asarray(_get(estimates_step3, "state_contemporaneous"))
    phi = float(param["phi"])
    rho_delta = float(param["rho_delta"])

    fci_star = state_contemp[:, fci_pos - 1]
    y_star = (
        state_contemp[:, y_n_pos]
        + state_contemp[:, g_pos]
        + state_contemp[:, delta_pos - 1]
        - rho_delta * state_contemp[:, delta_pos]
        + phi * x_mat_step3[:, covid_d_pos_step3[0]]
    )
    y_gap = gdp_tot[start_idx : end_idx + 1] - y_star
    fci_gap = fci_use[start_idx : end_idx + 1] - fci_star

    return pd.DataFrame(
        {
            "date": pd.to_datetime(date.iloc[start_idx : end_idx + 1]).dt.strftime("%Y-%m-%d"),
            "fci": fci_use[start_idx : end_idx + 1],
            "fcistar": fci_star,
            "fci_gap": fci_gap,
            "y_gap": y_gap,
        }
    )


def save_se_outputs(se_output: dict[str, Any] | None) -> None:
    """Write the parameter and state SE summaries used for diagnostics."""
    if se_output is None:
        return
    out_dir = output_root()
    out_dir.mkdir(parents=True, exist_ok=True)

    se_vec = np.asarray(se_output["se_vec"], dtype=float)
    t_stat_vec = np.asarray(se_output["t_stat_vec"], dtype=float)
    used_draws = int(np.asarray(se_output["used_draws"]).reshape(-1)[0])

    theta_names = [
        "a_y1", "a_y2", "b_pi", "b_y", "sigma_yt", "sigma_pi", "sigma_ys", "sigma_pi_e",
        "sigma_delta", "phi", "kappa_2020", "kappa_2021", "kappa_2022", "rho_pi", "rho_delta", "eta", "theta17",
    ]
    payload = {"used_draws": used_draws, "parameters": []}
    for i, name in enumerate(theta_names):
        payload["parameters"].append(
            {
                "index": i + 1,
                "name": name,
                "se": None if np.isnan(se_vec[i]) else float(se_vec[i]),
                "t_stat": None if np.isnan(t_stat_vec[i]) else float(t_stat_vec[i]),
            }
        )

    with open(out_dir / "theta_se.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    pd.DataFrame(np.asarray(se_output["state_ses"], dtype=float)).to_csv(out_dir / "state_ses.csv", index=False)


def write_outputs(fcistar_df: pd.DataFrame, theta_opt3: np.ndarray) -> None:
    """Write current frontend files and archive dated vintage copies."""
    out_dir = output_root()
    vint_dir = vintages_root()

    out_dir.mkdir(parents=True, exist_ok=True)
    vint_dir.mkdir(parents=True, exist_ok=True)

    run_date = pd.Timestamp.today().strftime("%Y-%m-%d")

    metadata = {
        "last_updated": run_date,
        "sample_start": str(fcistar_df["date"].iloc[0]),
        "sample_end": str(fcistar_df["date"].iloc[-1]),
    }

    # Current files read by the frontend
    fcistar_df.to_csv(out_dir / "fcistar.csv", index=False)
    with open(out_dir / "metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)
    with open(out_dir / "theta_opt3.json", "w", encoding="utf-8") as f:
        json.dump({"theta_opt3": theta_opt3.tolist()}, f, indent=2)

    # Archived vintage files keyed by pipeline run date
    vintage_fcistar_path = vint_dir / f"fcistar_{run_date}.csv"
    vintage_metadata_path = vint_dir / f"metadata_{run_date}.json"

    if vintage_fcistar_path.exists() or vintage_metadata_path.exists():
        raise FileExistsError(
            f"Vintage files for run date {run_date} already exist in {vint_dir}. "
            "Refusing to overwrite archived vintages."
        )

    fcistar_df.to_csv(vintage_fcistar_path, index=False, mode="x")
    with open(vintage_metadata_path, "x", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--recent-data", type=int, default=3)
    p.add_argument("--set-covid-2022", type=int, choices=[0, 1], default=1)
    p.add_argument("--restrict-af", type=int, default=1)
    p.add_argument("--sample-end-decimal", type=float, default=None)
    p.add_argument("--n-starts", type=int, default=4)
    p.add_argument("--compute-se", type=int, choices=[0, 1], default=0)
    p.add_argument("--se-draws", type=int, default=5000)
    return p.parse_args()


def main() -> None:
    logger = get_logger()
    args = parse_args()
    logger.info("Loading assembled inputs")
    data_hlm, covid_dummies = load_inputs()

    logger.info("Building step 3 matrices")
    inputs = build_step3_matrices(
        data_hlm=data_hlm,
        covid_dummies_df=covid_dummies,
        recent_data=args.recent_data,
        set_covid_2022=bool(args.set_covid_2022),
        sample_end_decimal=args.sample_end_decimal,
    )

    logger.info("Running step 3 estimation")
    theta_opt3, estimates_step3, fval, se_output = estimate_step3(
        y_mat=inputs["y_mat"],
        x_mat_step3=inputs["x_mat_step3"],
        covid_dummies=inputs["covid_dummies"],
        output_nat_init=inputs["output_nat_init"],
        nat_output=inputs["nat_output"],
        start_idx=inputs["start_idx"],
        restrict_af=args.restrict_af,
        set_covid_2022=bool(args.set_covid_2022),
        logger=logger,
        n_starts=args.n_starts,
        compute_se=bool(args.compute_se),
        se_draws=args.se_draws,
    )

    logger.info("Building website output series")
    fcistar_df = build_output_series(theta_opt3=theta_opt3, estimates_step3=estimates_step3, inputs=inputs)

    logger.info("Writing backend/data/output files")
    write_outputs(fcistar_df, theta_opt3)

    if se_output is not None:
        logger.info("Writing SE outputs")
        save_se_outputs(se_output)

    logger.info("estimate_CCS_struct completed successfully. Final objective = %.6f", fval)


if __name__ == "__main__":
    main()
