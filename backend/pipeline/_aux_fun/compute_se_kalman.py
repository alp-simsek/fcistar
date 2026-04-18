from __future__ import annotations

from typing import Any, Callable

import numpy as np

from _aux_fun.compute_opg import compute_opg


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def compute_se_kalman(
    theta_opt: np.ndarray,
    constraints,
    evaluate_kalman: Callable[[np.ndarray], Any],
    N_draws: int,
) -> dict[str, np.ndarray]:
    """
    Compute parameter and state standard errors.

    The procedure follows the MATLAB code:
    - use the OPG to compute the covariance matrix of the MLE
    - draw parameters around theta_opt
    - rerun the Kalman filter/smoother at each accepted draw
    - combine filter uncertainty and parameter uncertainty for the states
    """

    theta_opt = np.asarray(theta_opt, dtype=float).reshape(-1)

    UB = np.asarray(_get(constraints, "UB"), dtype=float).reshape(-1)
    LB = np.asarray(_get(constraints, "LB"), dtype=float).reshape(-1)

    try:
        const_f = _get(constraints, "const_f")
    except Exception:
        const_f = lambda x: np.array([-10.0], dtype=float)

    cons_vec = (theta_opt == LB).astype(int) + (theta_opt == UB).astype(int)

    unconst_cols = np.where(cons_vec == 0)[0]
    const_cols = np.where(cons_vec == 1)[0]
    N_param = len(theta_opt)

    out_base = evaluate_kalman(theta_opt)
    state_base = np.asarray(_get(out_base, "state_contemporaneous"), dtype=float)
    T, N_states = state_base.shape

    # Variance covariance matrix ONLY for the unconstrained parameters.
    opg_out = compute_opg(theta_opt, constraints, evaluate_kalman)
    if isinstance(opg_out, dict):
        V_mat = np.asarray(opg_out["V_mat"], dtype=float)
    else:
        V_mat = np.asarray(opg_out, dtype=float)

    # MATLAB: C_mat_aux = chol(V_mat,'lower')
    C_mat_aux = np.linalg.cholesky(V_mat)

    C_mat = np.zeros((N_param, N_param), dtype=float)
    C_mat[np.ix_(unconst_cols, unconst_cols)] = C_mat_aux

    se_vec = np.full(N_param, np.nan, dtype=float)
    if len(unconst_cols) > 0:
        se_vec[unconst_cols] = np.sqrt(np.diag(V_mat))

    t_stat_vec = np.full(N_param, np.nan, dtype=float)
    if len(unconst_cols) > 0:
        t_stat_vec[unconst_cols] = theta_opt[unconst_cols] / se_vec[unconst_cols]

    # Parameter draws around theta_opt.
    theta_draws = theta_opt[:, None] + C_mat @ np.random.randn(N_param, N_draws)

    filter_unc_mat = np.zeros((T, N_states, N_states), dtype=float)
    param_unc_mat = np.zeros((T, N_states, N_states), dtype=float)

    ii = 0
    used_draws = 0

    while ii < N_draws:
        theta_i = theta_draws[:, ii]

        # Skip draws that violate the constraints.
        violates = np.any(theta_i > UB) or np.any(theta_i < LB) or np.any(
            np.asarray(const_f(theta_i), dtype=float).reshape(-1) > 0
        )

        if violates:
            ii += 1
            continue

        out_i = evaluate_kalman(theta_i)

        # Smoothed covariance: uncertainty coming from the filter.
        state_smoothed_pred_i = np.asarray(_get(out_i, "state_smoothed_pred"), dtype=float)
        filter_unc_mat = filter_unc_mat + state_smoothed_pred_i

        # Difference of the smoothed states relative to the baseline.
        state_smoothed_i = np.asarray(_get(out_i, "state_smoothed"), dtype=float)
        state_smoothed_base = np.asarray(_get(out_base, "state_smoothed"), dtype=float)
        state_dif = state_smoothed_i - state_smoothed_base

        for t in range(T):
            param_unc_mat[t, :, :] = (
                param_unc_mat[t, :, :] + np.outer(state_dif[t, :], state_dif[t, :])
            )

        used_draws += 1
        ii += 1

    if used_draws == 0:
        raise RuntimeError("No accepted parameter draws in compute_se_kalman.")

    filter_unc = filter_unc_mat / used_draws
    param_unc = param_unc_mat / used_draws

    # Total state uncertainty: baseline smoothed covariance + parameter uncertainty.
    state_unc_mat = np.asarray(_get(out_base, "state_smoothed_pred"), dtype=float) + param_unc

    state_ses = np.full((T, N_states), np.nan, dtype=float)
    state_ses_no_param = np.full((T, N_states), np.nan, dtype=float)
    state_ses_contemp = np.full((T, N_states), np.nan, dtype=float)

    out_base_state_contemp_pred = np.asarray(_get(out_base, "state_contemporaneous_pred"), dtype=float)

    for t in range(T):
        state_ses[t, :] = np.sqrt(np.diag(state_unc_mat[t, :, :]))
        state_ses_no_param[t, :] = np.sqrt(np.diag(filter_unc[t, :, :]))
        state_ses_contemp[t, :] = np.sqrt(np.diag(out_base_state_contemp_pred[t, :, :]))

    return {
        "V_mat": V_mat,
        "se_vec": se_vec,
        "t_stat_vec": t_stat_vec,
        "filter_unc": filter_unc,
        "param_unc": param_unc,
        "state_unc_mat": state_unc_mat,
        "state_ses": state_ses,
        "state_ses_no_param": state_ses_no_param,
        "state_ses_contemp": state_ses_contemp,
        "used_draws": np.array([used_draws], dtype=int),
        "unconst_cols": unconst_cols,
        "const_cols": const_cols,
    }
