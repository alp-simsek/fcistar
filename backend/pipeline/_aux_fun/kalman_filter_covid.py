from __future__ import annotations

from typing import Any

import numpy as np
from scipy.linalg import solve_discrete_lyapunov


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def solve_ricatti(F: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Fallback used when P_0 is not supplied."""
    return solve_discrete_lyapunov(F, Q)


def kalman_filter_covid(
    y: np.ndarray,
    x: np.ndarray,
    F: np.ndarray,
    A: np.ndarray,
    H: np.ndarray,
    Q: np.ndarray,
    R_base: np.ndarray,
    kappa_seq: np.ndarray,
    settings,
) -> dict[str, np.ndarray]:
    """
    Run the Kalman filter.

    State equation
        xi_{t+1} = F xi_t + v_{t+1}

    Measurement equation
        y_t = A' x_t + H' xi_t + w_t
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    F = np.asarray(F, dtype=float)
    A = np.asarray(A, dtype=float)
    H = np.asarray(H, dtype=float)
    Q = np.asarray(Q, dtype=float)
    R_base = np.asarray(R_base, dtype=float)
    kappa_seq = np.asarray(kappa_seq, dtype=float).reshape(-1)

    T_y, n_y = y.shape
    T_x, _ = x.shape
    n_s = F.shape[0]

    if T_y != T_x:
        raise ValueError("y and x have different number of observations")

    T = T_y

    # xi_{t|t}, contemporaneous prediction
    xi_mat_c = np.zeros((T, n_s))
    # xi_{t+1|t}, one-step-ahead prediction
    xi_mat = np.zeros((T, n_s))
    # xi_{t|T}, smoothed
    xi_smoothed_mat = np.zeros((T, n_s))

    # P_{t|t}
    P_mat_c = np.zeros((T, n_s, n_s))
    # P_{t+1|t}
    P_mat = np.zeros((T, n_s, n_s))
    # P_{t|T}, smoothed
    P_smoothed_mat = np.zeros((T, n_s, n_s))

    # J_t
    J_mat = np.zeros((T, n_s, n_s))
    loglik_vec = np.zeros(T)

    # Initial value of the state
    try:
        xi_0 = np.asarray(_get(settings, "xi_0"), dtype=float).reshape(-1)
    except Exception:
        xi_0 = np.zeros(n_s)

    # Initial value of the covariance
    try:
        P_0 = np.asarray(_get(settings, "P_0"), dtype=float)
    except Exception:
        P_0 = solve_ricatti(F, Q)

    for t in range(T):
        if t == 0:
            P_tt = P_0
            xi_tt = xi_0
        else:
            P_tt = P_mat_c[t - 1, :, :]
            xi_tt = xi_mat_c[t - 1, :]

        # Time-varying measurement covariance
        R = R_base * (kappa_seq[t] ** 2)

        # One-step-ahead prediction
        xi_pred = F @ xi_tt
        xi_mat[t, :] = xi_pred

        P_ttm1 = F @ P_tt @ F.T + Q
        P_mat[t, :, :] = P_ttm1

        # Prediction error
        pred_error = y[t, :] - (A.T @ x[t, :]) - (H.T @ xi_pred)
        HPHR = H.T @ P_ttm1 @ H + R

        # Likelihood contribution of observation t
        sign, logdet = np.linalg.slogdet(HPHR)
        if sign <= 0:
            raise np.linalg.LinAlgError("Measurement covariance HPHR is not positive definite.")

        HPHR_inv_pred_error = np.linalg.solve(HPHR, pred_error)
        loglik_vec[t] = (
            -(n_y / 2.0) * np.log(2.0 * np.pi)
            - 0.5 * logdet
            - 0.5 * pred_error.T @ HPHR_inv_pred_error
        )

        # Contemporaneous prediction
        xi_upd = xi_pred + (P_ttm1 @ H) @ HPHR_inv_pred_error
        xi_mat_c[t, :] = xi_upd

        # Contemporaneous prediction error
        P_upd = P_ttm1 - P_ttm1 @ H @ np.linalg.solve(HPHR, H.T @ P_ttm1)
        P_mat_c[t, :, :] = P_upd

    # Run the smoother if requested
    if int(_get(settings, "smoother")) == 1:
        states_smoother = np.asarray(_get(settings, "states_smoother"), dtype=int)

        # MATLAB indices are 1-based; convert when passed through directly.
        if states_smoother.min() >= 1:
            states_smoother = states_smoother - 1

        xi_smoothed_mat[T - 1, states_smoother] = xi_mat_c[T - 1, states_smoother]
        P_smoothed_mat[np.ix_([T - 1], states_smoother, states_smoother)] = (
            P_mat_c[T - 1][np.ix_(states_smoother, states_smoother)]
        )

        for t in range(T - 2, -1, -1):
            P_tt = P_mat_c[t][np.ix_(states_smoother, states_smoother)]
            P_tp1t = P_mat[t + 1][np.ix_(states_smoother, states_smoother)]

            J_t = P_tt @ F[np.ix_(states_smoother, states_smoother)].T @ np.linalg.pinv(P_tp1t)
            J_mat[t][np.ix_(states_smoother, states_smoother)] = J_t

            xi_tt = xi_mat_c[t, states_smoother]
            xi_tp1_pred = xi_mat[t + 1, states_smoother]
            xi_tp1_smooth = xi_smoothed_mat[t + 1, states_smoother]

            xi_tT = xi_tt + J_t @ (xi_tp1_smooth - xi_tp1_pred)
            P_tT = (
                P_tt
                + J_t
                @ (
                    P_smoothed_mat[t + 1][np.ix_(states_smoother, states_smoother)]
                    - P_tp1t
                )
                @ J_t.T
            )

            xi_smoothed_mat[t, states_smoother] = xi_tT
            P_smoothed_mat[t][np.ix_(states_smoother, states_smoother)] = P_tT

    return {
        "state_contemporaneous": xi_mat_c,
        "state_contemporaneous_pred": P_mat_c,
        "state_one_ahead": xi_mat,
        "state_one_ahead_pred": P_mat,
        "state_smoothed": xi_smoothed_mat,
        "state_smoothed_pred": P_smoothed_mat,
        "J_mat": J_mat,
        "loglik_vec": loglik_vec,
    }