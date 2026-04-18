from __future__ import annotations

from typing import Any

import numpy as np

from _aux_fun.kalman_filter_covid import kalman_filter_covid


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def evaluate_kalman_covid(
    y: np.ndarray,
    x: np.ndarray,
    theta: np.ndarray,
    map_theta_to_mats,
    settings,
):
    """
    Evaluate the Kalman filter for a given parameter vector theta.

    This function:
    - maps theta to the state-space matrices
    - detrends the data if needed
    - runs the Kalman filter
    - retrends the state estimates
    """

    # Work on copies to avoid mutating the caller inputs.
    y_work = np.array(y, dtype=float, copy=True)
    x_work = np.array(x, dtype=float, copy=True)
    theta = np.asarray(theta, dtype=float)

    fixed_param = _get(settings, "fixed_param")
    covid_dummies = _get(settings, "covid_dummies")
    matrices = map_theta_to_mats(theta, fixed_param, covid_dummies)

    detrend_y = _get(settings, "detrend_y")
    T = y_work.shape[0]

    # MATLAB theta(5) -> Python theta[4]
    g = theta[4]

    if detrend_y == 1:
        y_work[:, 0] = y_work[:, 0] - g * np.arange(1, T + 1)
        x_work[:, 0] = x_work[:, 0] - g * np.arange(0, T)
        x_work[:, 1] = x_work[:, 1] - g * np.arange(-1, T - 1)

    elif detrend_y == 2:
        y_work[:, 0] = y_work[:, 0] - g * np.arange(1, T + 1)
        x_work[:, 0] = x_work[:, 0] - g * np.arange(0, T)
        x_work[:, 1] = x_work[:, 1] - g * np.arange(-1, T - 1)
        x_work[:, 2] = x_work[:, 2] - g * np.arange(-2, T - 2)

    F = _get(matrices, "F")
    A = _get(matrices, "A")
    H = _get(matrices, "H")
    Q = _get(matrices, "Q")
    R = _get(matrices, "R")
    kappa_seq = _get(matrices, "kappa_seq")

    # Run the Kalman filter.
    output = kalman_filter_covid(
        y=y_work,
        x=x_work,
        F=F,
        A=A,
        H=H,
        Q=Q,
        R_base=R,
        kappa_seq=kappa_seq,
        settings=settings,
    )

    # Retrend the state estimates.
    if detrend_y == 1:
        output["state_contemporaneous"][:, 0] += g * np.arange(1, T + 1)
        output["state_contemporaneous"][:, 1] += g * np.arange(0, T)
        output["state_contemporaneous"][:, 2] += g * np.arange(-1, T - 1)

        output["state_one_ahead"][:, 0] += g * np.arange(2, T + 2)
        output["state_one_ahead"][:, 1] += g * np.arange(1, T + 1)
        output["state_one_ahead"][:, 2] += g * np.arange(0, T)

        output["state_smoothed"][:, 0] += g * np.arange(1, T + 1)
        output["state_smoothed"][:, 1] += g * np.arange(0, T)
        output["state_smoothed"][:, 2] += g * np.arange(-1, T - 1)

    elif detrend_y == 2:
        output["state_contemporaneous"][:, 0] += g * np.arange(1, T + 1)
        output["state_contemporaneous"][:, 1] += g * np.arange(0, T)
        output["state_contemporaneous"][:, 2] += g * np.arange(-1, T - 1)
        output["state_contemporaneous"][:, 3] += g * np.arange(-2, T - 2)

        output["state_one_ahead"][:, 0] += g * np.arange(2, T + 2)
        output["state_one_ahead"][:, 1] += g * np.arange(1, T + 1)
        output["state_one_ahead"][:, 2] += g * np.arange(0, T)
        output["state_one_ahead"][:, 3] += g * np.arange(-1, T - 1)

        output["state_smoothed"][:, 0] += g * np.arange(1, T + 1)
        output["state_smoothed"][:, 1] += g * np.arange(0, T)
        output["state_smoothed"][:, 2] += g * np.arange(-1, T - 1)
        output["state_smoothed"][:, 3] += g * np.arange(-2, T - 2)

    return output