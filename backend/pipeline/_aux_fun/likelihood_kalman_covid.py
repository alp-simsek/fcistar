from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np

from _aux_fun.kalman_filter_covid import kalman_filter_covid


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def _set(obj: Any, name: str, value: Any) -> Any:
    if isinstance(obj, dict):
        obj[name] = value
    else:
        setattr(obj, name, value)
    return obj


def likelihood_kalman_covid(
    y: np.ndarray,
    x: np.ndarray,
    theta: np.ndarray,
    map_theta_to_mats,
    settings,
) -> float:
    """
    Compute the log-likelihood of the state-space model.

    This function:
    - maps theta to the state-space matrices
    - detrends the data if needed
    - runs the Kalman filter
    - returns the sum of loglik_vec
    """
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

    # Run the Kalman filter with smoother turned off.
    settings_local = deepcopy(settings)
    _set(settings_local, "smoother", 0)

    output = kalman_filter_covid(
        y=y_work,
        x=x_work,
        F=F,
        A=A,
        H=H,
        Q=Q,
        R_base=R,
        kappa_seq=kappa_seq,
        settings=settings_local,
    )

    loglik_vec = _get(output, "loglik_vec")
    return float(np.sum(loglik_vec))