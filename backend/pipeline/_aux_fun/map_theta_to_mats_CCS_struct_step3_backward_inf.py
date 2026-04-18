from __future__ import annotations

import numpy as np

from _aux_fun.unpack_params_struct_step3 import unpack_params_struct_step3


def map_theta_to_mats_CCS_struct_step3_backward_inf(
    theta: np.ndarray,
    fixed_param,
    covid_dummies: np.ndarray,
):
    """
    Map theta to the state-space matrices used in step 3.

    State order:
    xi_t = [
        y_t^*, y_{t-1}^*, y_{t-2}^*,
        g_t, g_{t-1}, g_{t-2},
        FCI_t^*, FCI_{t-1}^*,
        delta_t, delta_{t-1}, delta_{t-2}
    ]
    """
    theta = np.asarray(theta, dtype=float)
    covid_dummies = np.asarray(covid_dummies, dtype=float)

    # Unpack parameters
    param = unpack_params_struct_step3(theta)

    a_y1 = param["a_y1"]
    a_y2 = param["a_y2"]
    a_y3 = a_y2
    a_f = param["a_f"]
    b_pi = param["b_pi"]
    b_y = param["b_y"]
    sigma_yt = param["sigma_yt"]
    sigma_ys = param["sigma_ys"]
    sigma_pi = param["sigma_pi"]
    sigma_pi_e = param["sigma_pi_e"]
    sigma_delta = param["sigma_delta"]
    phi = param["phi"]
    kappa_vec = np.asarray(param["kappa_vec"], dtype=float).reshape(-1)
    rho_delta = param["rho_delta"]
    rho_pi = param["rho_pi"]
    eta = param["eta"]
    c_g = 1.0

    # MATLAB fixed_param(1) -> Python fixed_param[0]
    lambda_g = float(np.asarray(fixed_param).reshape(-1)[0])
    sigma_g = lambda_g * sigma_ys

    # Transition matrix
    F = np.zeros((11, 11), dtype=float)

    # y* equation and lags
    F[0, 0] = 1.0
    F[0, 3] = 1.0
    F[1, 0] = 1.0
    F[2, 1] = 1.0

    # g equation and lags
    F[3, 3] = 1.0
    F[4, 3] = 1.0
    F[5, 4] = 1.0

    # delta equation and lags
    F[8, 8] = rho_delta
    F[9, 8] = 1.0
    F[10, 9] = 1.0

    # FCI* equation
    # In the code a_f is defined with a minus, so the signs differ from the paper.
    F[6, 0] = -eta / a_f
    F[6, 1] = eta / a_f
    F[6, 3] = c_g * (1.0 - eta) / a_f
    F[6, 4] = c_g * eta / a_f
    F[6, 6] = eta
    F[6, 8] = (1.0 / a_f) * (eta + (1.0 - rho_delta) * rho_delta)
    F[6, 9] = -(1.0 / a_f) * eta * rho_delta
    F[7, 6] = 1.0

    # State innovation covariance
    Q = np.zeros((11, 11), dtype=float)

    Q[0, 0] = sigma_ys**2
    Q[3, 3] = sigma_g**2
    Q[8, 8] = sigma_delta**2

    # Covariances of FCI* innovations
    Q[6, 0] = (1.0 / a_f) * sigma_ys**2
    Q[0, 6] = Q[6, 0]

    Q[6, 3] = c_g * (1.0 / a_f) * sigma_g**2
    Q[3, 6] = Q[6, 3]

    Q[6, 8] = -(1.0 / a_f) * (eta + rho_delta) * sigma_delta**2
    Q[8, 6] = Q[6, 8]

    Q[6, 6] = (1.0 / a_f) ** 2 * (
        (c_g**2) * sigma_g**2
        + sigma_ys**2
        + ((eta + rho_delta) ** 2) * sigma_delta**2
    )

    # Measurement matrix
    H_t = np.zeros((2, 11), dtype=float)

    # Output relation
    H_t[0, 1] = 1.0
    H_t[0, 2] = -1.0
    H_t[0, 4] = 1.0
    H_t[0, 5] = -1.0
    H_t[0, 6] = -a_f
    H_t[0, 8] = 1.0
    H_t[0, 9] = -(1.0 + rho_delta)
    H_t[0, 10] = rho_delta

    # Phillips curve
    H_t[1, 2] = -b_y
    H_t[1, 5] = -b_y
    H_t[1, 9] = -b_y
    H_t[1, 10] = b_y * rho_delta

    H = H_t.T

    # Exogenous regressor loadings
    # x_t = [y_{t-1}, FCI_{t-1}^m, pi_{t-1}, pi_{t-2:4}, d_t, d_{t-1}]
    A_t = np.zeros((2, 6), dtype=float)

    A_t[0, 0] = 1.0
    A_t[0, 1] = a_f
    A_t[0, 4] = phi
    A_t[0, 5] = -phi

    A_t[1, 0] = b_y
    A_t[1, 2] = b_pi
    A_t[1, 3] = 1.0 - b_pi
    A_t[1, 5] = -b_y * phi

    A = A_t.T

    # Measurement error covariance
    R = np.array(
        [
            [sigma_yt**2, 0.0],
            [0.0, sigma_pi**2],
        ],
        dtype=float,
    )

    # Time-varying scaling of the measurement covariance
    kappa_seq = covid_dummies @ kappa_vec

    return {
        "F": F,
        "A": A,
        "H": H,
        "Q": Q,
        "R": R,
        "kappa_seq": kappa_seq,
    }