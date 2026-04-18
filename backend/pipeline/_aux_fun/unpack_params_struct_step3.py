from __future__ import annotations

import numpy as np


def unpack_params_struct_step3(theta: np.ndarray) -> dict:
    """
    Unpack the step-3 parameter vector.
    """
    theta = np.asarray(theta, dtype=float).reshape(-1)

    param = {}

    # Unpack parameters
    # MATLAB theta(1) -> Python theta[0]
    param["a_y1"] = theta[0]
    param["a_y2"] = theta[1]

    param["b_pi"] = theta[2]
    param["b_y"] = theta[3]
    param["sigma_yt"] = theta[4]
    param["sigma_pi"] = theta[5]
    param["sigma_ys"] = theta[6]
    param["sigma_pi_e"] = theta[7]
    param["sigma_delta"] = theta[8]
    param["phi"] = theta[9]

    # MATLAB: [1; theta(11:13)]
    param["kappa_vec"] = np.concatenate(([1.0], theta[10:13]))

    param["rho_pi"] = theta[13]
    param["rho_delta"] = theta[14]
    param["eta"] = theta[15]

    # MATLAB:
    # a_f = -theta(17)*(1-theta(16))/(1-theta(16)^4)
    param["a_f"] = -theta[16] * (1.0 - theta[15]) / (1.0 - theta[15] ** 4)

    param["alpha"] = 1.0
    param["c_g"] = 1.0

    return param
