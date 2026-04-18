from __future__ import annotations

from typing import Any, Callable

import numpy as np


def _get(obj: Any, name: str) -> Any:
    if isinstance(obj, dict):
        return obj[name]
    return getattr(obj, name)


def compute_opg(
    theta_opt: np.ndarray,
    constraints,
    evaluate_kalman: Callable[[np.ndarray], Any],
) -> np.ndarray:
    """
    Compute the OPG covariance matrix for the unconstrained parameters.

    This follows the MATLAB code:
    - detect constrained parameters from theta_opt == LB or theta_opt == UB
    - compute numerical derivatives using a forward difference
    - build the information matrix from the outer product of gradients
    - return inv(I_mat_free) / T
    """

    theta_opt = np.asarray(theta_opt, dtype=float).reshape(-1)

    UB = np.asarray(_get(constraints, "UB"), dtype=float).reshape(-1)
    LB = np.asarray(_get(constraints, "LB"), dtype=float).reshape(-1)

    cons_vec = (theta_opt == LB).astype(int) + (theta_opt == UB).astype(int)

    unconst_cols = np.where(cons_vec == 0)[0]
    const_cols = np.where(cons_vec == 1)[0]

    delta = 1e-5  # step size for numerical derivatives

    out_base = evaluate_kalman(theta_opt)
    lik_vec_base = np.asarray(_get(out_base, "loglik_vec"), dtype=float).reshape(-1)

    T = len(lik_vec_base)
    N_param = len(theta_opt)

    lik_gradient = np.full((T, N_param), np.nan, dtype=float)

    for i in range(N_param):
        if cons_vec[i] == 0:
            theta_prime = theta_opt.copy()
            theta_prime[i] = theta_opt[i] + delta

            output_prime = evaluate_kalman(theta_prime)
            lik_vec_prime = np.asarray(_get(output_prime, "loglik_vec"), dtype=float).reshape(-1)

            lik_gradient[:, i] = (lik_vec_prime - lik_vec_base) / delta

    # Information matrix from the outer product of gradients.
    I_mat = (lik_gradient.T @ lik_gradient) / T
    I_mat_aux = I_mat[np.ix_(unconst_cols, unconst_cols)]

    V_mat_aux = np.linalg.inv(I_mat_aux) / T
    return V_mat_aux