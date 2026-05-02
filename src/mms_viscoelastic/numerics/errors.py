"""Error norms and convergence-rate utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def l2_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the discrete L2 error between numerical and exact arrays."""
    numerical = np.asarray(numerical, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    if numerical.shape != exact.shape:
        raise ValueError("numerical and exact arrays must have the same shape.")
    return float(np.sqrt(np.mean((numerical - exact) ** 2)))


def relative_l2_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the relative discrete L2 error."""
    denominator = l2_error(np.zeros_like(exact, dtype=np.float64), exact)
    if denominator == 0:
        raise ValueError("Cannot compute relative error because exact field has zero L2 norm.")
    return l2_error(numerical, exact) / denominator


def linf_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the maximum absolute error."""
    numerical = np.asarray(numerical, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    if numerical.shape != exact.shape:
        raise ValueError("numerical and exact arrays must have the same shape.")
    return float(np.max(np.abs(numerical - exact)))


def estimate_order(h_values: FloatArray, errors: FloatArray) -> float:
    """Estimate convergence order from mesh sizes and errors using log-log fit."""
    h = np.asarray(h_values, dtype=np.float64)
    err = np.asarray(errors, dtype=np.float64)
    if h.ndim != 1 or err.ndim != 1 or h.size != err.size:
        raise ValueError("h_values and errors must be one-dimensional arrays with same size.")
    if h.size < 2:
        raise ValueError("At least two mesh levels are required.")
    if np.any(h <= 0) or np.any(err <= 0):
        raise ValueError("Mesh sizes and errors must be positive.")
    slope, _ = np.polyfit(np.log(h), np.log(err), deg=1)
    return float(slope)
