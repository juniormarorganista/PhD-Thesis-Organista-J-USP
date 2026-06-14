"""Finite-difference derivative operators.

This module starts with robust second-order operators.  The compact high-order
operators from the original Fortran implementation should be ported here in a
later phase and validated against manufactured solutions.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def _as_float_array(values: FloatArray) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("Expected a two-dimensional array with shape (ny, nx).")
    return array


def derivative_x(values: FloatArray, dx: float) -> FloatArray:
    """Return the numerical derivative with respect to x.

    Axis convention: x corresponds to axis 1.
    """
    if dx <= 0:
        raise ValueError("dx must be positive.")
    f = _as_float_array(values)
    out = np.empty_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dx)
    out[:, 0] = (-3.0 * f[:, 0] + 4.0 * f[:, 1] - f[:, 2]) / (2.0 * dx)
    out[:, -1] = (3.0 * f[:, -1] - 4.0 * f[:, -2] + f[:, -3]) / (2.0 * dx)
    return out


def derivative_y(values: FloatArray, dy: float) -> FloatArray:
    """Return the numerical derivative with respect to y.

    Axis convention: y corresponds to axis 0.
    """
    if dy <= 0:
        raise ValueError("dy must be positive.")
    f = _as_float_array(values)
    out = np.empty_like(f)
    out[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * dy)
    out[0, :] = (-3.0 * f[0, :] + 4.0 * f[1, :] - f[2, :]) / (2.0 * dy)
    out[-1, :] = (3.0 * f[-1, :] - 4.0 * f[-2, :] + f[-3, :]) / (2.0 * dy)
    return out


def second_derivative_x(values: FloatArray, dx: float) -> FloatArray:
    """Return the second-order finite-difference second derivative in x.

    Axis convention: x corresponds to axis 1. Boundary columns (which are not
    needed by the interior-point update in
    :mod:`mms_viscoelastic.numerics.time_stepping`) are set to zero.
    """
    if dx <= 0:
        raise ValueError("dx must be positive.")
    f = _as_float_array(values)
    out = np.zeros_like(f)
    out[:, 1:-1] = (f[:, 2:] - 2.0 * f[:, 1:-1] + f[:, :-2]) / dx**2
    return out


def second_derivative_y(values: FloatArray, dy: float) -> FloatArray:
    """Return the second-order finite-difference second derivative in y.

    Axis convention: y corresponds to axis 0. Boundary rows (which are not
    needed by the interior-point update in
    :mod:`mms_viscoelastic.numerics.time_stepping`) are set to zero.
    """
    if dy <= 0:
        raise ValueError("dy must be positive.")
    f = _as_float_array(values)
    out = np.zeros_like(f)
    out[1:-1, :] = (f[2:, :] - 2.0 * f[1:-1, :] + f[:-2, :]) / dy**2
    return out


def laplacian(values: FloatArray, dx: float, dy: float) -> FloatArray:
    """Return the second-order finite-difference Laplacian."""
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be positive.")
    f = _as_float_array(values)
    out = np.zeros_like(f)
    out[1:-1, 1:-1] = (
        (f[1:-1, 2:] - 2.0 * f[1:-1, 1:-1] + f[1:-1, :-2]) / dx**2
        + (f[2:, 1:-1] - 2.0 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / dy**2
    )
    return out
