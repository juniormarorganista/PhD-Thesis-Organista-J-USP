"""Lele-type compact spatial filter.

Ports the filter applied to vorticity and the polymeric stress components
after each time step, from
``src/fortran_original/4Order_LPOG_Euler/loop.f`` (``rhs_tridfx_wz``,
``rhs_tridfx_T``, ``lhs_tridfx``, and their y-direction counterparts).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .derivatives import _as_float_array, solve_tridiagonal

FloatArray = NDArray[np.float64]

#: Default filter parameter, matching ``alphaf`` in
#: ``src/fortran_original/4Order_LPOG_Euler/par.nn``. ``alphaf = 0.5``
#: makes the filter the identity (no filtering).
DEFAULT_ALPHA_F = 0.48


def filter_coefficients(alpha_f: float = DEFAULT_ALPHA_F) -> tuple[float, float, float, float]:
    """Return the explicit filter weights ``(af, bf, cf, df)`` for ``alpha_f``.

    Matches the ``af``/``bf``/``cf``/``df`` parameters in
    ``src/fortran_original/4Order_LPOG_Euler/par.nn``.
    """
    af = (11.0 + 10.0 * alpha_f) / 16.0
    bf = (15.0 + 34.0 * alpha_f) / 64.0
    cf = (-3.0 + 6.0 * alpha_f) / 32.0
    df = (1.0 - 2.0 * alpha_f) / 64.0
    return af, bf, cf, df


def _filter_lhs_coefficients(n: int, alpha_f: float) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Tridiagonal LHS coefficients for the filter (``lhs_tridfx``/``lhs_tridfy``).

    Identity at the three outermost points on each side; ``(alpha_f, 1,
    alpha_f)`` at the interior points (0-indexed ``i = 3, ..., n - 4``).
    """
    if n < 7:
        raise ValueError("The Lele filter requires at least 7 points along the filtered axis.")
    a = np.zeros(n, dtype=np.float64)
    b = np.ones(n, dtype=np.float64)
    c = np.zeros(n, dtype=np.float64)
    a[3:-3] = alpha_f
    c[3:-3] = alpha_f
    return a, b, c


def _filter_rhs(f: FloatArray, alpha_f: float) -> FloatArray:
    """RHS of the filter along axis 0 of ``f`` (``rhs_tridfx_wz``/``rhs_tridfx_T``).

    The same formula applies to vorticity and to each stress component in
    the original Fortran (``rhs_tridfx_wz`` and ``rhs_tridfx_T`` share the
    identical weights); this helper is therefore generic over the filtered
    field.
    """
    n = f.shape[0]
    if n < 7:
        raise ValueError("The Lele filter requires at least 7 points along the filtered axis.")
    af, bf, cf, df = filter_coefficients(alpha_f)
    rhs = np.empty_like(f)
    rhs[0] = (15.0 * f[0] + 4.0 * f[1] - 6.0 * f[2] + 4.0 * f[3] - f[4]) / 16.0
    rhs[1] = (f[0] + 12.0 * f[1] + 6.0 * f[2] - 4.0 * f[3] + f[4]) / 16.0
    rhs[2] = (-f[0] + 4.0 * f[1] + 10.0 * f[2] + 4.0 * f[3] - f[4]) / 16.0
    rhs[3:-3] = (
        af * f[3:-3]
        + bf * (f[4:-2] + f[2:-4])
        + cf * (f[5:-1] + f[1:-5])
        + df * (f[6:] + f[:-6])
    )
    rhs[-3] = (-f[-1] + 4.0 * f[-2] + 10.0 * f[-3] + 4.0 * f[-4] - f[-5]) / 16.0
    rhs[-2] = (f[-1] + 12.0 * f[-2] + 6.0 * f[-3] - 4.0 * f[-4] + f[-5]) / 16.0
    rhs[-1] = (15.0 * f[-1] + 4.0 * f[-2] - 6.0 * f[-3] + 4.0 * f[-4] - f[-5]) / 16.0
    return rhs


def lele_filter_x(values: FloatArray, alpha_f: float = DEFAULT_ALPHA_F) -> FloatArray:
    """Apply the compact filter along x.

    Ports ``rhs_tridfx_wz``/``rhs_tridfx_T`` + ``lhs_tridfx`` + ``tridx``
    from ``src/fortran_original/4Order_LPOG_Euler/loop.f``. Axis
    convention: x corresponds to axis 1, matching
    :func:`mms_viscoelastic.numerics.derivatives.compact_derivative_x`.
    """
    f = _as_float_array(values)
    f_t = f.T
    a, b, c = _filter_lhs_coefficients(f_t.shape[0], alpha_f)
    rhs = _filter_rhs(f_t, alpha_f)
    return solve_tridiagonal(a, b, c, rhs).T


def lele_filter_y(values: FloatArray, alpha_f: float = DEFAULT_ALPHA_F) -> FloatArray:
    """Apply the compact filter along y.

    Ports ``rhs_tridfy_wz``/``rhs_tridfy_T`` + ``lhs_tridfy`` + ``tridy``
    from ``src/fortran_original/4Order_LPOG_Euler/loop.f``. Axis
    convention: y corresponds to axis 0, matching
    :func:`mms_viscoelastic.numerics.derivatives.compact_derivative_y`.
    """
    f = _as_float_array(values)
    a, b, c = _filter_lhs_coefficients(f.shape[0], alpha_f)
    rhs = _filter_rhs(f, alpha_f)
    return solve_tridiagonal(a, b, c, rhs)


__all__ = [
    "DEFAULT_ALPHA_F",
    "filter_coefficients",
    "lele_filter_x",
    "lele_filter_y",
]
