"""Finite-difference derivative operators.

This module provides the robust second-order operators used by the early
phases of the port, and the fourth-order compact (Padé-type) operators that
reproduce the scheme actually used in the thesis, ported from
``src/fortran_original/4Order_LPOG_Euler/derivs.f`` (``coef_4a``, ``coeff_4a``,
``rhsx_4a``, ``rhsxx_4a``, ``rhsy_4a``, ``rhsyy_4a``) and the Thomas-algorithm
solver in ``solvers.f`` (``tridx``, ``tridy``).
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


def solve_tridiagonal(a: FloatArray, b: FloatArray, c: FloatArray, rhs: FloatArray) -> FloatArray:
    """Solve ``A x = rhs`` for a tridiagonal ``A`` via the Thomas algorithm (TDMA).

    Mirrors ``tridx``/``tridy`` in
    ``src/fortran_original/4Order_LPOG_Euler/solvers.f``. ``a``, ``b`` and
    ``c`` are the sub-, main- and super-diagonal of ``A`` (each of length
    ``n``, with ``a[0]`` and ``c[-1]`` unused). ``rhs`` has shape ``(n, m)``:
    each of the ``m`` columns is an independent right-hand side sharing the
    same ``A`` (the same way the Fortran routine loops over ``j``/``i`` while
    reusing ``a, b, c``).
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64)
    n = b.shape[0]
    if a.shape[0] != n or c.shape[0] != n:
        raise ValueError("a, b and c must have the same length.")
    rhs = np.asarray(rhs, dtype=np.float64)
    if rhs.ndim != 2 or rhs.shape[0] != n:
        raise ValueError("rhs must have shape (n, m) with n matching a/b/c.")

    gam = np.empty(n, dtype=np.float64)
    u = np.empty_like(rhs)

    bet = b[0]
    u[0] = rhs[0] / bet
    for i in range(1, n):
        gam[i] = c[i - 1] / bet
        bet = b[i] - a[i] * gam[i]
        u[i] = (rhs[i] - a[i] * u[i - 1]) / bet

    for i in range(n - 2, -1, -1):
        u[i] = u[i] - gam[i + 1] * u[i + 1]

    return u


def _compact_first_derivative_coefficients(n: int) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Tridiagonal LHS coefficients for the compact 1st derivative (``coef_4a``)."""
    a = np.ones(n, dtype=np.float64)
    b = np.full(n, 4.0, dtype=np.float64)
    c = np.ones(n, dtype=np.float64)
    a[0], b[0], c[0] = 0.0, 1.0, 3.0
    a[-1], b[-1], c[-1] = 3.0, 1.0, 0.0
    return a, b, c


def _compact_second_derivative_coefficients(n: int) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Tridiagonal LHS coefficients for the compact 2nd derivative (``coeff_4a``)."""
    a = np.ones(n, dtype=np.float64)
    b = np.full(n, 10.0, dtype=np.float64)
    c = np.ones(n, dtype=np.float64)
    a[0], b[0], c[0] = 0.0, 1.0, 11.0
    a[-1], b[-1], c[-1] = 11.0, 1.0, 0.0
    return a, b, c


def _compact_first_derivative_rhs(f: FloatArray, h: float) -> FloatArray:
    """RHS for the compact 1st derivative along axis 0 of ``f`` (``rhsx_4a``/``rhsy_4a``)."""
    n = f.shape[0]
    if n < 4:
        raise ValueError("The compact scheme requires at least 4 points along the derivative axis.")
    rhs = np.empty_like(f)
    rhs[0] = (-17.0 * f[0] + 9.0 * f[1] + 9.0 * f[2] - f[3]) / (6.0 * h)
    rhs[1:-1] = 3.0 * (f[2:] - f[:-2]) / h
    rhs[-1] = (-17.0 * f[-1] + 9.0 * f[-2] + 9.0 * f[-3] - f[-4]) / (-6.0 * h)
    return rhs


def _compact_second_derivative_rhs(f: FloatArray, h: float) -> FloatArray:
    """RHS for the compact 2nd derivative along axis 0 of ``f`` (``rhsxx_4a``/``rhsyy_4a``)."""
    n = f.shape[0]
    if n < 4:
        raise ValueError("The compact scheme requires at least 4 points along the derivative axis.")
    hh = h * h
    rhs = np.empty_like(f)
    rhs[0] = (39.0 * f[0] - 81.0 * f[1] + 45.0 * f[2] - 3.0 * f[3]) / (3.0 * hh)
    rhs[1:-1] = (12.0 * (f[2:] + f[:-2]) - 24.0 * f[1:-1]) / hh
    rhs[-1] = (39.0 * f[-1] - 81.0 * f[-2] + 45.0 * f[-3] - 3.0 * f[-4]) / (3.0 * hh)
    return rhs


def compact_derivative_x(values: FloatArray, dx: float) -> FloatArray:
    """Return the fourth-order compact derivative with respect to x.

    Ports ``derx`` (``coef_4a`` + ``rhsx_4a`` + ``tridx``) from
    ``src/fortran_original/4Order_LPOG_Euler/derivs.f``. Axis convention: x
    corresponds to axis 1, matching :func:`derivative_x`.
    """
    if dx <= 0:
        raise ValueError("dx must be positive.")
    f = _as_float_array(values)
    f_t = f.T
    a, b, c = _compact_first_derivative_coefficients(f_t.shape[0])
    rhs = _compact_first_derivative_rhs(f_t, dx)
    return solve_tridiagonal(a, b, c, rhs).T


def compact_derivative_y(values: FloatArray, dy: float) -> FloatArray:
    """Return the fourth-order compact derivative with respect to y.

    Ports ``dery`` (``coef_4a`` + ``rhsy_4a`` + ``tridy``) from
    ``src/fortran_original/4Order_LPOG_Euler/derivs.f``. Axis convention: y
    corresponds to axis 0, matching :func:`derivative_y`.
    """
    if dy <= 0:
        raise ValueError("dy must be positive.")
    f = _as_float_array(values)
    a, b, c = _compact_first_derivative_coefficients(f.shape[0])
    rhs = _compact_first_derivative_rhs(f, dy)
    return solve_tridiagonal(a, b, c, rhs)


def compact_second_derivative_x(values: FloatArray, dx: float) -> FloatArray:
    """Return the fourth-order compact second derivative with respect to x.

    Ports ``derxx`` (``coeff_4a`` + ``rhsxx_4a`` + ``tridx``) from
    ``src/fortran_original/4Order_LPOG_Euler/derivs.f``.

    The interior stencil (``f_{i-1}''+10f_i''+f_{i+1}''=[12(f_{i+1}+f_{i-1})
    -24f_i]/h^2``) is fourth-order accurate. The boundary closure
    (``f_1''+11f_2''=[39f_1-81f_2+45f_3-3f_4]/(3h^2)``) is only third-order
    accurate: substituting the exact Taylor expansion gives a residual of
    ``(1/12) h^3 f^(5)``, one order below the interior (verified
    symbolically; the analogous first-derivative closure in
    :func:`compact_derivative_x` has residual ``(1/20) h^4 f^(5)``, i.e. it
    matches the interior order). This is a property of the boundary closure
    itself, not a porting defect: both points and the formula above are
    transcribed from ``coeff_4a``/``rhsxx_4a`` unchanged. It does not
    degrade the convergence order reported for Case 12, because the
    boundary rows of the marched fields are subsequently overwritten by the
    exact manufactured boundary condition at every time step, never read
    from this operator's own boundary output.
    """
    if dx <= 0:
        raise ValueError("dx must be positive.")
    f = _as_float_array(values)
    f_t = f.T
    a, b, c = _compact_second_derivative_coefficients(f_t.shape[0])
    rhs = _compact_second_derivative_rhs(f_t, dx)
    return solve_tridiagonal(a, b, c, rhs).T


def compact_second_derivative_y(values: FloatArray, dy: float) -> FloatArray:
    """Return the fourth-order compact second derivative with respect to y.

    Ports ``deryy`` (``coeff_4a`` + ``rhsyy_4a`` + ``tridy``) from
    ``src/fortran_original/4Order_LPOG_Euler/derivs.f``. See
    :func:`compact_second_derivative_x` for the boundary-closure truncation
    order (third-order, not fourth-order, at the two points nearest each
    boundary).
    """
    if dy <= 0:
        raise ValueError("dy must be positive.")
    f = _as_float_array(values)
    a, b, c = _compact_second_derivative_coefficients(f.shape[0])
    rhs = _compact_second_derivative_rhs(f, dy)
    return solve_tridiagonal(a, b, c, rhs)
