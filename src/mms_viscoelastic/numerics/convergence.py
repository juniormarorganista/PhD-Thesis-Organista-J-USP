"""Grid-convergence verification for MMS case 12.

This module verifies that the fourth-order compact finite-difference
operators in :mod:`mms_viscoelastic.numerics.derivatives`
(``compact_derivative_x``/``compact_derivative_y``, ported from
``derivs.f``'s ``coef_4a``/``rhsx_4a`` family) reproduce the differential
relations satisfied by the exact manufactured fields of case 12
(:mod:`mms_viscoelastic.manufactured_solutions.case_12``) at an order
consistent with the thesis and the JNNFM paper (observed order
approaching, and in places exceeding, 4 for the mesh range considered; see
the pre-asymptotic discussion in the paper for why the two-point order
estimator can lie above or below 4 at finite mesh spacing).

The stream-function formulation used by the manufactured solution defines

.. math::

    u = \\frac{\\partial \\psi}{\\partial y}, \\qquad
    v = -\\frac{\\partial \\psi}{\\partial x}, \\qquad
    \\omega_z = \\frac{\\partial u}{\\partial y} - \\frac{\\partial v}{\\partial x}.

Each relation is evaluated with the discrete operators ``compact_derivative_x``
and ``compact_derivative_y`` applied to the exact fields, and the result is
compared to the corresponding exact field using the discrete :math:`L_2` and
:math:`L_\\infty` norms from :mod:`mms_viscoelastic.numerics.errors`. The
observed order of accuracy ``p`` is estimated from a log-log fit of the
error against the mesh size ``h`` across the requested grids.

This module does not exercise the time-marching pipeline
(:mod:`mms_viscoelastic.numerics.time_stepping`,
:mod:`mms_viscoelastic.numerics.time_marching`) or the Poisson solver
(:mod:`mms_viscoelastic.numerics.poisson`), which still use second-order
operators pending a separate port of the pentadiagonal Poisson discretization
and of the compact derivatives/filter into the time-marching steps.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions.case_12 import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.numerics.derivatives import compact_derivative_x, compact_derivative_y
from mms_viscoelastic.numerics.errors import estimate_order, l2_error, linf_error

FloatArray = NDArray[np.float64]

QUANTITIES: tuple[str, ...] = ("u", "v", "wz")
"""Field components checked by the convergence study, in table order."""


@dataclass(frozen=True)
class ConvergenceLevel:
    """Errors of a single grid in the convergence study.

    Attributes
    ----------
    n:
        Number of cells along each direction (the grid has ``n + 1`` points).
    h:
        Mesh size ``dx`` (equal to ``dy`` for the uniform unit-square grid).
    l2_errors:
        Discrete :math:`L_2` error of each quantity in :data:`QUANTITIES`.
    linf_errors:
        Discrete :math:`L_\\infty` error of each quantity in
        :data:`QUANTITIES`.
    """

    n: int
    h: float
    l2_errors: dict[str, float]
    linf_errors: dict[str, float]


@dataclass(frozen=True)
class ConvergenceReport:
    """Result of :func:`run_convergence_study`.

    Attributes
    ----------
    levels:
        One :class:`ConvergenceLevel` per requested grid, ordered as given in
        ``grid_sizes``.
    l2_orders:
        Estimated order of accuracy ``p`` for each quantity, computed from the
        :math:`L_2` errors of all levels via a log-log fit.
    linf_orders:
        Estimated order of accuracy ``p`` for each quantity, computed from the
        :math:`L_\\infty` errors of all levels via a log-log fit.
    """

    levels: list[ConvergenceLevel]
    l2_orders: dict[str, float]
    linf_orders: dict[str, float]

    def to_rows(self) -> list[dict[str, float | int]]:
        """Flatten the report into a list of dictionaries (one per grid).

        Each row contains ``n``, ``h``, and ``l2_<q>`` / ``linf_<q>`` for
        every quantity ``q`` in :data:`QUANTITIES`. This is convenient for
        building a :class:`pandas.DataFrame` without adding a hard
        dependency on pandas in this module.
        """
        rows: list[dict[str, float | int]] = []
        for level in self.levels:
            row: dict[str, float | int] = {"n": level.n, "h": level.h}
            for quantity in QUANTITIES:
                row[f"l2_{quantity}"] = level.l2_errors[quantity]
                row[f"linf_{quantity}"] = level.linf_errors[quantity]
            rows.append(row)
        return rows


def _stream_function_derived_fields(
    psi: FloatArray, u: FloatArray, v: FloatArray, dx: float, dy: float
) -> dict[str, FloatArray]:
    """Return numerical ``u``, ``v`` and ``wz`` derived from ``psi``/``u``/``v``."""
    u_num = compact_derivative_y(psi, dy)
    v_num = -compact_derivative_x(psi, dx)
    wz_num = compact_derivative_y(u, dy) - compact_derivative_x(v, dx)
    return {"u": u_num, "v": v_num, "wz": wz_num}


def run_convergence_study(
    grid_sizes: list[int],
    *,
    t: float = 0.0,
    parameters: Case12Parameters | None = None,
) -> ConvergenceReport:
    """Run the grid-convergence study for MMS case 12.

    Parameters
    ----------
    grid_sizes:
        Numbers of cells ``n`` along each direction. Each grid has
        ``nx = ny = n + 1`` points and mesh size ``h = 1 / n`` on the unit
        square. Must contain at least two values, in any order; the order of
        the resulting :attr:`ConvergenceReport.levels` follows ``grid_sizes``.
    t:
        Physical time at which the exact fields are evaluated.
    parameters:
        Case-12 parameters. Defaults to :class:`Case12Parameters` defaults.

    Returns
    -------
    ConvergenceReport
        Per-grid errors and estimated orders of accuracy for ``u``, ``v`` and
        ``wz`` (see module docstring for the differential relations checked).

    Raises
    ------
    ValueError
        If ``grid_sizes`` has fewer than two elements.
    """
    if len(grid_sizes) < 2:
        raise ValueError("grid_sizes must contain at least two grid resolutions.")

    levels: list[ConvergenceLevel] = []
    for n in grid_sizes:
        grid = Grid2D.from_bounds(nx=n + 1, ny=n + 1)
        field = evaluate_case_12_exact(grid, t=t, parameters=parameters)
        assert field.psi is not None
        assert field.u is not None
        assert field.v is not None
        derived = _stream_function_derived_fields(
            field.psi, field.u, field.v, grid.dx, grid.dy
        )

        l2_errors = {
            quantity: l2_error(derived[quantity], getattr(field, quantity))
            for quantity in QUANTITIES
        }
        linf_errors = {
            quantity: linf_error(derived[quantity], getattr(field, quantity))
            for quantity in QUANTITIES
        }
        levels.append(
            ConvergenceLevel(n=n, h=grid.dx, l2_errors=l2_errors, linf_errors=linf_errors)
        )

    h_values = np.array([level.h for level in levels], dtype=np.float64)
    l2_orders = {
        quantity: estimate_order(
            h_values, np.array([level.l2_errors[quantity] for level in levels])
        )
        for quantity in QUANTITIES
    }
    linf_orders = {
        quantity: estimate_order(
            h_values, np.array([level.linf_errors[quantity] for level in levels])
        )
        for quantity in QUANTITIES
    }

    return ConvergenceReport(levels=levels, l2_orders=l2_orders, linf_orders=linf_orders)


__all__ = [
    "QUANTITIES",
    "ConvergenceLevel",
    "ConvergenceReport",
    "run_convergence_study",
]
