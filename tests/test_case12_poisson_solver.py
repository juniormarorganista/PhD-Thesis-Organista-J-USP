"""Tests for the direct Poisson solver and its MMS case 12 verification."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions.case_12 import evaluate_case_12_exact
from mms_viscoelastic.numerics import (
    run_case12_poisson_convergence_study,
    solve_poisson_dirichlet,
)

GRID_SIZES = [8, 16, 32, 64]

L2_ORDER_MIN = 1.8
LINF_ORDER_MIN = 1.8


def test_solve_poisson_dirichlet_recovers_known_quadratic() -> None:
    """phi = x^2 + y^2 has laplacian = 4 everywhere; Dirichlet BC = exact phi."""
    grid = Grid2D.from_bounds(nx=9, ny=9)
    phi_exact = grid.X**2 + grid.Y**2
    source = np.full_like(phi_exact, 4.0)

    phi_numerical = solve_poisson_dirichlet(source, grid, phi_exact)

    np.testing.assert_allclose(phi_numerical, phi_exact, atol=1.0e-10)


def test_solve_poisson_dirichlet_rejects_shape_mismatch() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=5)
    source = np.zeros((5, 5))
    boundary_values = np.zeros((4, 5))

    with pytest.raises(ValueError):
        solve_poisson_dirichlet(source, grid, boundary_values)


def test_case12_poisson_convergence_report_has_one_level_per_grid() -> None:
    report = run_case12_poisson_convergence_study(GRID_SIZES)

    assert [level.n for level in report.levels] == GRID_SIZES


def test_case12_poisson_errors_decrease_monotonically() -> None:
    report = run_case12_poisson_convergence_study(GRID_SIZES)

    errors = [level.l2_error for level in report.levels]
    assert all(later < earlier for earlier, later in zip(errors, errors[1:]))


def test_case12_poisson_l2_order_is_second_order() -> None:
    report = run_case12_poisson_convergence_study(GRID_SIZES)

    assert report.l2_order >= L2_ORDER_MIN


def test_case12_poisson_linf_order_is_second_order() -> None:
    report = run_case12_poisson_convergence_study(GRID_SIZES)

    assert report.linf_order >= LINF_ORDER_MIN


def test_case12_poisson_reproduces_psi_at_high_resolution() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_12_exact(grid, t=0.0)
    assert field.wz is not None
    assert field.psi is not None

    psi_numerical = solve_poisson_dirichlet(field.wz, grid, field.psi)

    assert np.max(np.abs(psi_numerical - field.psi)) < 1.0e-3


def test_run_case12_poisson_convergence_study_requires_at_least_two_grids() -> None:
    with pytest.raises(ValueError):
        run_case12_poisson_convergence_study([16])
