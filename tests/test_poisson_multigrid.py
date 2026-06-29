"""Tests for the 4th-order compact multigrid Poisson solver (Etapa 17, ``poi_psi4a``)."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.numerics.errors import estimate_order
from mms_viscoelastic.numerics.poisson_multigrid import solve_poisson_multigrid_4th_order


def test_solve_poisson_multigrid_preserves_boundary() -> None:
    grid = Grid2D.from_bounds(nx=17, ny=17)
    fields = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters())
    assert fields.wz is not None
    assert fields.psi is not None

    psi_numerical = solve_poisson_multigrid_4th_order(fields.wz, grid, fields.psi)

    assert np.array_equal(psi_numerical[0, :], fields.psi[0, :])
    assert np.array_equal(psi_numerical[-1, :], fields.psi[-1, :])
    assert np.array_equal(psi_numerical[:, 0], fields.psi[:, 0])
    assert np.array_equal(psi_numerical[:, -1], fields.psi[:, -1])


@pytest.mark.parametrize("nx", [9, 17, 33])
def test_solve_poisson_multigrid_zero_source_is_boundary_harmonic_extension(nx: int) -> None:
    grid = Grid2D.from_bounds(nx=nx, ny=nx)
    zero = np.zeros((nx, nx))
    boundary = np.zeros((nx, nx))
    boundary[0, :] = 1.0

    psi = solve_poisson_multigrid_4th_order(zero, grid, boundary)

    assert np.array_equal(psi[0, :], boundary[0, :])
    assert np.all(np.isfinite(psi))


def test_solve_poisson_multigrid_converges_at_fourth_order() -> None:
    grid_sizes = [8, 16, 32, 64]
    linf_errors = []
    for n in grid_sizes:
        grid = Grid2D.from_bounds(nx=n + 1, ny=n + 1)
        fields = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters())
        assert fields.wz is not None
        assert fields.psi is not None

        psi_numerical = solve_poisson_multigrid_4th_order(fields.wz, grid, fields.psi)
        linf_errors.append(np.max(np.abs(psi_numerical - fields.psi)))

    h_values = np.array([1.0 / n for n in grid_sizes])
    order = estimate_order(h_values, np.array(linf_errors))

    assert order > 3.8
