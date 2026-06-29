"""Tests for the log-conformation explicit Euler time-stepping (Etapa 15)."""

from __future__ import annotations

import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.numerics.time_stepping import EulerStepParameters
from mms_viscoelastic.numerics.time_stepping_log import euler_step_log_conformation


def test_euler_step_log_conformation_preserves_boundary_and_is_finite() -> None:
    n = 9
    grid = Grid2D.from_bounds(nx=n, ny=n)
    rng = np.random.default_rng(5)

    wz = rng.uniform(-1.0, 1.0, (n, n))
    psixx = rng.uniform(-0.5, 0.5, (n, n))
    psixy = rng.uniform(-0.5, 0.5, (n, n))
    psiyy = rng.uniform(-0.5, 0.5, (n, n))
    ux = rng.uniform(-1.0, 1.0, (n, n))
    uy = rng.uniform(-1.0, 1.0, (n, n))
    zero = np.zeros((n, n))
    parameters = EulerStepParameters(reynolds=100.0, wi=0.05, beta=0.5, dt=1.0e-5)

    wz_new, psixx_new, psixy_new, psiyy_new = euler_step_log_conformation(
        wz, psixx, psixy, psiyy, ux, uy, zero, zero, zero, zero, grid, parameters
    )

    for new, original in (
        (wz_new, wz),
        (psixx_new, psixx),
        (psixy_new, psixy),
        (psiyy_new, psiyy),
    ):
        assert np.array_equal(new[0, :], original[0, :])
        assert np.array_equal(new[-1, :], original[-1, :])
        assert np.array_equal(new[:, 0], original[:, 0])
        assert np.array_equal(new[:, -1], original[:, -1])
        assert np.all(np.isfinite(new))
