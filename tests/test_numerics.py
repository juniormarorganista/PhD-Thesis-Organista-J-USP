import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.numerics import derivative_x, derivative_y, estimate_order, l2_error, linf_error


def test_derivatives_for_linear_field() -> None:
    grid = Grid2D.from_bounds(nx=8, ny=6)
    f = 2.0 * grid.X - 3.0 * grid.Y + 1.0

    dfdx = derivative_x(f, grid.dx)
    dfdy = derivative_y(f, grid.dy)

    assert np.allclose(dfdx, 2.0)
    assert np.allclose(dfdy, -3.0)


def test_error_norms_and_order() -> None:
    exact = np.array([1.0, 2.0, 3.0])
    numerical = np.array([1.0, 2.0, 4.0])

    assert l2_error(numerical, exact) > 0.0
    assert linf_error(numerical, exact) == 1.0

    h = np.array([0.5, 0.25, 0.125])
    err = h**2
    assert np.isclose(estimate_order(h, err), 2.0)
