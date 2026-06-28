"""Grid-convergence tests for the case 12 compact finite-difference operators."""

from __future__ import annotations

from mms_viscoelastic.numerics import QUANTITIES, run_convergence_study

GRID_SIZES = [8, 16, 32, 64]

# Observed orders for this grid range: L2 in [4.13, 4.49], Linf in [3.68, 4.04]
# (wz lags slightly behind u/v, consistent with the JNNFM paper). Thresholds
# below keep a margin under the worst observed quantity rather than pinning
# the exact value, so the test does not become flaky under minor changes.
L2_ORDER_MIN = 3.9
LINF_ORDER_MIN = 3.4


def test_convergence_report_has_one_level_per_grid() -> None:
    report = run_convergence_study(GRID_SIZES)

    assert [level.n for level in report.levels] == GRID_SIZES
    for level in report.levels:
        assert set(level.l2_errors) == set(QUANTITIES)
        assert set(level.linf_errors) == set(QUANTITIES)


def test_l2_errors_decrease_monotonically() -> None:
    report = run_convergence_study(GRID_SIZES)

    for quantity in QUANTITIES:
        errors = [level.l2_errors[quantity] for level in report.levels]
        assert all(later < earlier for earlier, later in zip(errors, errors[1:]))


def test_l2_order_is_fourth_order() -> None:
    report = run_convergence_study(GRID_SIZES)

    for quantity in QUANTITIES:
        assert report.l2_orders[quantity] >= L2_ORDER_MIN


def test_linf_order_is_close_to_fourth_order() -> None:
    report = run_convergence_study(GRID_SIZES)

    for quantity in QUANTITIES:
        assert report.linf_orders[quantity] >= LINF_ORDER_MIN


def test_run_convergence_study_requires_at_least_two_grids() -> None:
    try:
        run_convergence_study([16])
    except ValueError:
        pass
    else:
        raise AssertionError("Expected ValueError for a single grid size.")
