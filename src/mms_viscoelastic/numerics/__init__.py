"""Numerical operators and verification utilities."""

from .convergence import (
    QUANTITIES,
    ConvergenceLevel,
    ConvergenceReport,
    run_convergence_study,
)
from .derivatives import (
    compact_derivative_x,
    compact_derivative_y,
    compact_second_derivative_x,
    compact_second_derivative_y,
    derivative_x,
    derivative_y,
    laplacian,
    second_derivative_x,
    second_derivative_y,
    solve_tridiagonal,
)
from .errors import estimate_order, l2_error, linf_error, relative_l2_error
from .filters import DEFAULT_ALPHA_F, filter_coefficients, lele_filter_x, lele_filter_y
from .poisson import (
    PoissonConvergenceLevel,
    PoissonConvergenceReport,
    run_case12_poisson_convergence_study,
    solve_poisson_dirichlet,
)
from .poisson_multigrid import (
    DEFAULT_ERROR_TOL,
    DEFAULT_MAX_VCYCLES,
    DEFAULT_N_LEVELS,
    solve_poisson_multigrid_4th_order,
)
from .time_marching import TimeMarchingState, case12_multistep_march, case12_time_step
from .time_stepping import (
    EulerStepParameters,
    euler_step_stress,
    euler_step_wz,
    stress_time_derivatives,
    wz_time_derivative,
)

__all__ = [
    "DEFAULT_ALPHA_F",
    "DEFAULT_ERROR_TOL",
    "DEFAULT_MAX_VCYCLES",
    "DEFAULT_N_LEVELS",
    "QUANTITIES",
    "ConvergenceLevel",
    "ConvergenceReport",
    "EulerStepParameters",
    "PoissonConvergenceLevel",
    "PoissonConvergenceReport",
    "TimeMarchingState",
    "case12_multistep_march",
    "case12_time_step",
    "compact_derivative_x",
    "compact_derivative_y",
    "compact_second_derivative_x",
    "compact_second_derivative_y",
    "derivative_x",
    "derivative_y",
    "estimate_order",
    "euler_step_stress",
    "euler_step_wz",
    "filter_coefficients",
    "l2_error",
    "laplacian",
    "lele_filter_x",
    "lele_filter_y",
    "linf_error",
    "relative_l2_error",
    "run_case12_poisson_convergence_study",
    "run_convergence_study",
    "second_derivative_x",
    "second_derivative_y",
    "solve_poisson_dirichlet",
    "solve_poisson_multigrid_4th_order",
    "solve_tridiagonal",
    "stress_time_derivatives",
    "wz_time_derivative",
]
