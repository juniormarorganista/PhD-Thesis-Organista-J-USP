"""Numerical operators and verification utilities."""

from .convergence import (
    QUANTITIES,
    ConvergenceLevel,
    ConvergenceReport,
    run_convergence_study,
)
from .derivatives import (
    derivative_x,
    derivative_y,
    laplacian,
    second_derivative_x,
    second_derivative_y,
)
from .errors import estimate_order, l2_error, linf_error, relative_l2_error
from .poisson import (
    PoissonConvergenceLevel,
    PoissonConvergenceReport,
    run_case12_poisson_convergence_study,
    solve_poisson_dirichlet,
)
from .time_stepping import (
    EulerStepParameters,
    euler_step_stress,
    euler_step_wz,
    stress_time_derivatives,
    wz_time_derivative,
)

__all__ = [
    "QUANTITIES",
    "ConvergenceLevel",
    "ConvergenceReport",
    "EulerStepParameters",
    "PoissonConvergenceLevel",
    "PoissonConvergenceReport",
    "derivative_x",
    "derivative_y",
    "estimate_order",
    "euler_step_stress",
    "euler_step_wz",
    "l2_error",
    "laplacian",
    "linf_error",
    "relative_l2_error",
    "run_case12_poisson_convergence_study",
    "run_convergence_study",
    "second_derivative_x",
    "second_derivative_y",
    "solve_poisson_dirichlet",
    "stress_time_derivatives",
    "wz_time_derivative",
]
