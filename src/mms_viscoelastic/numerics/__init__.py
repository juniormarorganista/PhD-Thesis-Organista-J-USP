"""Numerical operators and verification utilities."""

from .convergence import (
    QUANTITIES,
    ConvergenceLevel,
    ConvergenceReport,
    run_convergence_study,
)
from .derivatives import derivative_x, derivative_y, laplacian
from .errors import estimate_order, l2_error, linf_error, relative_l2_error

__all__ = [
    "QUANTITIES",
    "ConvergenceLevel",
    "ConvergenceReport",
    "derivative_x",
    "derivative_y",
    "estimate_order",
    "l2_error",
    "laplacian",
    "linf_error",
    "relative_l2_error",
    "run_convergence_study",
]
