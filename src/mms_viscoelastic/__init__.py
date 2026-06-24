"""Python foundation for MMS-based verification of viscoelastic flow solvers.

The package starts as a faithful, testable Python porting layer around the
original Fortran thesis code.  The original Fortran files are preserved in the
repository and should be treated as the reference implementation during the
porting process.
"""

from ._version import __version__
from .fields import Field2D
from .grids import Grid2D
from .parameters import FlowParameters

__all__ = [
    "__version__",
    "Field2D",
    "FlowParameters",
    "Grid2D",
]
