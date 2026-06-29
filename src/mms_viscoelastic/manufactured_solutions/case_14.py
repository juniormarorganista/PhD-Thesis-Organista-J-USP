"""Manufactured solution: Fortran ``case_14.f`` exact fields.

Unlike cases 01-12, case 14 is **not** a method-of-manufactured-solutions
verification case: it is a smooth-startup variant of the lid-driven
cavity flow benchmark (the lid velocity ramps up smoothly in time via a
``tanh`` profile, rather than starting at a constant value as in case
13). Its ``exact_solutions_case_14(t)`` subroutine sets every field
(``uxe``, ``uye``, ``wze``, ``psie``, ``Txxe``, ``Txye``, ``Tyye``) to
``0.0`` identically -- there is no analytical solution to derive, which
is also why no Mathematica notebook exists for this case (unlike cases
01-12, which each have a corresponding ``.nb`` file under
``src/mathematica_original/``).

This module is ported only for structural parity with the other 12
cases (same function signature, same ``Field2D`` return type); it has no
verification value beyond confirming the fields are identically zero.
The case's only non-trivial content is its boundary condition (a
spatially-varying, smoothly-ramping lid velocity,
``ux(i, jmax) = 0.8 * (1 - x**2)**2 * x**2 * (1 + tanh(8*t - 4))``, set
in ``boundaries_condictions_case_14``), which -- following the same
convention used for every other case -- is not ported here.

Array convention
----------------
Fortran stores fields as ``field(i, j)``. NumPy stores the same physical
field as arrays with shape ``(ny, nx)``, where axis 0 is ``y`` and axis 1 is
``x``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case14Parameters:
    """Parameters for case 14. Unused: every field is identically zero."""

    def validate(self) -> None:
        """No parameters to validate; present for API parity with other cases."""
        return


def evaluate_case_14_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case14Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_14.f``.

    Every field is identically zero, matching the original Fortran
    subroutine exactly.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid. The returned arrays have shape
        ``(grid.ny, grid.nx)``.
    t:
        Physical time (unused; included for API parity).
    parameters:
        Case-14 parameters (unused; included for API parity).

    Returns
    -------
    Field2D
        All-zero fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy``
        and ``tyy``.
    """
    params = parameters or Case14Parameters()
    params.validate()

    zeros = np.zeros((grid.ny, grid.nx))

    return Field2D(
        u=zeros.copy(),
        v=zeros.copy(),
        psi=zeros.copy(),
        wz=zeros.copy(),
        txx=zeros.copy(),
        txy=zeros.copy(),
        tyy=zeros.copy(),
    )


def case_14_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 14."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
