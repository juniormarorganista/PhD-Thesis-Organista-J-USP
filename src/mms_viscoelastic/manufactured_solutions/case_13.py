"""Manufactured solution: Fortran ``case_13.f`` exact fields.

Unlike cases 01-12, case 13 is **not** a method-of-manufactured-solutions
verification case: it is the classic lid-driven cavity flow benchmark,
used to exercise the solver's qualitative behavior (vortex formation)
rather than to verify numerical order of accuracy. Its
``exact_solutions_case_13(t)`` subroutine sets every field
(``uxe``, ``uye``, ``wze``, ``psie``, ``Txxe``, ``Txye``, ``Tyye``) to
``0.0`` identically -- there is no analytical solution to derive, which
is also why no Mathematica notebook exists for this case (unlike cases
01-12, which each have a corresponding ``.nb`` file under
``src/mathematica_original/``).

This module is ported only for structural parity with the other 12
cases (same function signature, same ``Field2D`` return type); it has no
verification value beyond confirming the fields are identically zero.
The case's only non-trivial content is its boundary condition (a
constant unit lid velocity, ``ux(i, jmax) = 1.0``, set in
``boundaries_condictions_case_13``), which -- following the same
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
class Case13Parameters:
    """Parameters for case 13. Unused: every field is identically zero."""

    def validate(self) -> None:
        """No parameters to validate; present for API parity with other cases."""
        return


def evaluate_case_13_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case13Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_13.f``.

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
        Case-13 parameters (unused; included for API parity).

    Returns
    -------
    Field2D
        All-zero fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy``
        and ``tyy``.
    """
    params = parameters or Case13Parameters()
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


def case_13_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 13."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
