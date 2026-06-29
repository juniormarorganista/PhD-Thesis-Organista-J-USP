"""Demonstrate Cases 13 and 14 -- the two non-MMS benchmark configurations.

Unlike Cases 01-12, Cases 13 and 14 are not method-of-manufactured-
solutions verification cases: they are the classic and smooth-startup
variants of the lid-driven cavity flow benchmark, used to exercise the
solver's qualitative behavior rather than to verify numerical order of
accuracy. Their exact fields and source terms are identically zero in
the original Fortran code (and have no Mathematica notebook, since
there is nothing to symbolically derive). They are ported only for
structural parity with the other 12 cases.
"""

from __future__ import annotations

import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_13_exact, evaluate_case_14_exact
from mms_viscoelastic.source_terms import evaluate_case_13_source_terms, evaluate_case_14_source_terms


def main() -> None:
    grid = Grid2D.from_bounds(nx=9, ny=9)

    print("MMS Case 13 / 14 - non-MMS benchmark configurations")
    print("(lid-driven cavity flow: Case 13 constant lid, Case 14 smooth-startup lid)")
    print()

    for case_label, evaluate_exact, evaluate_terms in (
        ("Case 13", evaluate_case_13_exact, evaluate_case_13_source_terms),
        ("Case 14", evaluate_case_14_exact, evaluate_case_14_source_terms),
    ):
        field = evaluate_exact(grid, t=0.7)
        terms = evaluate_terms(grid, t=0.7)

        all_field_zero = all(
            np.all(arr == 0.0) for arr in (field.u, field.v, field.wz, field.psi, field.txx, field.txy, field.tyy)
        )
        all_terms_zero = all(np.all(arr == 0.0) for arr in (terms.tfwz, terms.tfTxx, terms.tfTxy, terms.tfTyy))

        print(f"{case_label}: exact fields all zero = {all_field_zero}, source terms all zero = {all_terms_zero}")


if __name__ == "__main__":
    main()
