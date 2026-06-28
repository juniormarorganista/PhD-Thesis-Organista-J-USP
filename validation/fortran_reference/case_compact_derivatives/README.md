# Compact 4th-order derivatives Fortran reference

This directory contains a standalone Fortran program that exports reference
values for the fourth-order compact (Padé-type) derivative operators:

```text
compact_derivative_x / compact_derivative_y
compact_second_derivative_x / compact_second_derivative_y
```

The algorithm mirrors the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/derivs.f   (coef_4a, coeff_4a, rhsx_4a, rhsxx_4a, rhsy_4a, rhsyy_4a)
src/fortran_original/4Order_LPOG_Euler/solvers.f  (tridx, tridy)
```

applied to the analytic test function `f(x, y) = sin(pi*x) * cos(pi*y)`,
independent of any MMS case, since this validates the differentiation
operator itself rather than a specific manufactured solution. The original
Fortran thesis code is not modified.

## Boundary-closure truncation order

Substituting the exact solution into the boundary-closure equations and
expanding in Taylor series (verified symbolically with `sympy`):

- First derivative (`f_1' + 3f_2' = (-17f_1+9f_2+9f_3-f_4)/(6h)`): residual
  `(1/20) h^4 f^(5)` -- fourth order, matching the interior stencil.
- Second derivative (`f_1'' + 11f_2'' = (39f_1-81f_2+45f_3-3f_4)/(3h^2)`):
  residual `(1/12) h^3 f^(5)` -- third order, one order below the interior
  `O(h^4)` stencil.

This is a property of the boundary closure itself (both formulas are
transcribed unchanged from `coef_4a`/`coeff_4a`/`rhsx_4a`/`rhsxx_4a`), not a
porting defect. It does not degrade the global convergence order reported
for Case 12, because field boundaries are subsequently overwritten by the
exact manufactured boundary condition at every time step, never read from
this operator's own boundary output. See
`tests/test_compact_derivatives_fortran_reference.py` for the numerical
confirmation (observed boundary order ~3 for the second derivative, ~4
everywhere for the first derivative).

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case_compact_derivatives\case_compact_derivatives_reference.f90 `
  -o validation\fortran_reference\case_compact_derivatives\case_compact_derivatives_reference.exe
```

## Generate reference CSV

```powershell
.\validation\fortran_reference\case_compact_derivatives\case_compact_derivatives_reference.exe `
  validation\reference_data\case_compact_derivatives_9x9_default.csv 9 9

.\validation\fortran_reference\case_compact_derivatives\case_compact_derivatives_reference.exe `
  validation\reference_data\case_compact_derivatives_17x17_default.csv 17 17
```

## Validate Python against reference

```powershell
python -m pytest tests\test_compact_derivatives_fortran_reference.py -v
```
