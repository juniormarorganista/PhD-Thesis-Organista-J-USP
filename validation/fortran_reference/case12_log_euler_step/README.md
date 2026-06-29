# Case 12 log-conformation explicit Euler step Fortran reference

This directory contains a standalone Fortran program that exports reference
values for one explicit Euler step of the MMS Case 12 log-conformation
formulation (``ty_sim == 22``): ``wz``, ``Psixx``, ``Psixy``, ``Psiyy``,
reported in the output CSV as the recovered ``(Txx, Txy, Tyy)`` for direct
comparison with `validation/fortran_reference/case12_euler_step/`.

This path is **not** discussed in the thesis text (zero occurrences of
"conformation" in `Thesis/tex/`); it is ported only for source-code fidelity
with `src/fortran_original/4Order_LPOG_Euler/` (Etapa 15, see
`docs/case-12-log-conformation.md`).

The right-hand sides mirror the original routines:

```text
src/fortran_original/4Order_LPOG_Euler/conform.f::decomp_matrix
src/fortran_original/4Order_LPOG_Euler/conform.f::multi_matrix
src/fortran_original/4Order_LPOG_Euler/conform.f::T_to_Psi
src/fortran_original/4Order_LPOG_Euler/conform.f::Psi_to_T
src/fortran_original/4Order_LPOG_Euler/conform.f::tf_conformation
src/fortran_original/4Order_LPOG_Euler/conform.f::kernel_conformation
src/fortran_original/4Order_LPOG_Euler/drvs.f::drv_Psi
src/fortran_original/4Order_LPOG_Euler/drvs.f::nlterms_logPsi
```

and the Euler update performed by the `ty_sim == 22` branch of:

```text
src/fortran_original/4Order_LPOG_Euler/euler_rk.f::solver_log_euler
```

The original Fortran thesis code is not modified.

## Scope

- `decomp_matrix` and `multi_matrix` are ported verbatim as the generic
  Jacobi-rotation algorithm (Numerical Recipes `jacobi`), specialized to
  `dimen = 2` (the only size the original code ever instantiates), rather
  than the closed-form rotation used by
  `mms_viscoelastic.models.log_conformation.decompose_symmetric_2x2`. This
  keeps the Fortran reference algorithmically independent from the Python
  implementation it validates.
- Derivatives use the same second-order finite-difference operators as
  `mms_viscoelastic.numerics.derivatives` (central differences with
  one-sided edges for first derivatives, 3-point central with zero boundary
  for second derivatives), **not** the 4th-order LPOG compact operators nor
  the Lele compact filter, mirroring the scope of
  `case12_euler_step`.
- Only the explicit-Euler variant (`ty_sim == 22`, `solver_log_euler`) is
  ported. The RK4 variant (`ty_sim == 2`, `solver_log_rk4`) is deferred —
  it re-evaluates the boundary conditions and `cc_Psi`/`tf_conformation`
  between stages, which is out of scope for this single-step reference.
- `cc_Psi` (the boundary-specific `Psi` re-projection) is not ported: only
  interior grid points are advanced; boundary values are copied unchanged,
  as in `case12_euler_step`.
- Default parameters (`Re=1`, `Wi=1`, `beta=0.5`, `alphaG=0`, `dt=1.0e-5`)
  match `mms_viscoelastic.numerics.time_stepping.EulerStepParameters()`.

## Generate the input CSV

The input CSV provides the Case 12 exact fields and source terms from the
already-validated Python implementations, evaluated at `t=0.3` (not `t=0`:
at `t=0` the Case 12 polymeric stress is identically zero everywhere, which
makes the conformation matrix the identity with repeated eigenvalues and
exercises neither the eigendecomposition nor the rotation/kernel terms):

```powershell
python validation\fortran_reference\case12_log_euler_step\generate_input_csv.py `
  validation\reference_data\case12_log_euler_step_input_9x9.csv 9 9 0.3
```

## Compile

From the repository root:

```powershell
gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.f90 `
  -o validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.exe
```

## Generate reference CSV

Default `9 x 9` grid with `dt = 1.0e-5`:

```powershell
.\validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.exe `
  validation\reference_data\case12_log_euler_step_input_9x9.csv `
  validation\reference_data\case12_log_euler_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5
```

## CSV schema

Input (`case12_log_euler_step_input_9x9.csv`), identical to
`case12_euler_step_input_9x9.csv`:

```text
x,y,u,v,wz,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy
```

Output (`case12_log_euler_step_output_9x9_dt1e-05.csv`):

```text
x,y,wz,txx,txy,tyy
```

`wz` is the field after one explicit Euler step in the log-conformation
formulation; `txx`, `txy`, `tyy` are the polymeric stress recovered
(`Psi_to_T`) from the post-step `Psi`. Interior points only; boundary values
equal the input.

## Validate Python against reference

```powershell
python -m pytest tests\test_case12_log_conformation_fortran_reference.py
```
