# MMS Case 12 - Log-conformation reformulation (`conform.f`)

This document describes the log-conformation reformulation of the
polymeric stress, implemented in
[`src/mms_viscoelastic/models/log_conformation.py`](../src/mms_viscoelastic/models/log_conformation.py)
(eigendecomposition, `Psi` conversions, `kernel_conformation`) and
[`src/mms_viscoelastic/numerics/time_stepping_log.py`](../src/mms_viscoelastic/numerics/time_stepping_log.py)
(the `drv_Psi`/`nlterms_logPsi` right-hand sides and the explicit Euler
update), tested in
[`tests/test_model_log_conformation.py`](../tests/test_model_log_conformation.py),
[`tests/test_time_stepping_log.py`](../tests/test_time_stepping_log.py) and
[`tests/test_case12_log_conformation_fortran_reference.py`](../tests/test_case12_log_conformation_fortran_reference.py).

This is **Etapa 15** of `.claude/PLANEJAMENTO.md`.

## 0. Scope note: not discussed in the thesis text

`grep -ri conformation Thesis/tex/` returns zero matches. This path exists in
`src/fortran_original/4Order_LPOG_Euler/` but is not part of the published
thesis results. It is ported here purely for **source-code fidelity** with
the original Fortran project, not to reproduce a result that appears in the
text. Readers looking for the formulation actually discussed in the thesis
should refer to `docs/case-12-euler-step-comparison.md` (the `Txx`/`Txy`/`Tyy`
formulation, `ty_sim == 1`/`12`).

## 1. Why log-conformation: motivation

The raw polymeric-stress formulation tracks `(Txx, Txy, Tyy)` directly. At
high Weissenberg number the corresponding conformation tensor
`A = I + T_p/etap` can lose positive-definiteness under discretization
error, a well-known failure mode in viscoelastic flow simulation. The
log-conformation reformulation (Fattal & Kupferman) instead evolves the
matrix logarithm

```text
Psi = Q . diag(log(eigval)) . Q^T,   A = Q . diag(eigval) . Q^T
```

which is unconstrained (no positive-definiteness constraint to violate) and
recovers `A` (hence `T_p = etap * (A - I)`) by exponentiating back.

## 2. Source formulation: `conform.f` and `drv_Psi`

### 2.1 Eigendecomposition: `decomp_matrix`, specialized to `dimen = 2`

`decomp_matrix` (`conform.f` lines 304-407) is a generic Jacobi-rotation
eigensolver (Numerical Recipes `jacobi`) for a symmetric `dimen x dimen`
matrix, but `dimen = 2` (`par.nn`) everywhere it is called. For a 2x2
symmetric matrix the algorithm performs **exactly one** Jacobi rotation
before the off-diagonal sum vanishes and it returns. Working through that
single rotation in closed form (including the
`dabs(h)+g.eq.dabs(h)` floating-point absorption guard used to select the
small-rotation branch) shows that the eigenvectors are always the columns of
a plane rotation:

```text
Q = [[cos, sin], [-sin, cos]]
```

`decompose_symmetric_2x2` in
[`models/log_conformation.py`](../src/mms_viscoelastic/models/log_conformation.py)
implements this closed form directly (not `numpy.linalg.eigh`, which would
not reproduce the same eigenvector sign/ordering convention bit-for-bit, and
not an iterative Jacobi loop, which is unnecessary for `dimen = 2`).

The Fortran reference program used to validate this port
(Section 5) instead keeps the **generic iterative Jacobi routine**
(specialized only by fixing `dimen = 2` at compile time), so the closed form
and the iterative algorithm are cross-checked against each other rather than
against themselves.

### 2.2 `T_to_Psi` / `Psi_to_T`

```text
A = [[Txx/etap + 1, Txy/etap], [Txy/etap, Tyy/etap + 1]]
A = Q . diag(l1, l2) . Q^T            (decompose_matrix)
Psi = Q . diag(log(l1), log(l2)) . Q^T

Psi = [[Psixx, Psixy], [Psixy, Psiyy]]
Psi = Q . diag(l1, l2) . Q^T          (decompose_matrix)
A = Q . diag(exp(l1), exp(l2)) . Q^T
Txx = etap*(A11 - 1), Txy = etap*A12, Tyy = etap*(A22 - 1)
```

ported as `stress_to_log_conformation` / `log_conformation_to_stress`.

### 2.3 `tf_conformation`

The manufactured source term is carried through unchanged, just rescaled by
`1/etap` (no nonlinear transform is applied to the source term itself):

```text
tfpsixx = tfTxx / etap
tfpsixy = tfTxy / etap
tfpsiyy = tfTyy / etap
```

ported as `source_term_to_log_conformation`.

### 2.4 `kernel_conformation`

Computes the rotation rate of the eigenbasis (`Omegalog`) and the tensors
`B`, `M` that enter the `Psi` transport equation:

```text
gradu = [[dudx, dudy], [dvdx, dvdy]]
M = Q^T . gradu . Q
omegal = (M12*l2 + M21*l1) / (l2 - l1)        [floored at 1e-14 in magnitude]

B = Q . diag(M11, M22) . diag(l1, l2) . diag(1/l1, 1/l2) . Q^T
M(lambda) = (I - Lambda) . (I + alphaG*(Lambda - I)),  Lambda = diag(l1, l2)
M = Q . M(lambda) . diag(1/l1, 1/l2) . Q^T
Omegalog = [Q . [[0, omegal], [-omegal, 0]] . Q^T]_12
```

`kernel_conformation` in `models/log_conformation.py` mirrors this exact
step-by-step matrix composition (rather than an algebraically simplified
closed form) to keep the port mechanically checkable against `conform.f`,
following the same "literal mechanical transcription over manual
simplification" practice used for the MMS source terms in Etapa 14
(`.claude/PLANEJAMENTO.md`, Case 04 lesson).

### 2.5 `drv_Psi` / `nlterms_logPsi`

```text
call Psi_to_T
call kernel_conformation(Omegalog, Bxx, Bxy, Byy, Mxx, Mxy, Myy)
uPsixx = ux*Psixx, vPsixx = uy*Psixx   (and similarly for Psixy, Psiyy)

dvz = <same formula as the non-log wz_time_derivative,
       using the Txx/Txy/Tyy recovered from Psi_to_T>

dvPsixx = -d(uPsixx)/dx - d(vPsixx)/dy + 2*Omegalog*Psixy + 2*Bxx + Mxx/Wi - tfpsixx
dvPsixy = -d(uPsixy)/dx - d(vPsixy)/dy + Omegalog*(Psiyy-Psixx) + 2*Bxy + Mxy/Wi - tfpsixy
dvPsiyy = -d(uPsiyy)/dx - d(vPsiyy)/dy - 2*Omegalog*Psixy + 2*Byy + Myy/Wi - tfpsiyy
```

`log_conformation_time_derivatives`
(`numerics/time_stepping_log.py`) implements this, reusing
`wz_time_derivative` from the non-log
[`time_stepping.py`](../src/mms_viscoelastic/numerics/time_stepping.py)
unchanged for `dvz` (the formula is identical; only the `Txx`/`Txy`/`Tyy`
inputs differ, being recovered from `Psi` rather than tracked directly).

## 3. Discretization and scope

As in Etapa 9, the second-order operators from
[`numerics/derivatives.py`](../src/mms_viscoelastic/numerics/derivatives.py)
are used, **not** the 4th-order LPOG compact operators nor the Lele compact
filter. Porting those for this path is deferred, mirroring the Etapa 9/17
scope decision for the non-log formulation.

Only the **explicit-Euler variant** (`ty_sim == 22`, `solver_log_euler`) is
ported (`euler_step_log_conformation`). The RK4 variant (`ty_sim == 2`,
`solver_log_rk4`) is **deferred**: it re-evaluates the source term, boundary
conditions and `cc_Psi` between each of its four stages, which would require
the full multi-step time-marching machinery (Etapa 10-style boundary
re-imposition), not just a single-step right-hand-side port.

`cc_Psi` (`conform.f` lines 50-122 — re-projects `Psi` specifically on the
boundary rows/columns) is not ported: only interior grid points are
advanced; boundary values are copied unchanged, as in the non-log Euler
step (Etapa 9).

## 4. Parameter choice

Reuses `EulerStepParameters` unchanged (`reynolds`, `wi`, `beta`, `alpha_g`,
`dt`); `etap = (1 - beta) / (reynolds * wi)` and `etapinv = 1/etap` are
derived from it, matching `par.nn`'s `etap`/`etapinv` definitions.

## 5. Fortran reference program

`validation/fortran_reference/case12_log_euler_step/case12_log_euler_step_reference.f90`
is a standalone program (it does not modify or depend on the original
thesis Fortran files) that:

1. Reads the same input CSV schema as `case12_euler_step`
   (`x,y,u,v,wz,txx,txy,tyy,tfwz,tftxx,tftxy,tftyy`), generated by
   `generate_input_csv.py` at **`t = 0.3`** (not `t = 0`: at `t = 0` the
   Case 12 polymeric stress is identically zero everywhere, making the
   conformation matrix the identity with a repeated eigenvalue, which
   exercises neither the eigendecomposition nor the rotation/kernel terms).
2. Converts `(Txx, Txy, Tyy)` to `(Psixx, Psixy, Psiyy)` (`T_to_Psi_point`).
3. Computes `dvz`, `dvPsixx`, `dvPsixy`, `dvPsiyy` using verbatim ports of
   `decomp_matrix`, `multi_matrix`, `Psi_to_T`, `kernel_conformation` (kept
   as the **generic iterative Jacobi algorithm**, specialized only to
   `dimen = 2`, deliberately not the Python closed form, to keep this
   reference algorithmically independent from the code it validates).
4. Advances the interior points by one Euler step of size `dt = 1.0e-5`.
5. Converts the post-step `Psi` back to `(Txx, Txy, Tyy)` (`Psi_to_T_point`)
   and writes an output CSV with columns `x,y,wz,txx,txy,tyy`, reusing the
   same schema as `case12_euler_step` for direct comparison.

See
[`validation/fortran_reference/case12_log_euler_step/README.md`](../validation/fortran_reference/case12_log_euler_step/README.md)
for compile and run instructions.

## 6. Observed results

`compare_case12_log_euler_step_reference` compares each of `wz`, `Txx`,
`Txy`, `Tyy` after one log-conformation Euler step (`dt = 1.0e-5`, `9 x 9`
grid, `t = 0.3`) between Python and the Fortran reference:

```text
  wz: max_abs_error=4.440892e-16 max_rel_error=2.220446e-16 passed=True
 txx: max_abs_error=1.776357e-15 max_rel_error=1.008560e-15 passed=True
 txy: max_abs_error=2.664535e-15 max_rel_error=9.996688e-16 passed=True
 tyy: max_abs_error=3.108624e-15 max_rel_error=1.348073e-15 passed=True
```

All four components agree to within machine precision, as expected since
both implementations apply the same eigendecomposition, kernel and
Euler-update formulas to identical input fields — the closed-form
eigensolver and the iterative Jacobi rotation agree exactly for `dimen = 2`.
[`tests/test_case12_log_conformation_fortran_reference.py`](../tests/test_case12_log_conformation_fortran_reference.py)
asserts `max_abs_error < 1e-12` and `max_rel_error < 1e-12` for all four
components.

`tests/test_model_log_conformation.py` additionally cross-checks
`decompose_symmetric_2x2` and `kernel_conformation` against independent
pointwise `numpy.linalg.eigh`-based linear algebra (not the Fortran
reference), and verifies the `stress_to_log_conformation` /
`log_conformation_to_stress` round-trip to machine precision over random
positive-definite conformation matrices.

## 7. Run

```powershell
python validation\fortran_reference\case12_log_euler_step\generate_input_csv.py `
  validation\reference_data\case12_log_euler_step_input_9x9.csv 9 9 0.3

gfortran -O2 -std=f2008 `
  validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.f90 `
  -o validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.exe

.\validation\fortran_reference\case12_log_euler_step\case12_log_euler_step_reference.exe `
  validation\reference_data\case12_log_euler_step_input_9x9.csv `
  validation\reference_data\case12_log_euler_step_output_9x9_dt1e-05.csv `
  9 9 1.0e-5

python -m pytest tests\test_case12_log_conformation_fortran_reference.py tests\test_model_log_conformation.py tests\test_time_stepping_log.py
```

## 8. Deferred work

- The RK4 variant (`ty_sim == 2`, `solver_log_rk4`) and the boundary
  re-projection `cc_Psi`, as noted in Section 3.
- A full multi-step `time_marching_log.py` analogous to Etapa 10's
  `case12_time_step`/`case12_multistep_march` (re-imposing boundary
  conditions and recomputing `psi`/`u`/`v` after each Euler step) is not
  built; this Etapa scoped to the single-step right-hand-side/Euler-update
  port only, consistent with the Etapa 9 scope before Etapa 10 extended it.
- Porting the 4th-order LPOG compact derivative operators and the Lele
  compact filter for this path remains out of scope (Etapa 17 precedent).
