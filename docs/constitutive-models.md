# Constitutive models: UCM and Oldroyd-B

This document derives the equations implemented in
[`src/mms_viscoelastic/models/ucm.py`](../src/mms_viscoelastic/models/ucm.py) and
[`src/mms_viscoelastic/models/oldroyd_b.py`](../src/mms_viscoelastic/models/oldroyd_b.py),
and the closed-form steady simple-shear (Couette) solutions used to validate
them in
[`tests/test_model_ucm.py`](../tests/test_model_ucm.py) and
[`tests/test_model_oldroyd_b.py`](../tests/test_model_oldroyd_b.py).

These are the first two constitutive models with real physics in the Python
port; the README's
[generic constitutive structure](../README.md#generic-viscoelastic-constitutive-structure)

```math
\mathrm{Wi}\,\mathcal{D}\boldsymbol{\tau}
+
\mathcal{F}(\boldsymbol{\tau}, \nabla \mathbf{u})
=
2(1-\beta)\mathbf{D}
```

is here specialized first to UCM (`\mathcal{F}(\boldsymbol{\tau}, \nabla
\mathbf{u}) = \boldsymbol{\tau}`) and then to Oldroyd-B, which adds a
Newtonian solvent contribution on top of a UCM polymeric stress.

## 1. The velocity-gradient and rate-of-deformation tensors

For a 2D velocity field $\mathbf{u} = (u, v)$, the velocity-gradient tensor
$\mathbf{L} = \nabla\mathbf{u}$ is defined componentwise as $L_{ij} =
\partial u_i / \partial x_j$, i.e.

```math
\mathbf{L}
=
\begin{bmatrix}
\dfrac{\partial u}{\partial x} & \dfrac{\partial u}{\partial y} \\[6pt]
\dfrac{\partial v}{\partial x} & \dfrac{\partial v}{\partial y}
\end{bmatrix}.
```

This is exactly the `VelocityGradient2D` dataclass in
[`models/ucm.py`](../src/mms_viscoelastic/models/ucm.py): fields `dudx,
dudy, dvdx, dvdy` correspond to $L_{11}, L_{12}, L_{21}, L_{22}$.

The rate-of-deformation tensor is the symmetric part of $\mathbf{L}$,

```math
\mathbf{D} = \frac{1}{2}\left(\mathbf{L} + \mathbf{L}^{T}\right)
=
\begin{bmatrix}
D_{xx} & D_{xy} \\
D_{xy} & D_{yy}
\end{bmatrix}
=
\begin{bmatrix}
\dfrac{\partial u}{\partial x}
&
\dfrac{1}{2}\left(\dfrac{\partial u}{\partial y} + \dfrac{\partial v}{\partial x}\right)
\\[6pt]
\dfrac{1}{2}\left(\dfrac{\partial u}{\partial y} + \dfrac{\partial v}{\partial x}\right)
&
\dfrac{\partial v}{\partial y}
\end{bmatrix},
```

matching the definition already in the README's
[Rate-of-deformation tensor](../README.md#rate-of-deformation-tensor) section
and implemented by `rate_of_deformation_tensor`.

## 2. The upper-convected derivative

The upper-convected (Oldroyd) derivative of a tensor $\boldsymbol{\tau}$ is

```math
\overset{\nabla}{\boldsymbol{\tau}}
=
\frac{D\boldsymbol{\tau}}{Dt}
-
\mathbf{L}\boldsymbol{\tau}
-
\boldsymbol{\tau}\mathbf{L}^{T},
\qquad
\frac{D\boldsymbol{\tau}}{Dt}
=
\frac{\partial \boldsymbol{\tau}}{\partial t}
+
\mathbf{u}\cdot\nabla\boldsymbol{\tau}.
```

This is the standard frame-invariant derivative used for the polymeric
extra-stress in Maxwell-type models (Bird, Armstrong & Hassager, *Dynamics of
Polymeric Liquids*, Vol. 1, 2nd ed., Wiley, 1987, Ch. 9; Oldroyd, *Proc. R.
Soc. Lond. A* 200 (1950) 523-541).

Writing $\boldsymbol{\tau} = \begin{bmatrix}T_{xx} & T_{xy}\\T_{xy} &
T_{yy}\end{bmatrix}$ and $\mathbf{L} = \begin{bmatrix}L_{11} &
L_{12}\\L_{21} & L_{22}\end{bmatrix}$, direct matrix multiplication gives

```math
\mathbf{L}\boldsymbol{\tau} + \boldsymbol{\tau}\mathbf{L}^{T}
=
\begin{bmatrix}
2\left(L_{11}T_{xx} + L_{12}T_{xy}\right)
&
L_{21}T_{xx} + (L_{11}+L_{22})T_{xy} + L_{12}T_{yy}
\\
L_{21}T_{xx} + (L_{11}+L_{22})T_{xy} + L_{12}T_{yy}
&
2\left(L_{21}T_{xy} + L_{22}T_{yy}\right)
\end{bmatrix}.
```

This componentwise formula is exactly what
`upper_convected_derivative` evaluates: it takes an (optional)
`material_derivative` for the $D\boldsymbol{\tau}/Dt$ term (defaulting to
zero) and subtracts the matrix above.

> **Steady, spatially uniform stress.** If $\boldsymbol{\tau}$ is steady
> ($\partial \boldsymbol{\tau}/\partial t = 0$) and spatially uniform
> ($\nabla\boldsymbol{\tau} = 0$, so $\mathbf{u}\cdot\nabla\boldsymbol{\tau} =
> 0$), then $D\boldsymbol{\tau}/Dt = 0$ exactly, and
> $\overset{\nabla}{\boldsymbol{\tau}} = -(\mathbf{L}\boldsymbol{\tau} +
> \boldsymbol{\tau}\mathbf{L}^{T})$. This is the regime of the simple-shear
> (Couette) solutions derived below, so `material_derivative=None` (the
> default) is exact for them - not an approximation.

## 3. The UCM equation

The dimensionless Upper-Convected Maxwell (UCM) equation for the polymeric
extra-stress $\boldsymbol{\tau}$ is

```math
\mathrm{Wi}\,\overset{\nabla}{\boldsymbol{\tau}} + \boldsymbol{\tau} = 2(1-\beta)\mathbf{D},
```

where $\mathrm{Wi} > 0$ is the Weissenberg number (dimensionless relaxation
time) and $0 \le \beta \le 1$ is the solvent-viscosity ratio, so $(1-\beta)$
is the dimensionless polymeric viscosity. `ucm_residual` evaluates the
left-hand side minus the right-hand side of this equation; for an exact UCM
solution the residual is zero (to floating-point round-off).

### 3.1 Steady simple shear (Couette flow)

Consider plane Couette flow, $u(y) = \dot\gamma\, y$, $v = 0$, with constant
shear rate $\dot\gamma$. Then

```math
\mathbf{L}
=
\begin{bmatrix}0 & \dot\gamma \\ 0 & 0\end{bmatrix},
\qquad
\mathbf{D}
=
\begin{bmatrix}0 & \dot\gamma/2 \\ \dot\gamma/2 & 0\end{bmatrix}.
```

Seek a steady, spatially uniform solution $\boldsymbol{\tau} =
\begin{bmatrix}T_{xx} & T_{xy}\\T_{xy} & T_{yy}\end{bmatrix}$ with constant
entries. By the remark in Section 2,
$\overset{\nabla}{\boldsymbol{\tau}} = -(\mathbf{L}\boldsymbol{\tau} +
\boldsymbol{\tau}\mathbf{L}^{T})$, and substituting $L_{11}=L_{22}=L_{21}=0$,
$L_{12}=\dot\gamma$ into the formula of Section 2 gives

```math
\mathbf{L}\boldsymbol{\tau} + \boldsymbol{\tau}\mathbf{L}^{T}
=
\begin{bmatrix}
2\dot\gamma T_{xy} & \dot\gamma T_{yy} \\
\dot\gamma T_{yy} & 0
\end{bmatrix}.
```

The UCM equation $\mathrm{Wi}\,\overset{\nabla}{\boldsymbol{\tau}} +
\boldsymbol{\tau} = 2(1-\beta)\mathbf{D}$ then decouples into three scalar
equations:

```math
T_{xx} - 2\,\mathrm{Wi}\,\dot\gamma\, T_{xy} = 0,
\qquad
T_{xy} - \mathrm{Wi}\,\dot\gamma\, T_{yy} = (1-\beta)\dot\gamma,
\qquad
T_{yy} = 0.
```

Solving from the last equation upward gives the closed-form steady solution

```math
\boxed{
T_{yy} = 0,
\qquad
T_{xy} = (1-\beta)\,\dot\gamma,
\qquad
T_{xx} = 2\,\mathrm{Wi}\,(1-\beta)\,\dot\gamma^{2}.
}
```

This is the textbook UCM simple-shear result: the shear stress $T_{xy}$ is
identical to the Newtonian (solvent-only) value scaled by the polymeric
viscosity $(1-\beta)$, the second normal-stress difference $T_{yy}$ vanishes,
and the first normal-stress difference $N_1 = T_{xx} - T_{yy} = T_{xx} =
2\,\mathrm{Wi}\,(1-\beta)\,\dot\gamma^2$ grows quadratically with the shear
rate (Bird, Armstrong & Hassager, 1987, Eq. 7.3-7.5 for the dimensional
analogue $N_1 = 2\eta_p\lambda\dot\gamma^2$, $\tau_{xy} = \eta_p\dot\gamma$,
$N_2 = 0$). This is exactly what `ucm_steady_simple_shear` returns, and
`tests/test_model_ucm.py` checks both the closed-form values for a grid of
$(\mathrm{Wi}, \beta, \dot\gamma)$ and that substituting the result back into
`ucm_residual` (with the corresponding simple-shear velocity gradient)
returns zero to within `1e-12`.

## 4. The Oldroyd-B equation

Oldroyd-B decomposes the total extra-stress into a polymeric part obeying
UCM and a Newtonian solvent part:

```math
\boldsymbol{\tau} = \boldsymbol{\tau}_p + \boldsymbol{\tau}_s,
\qquad
\boldsymbol{\tau}_s = 2\beta\,\mathbf{D},
\qquad
\mathrm{Wi}\,\overset{\nabla}{\boldsymbol{\tau}_p} + \boldsymbol{\tau}_p = 2(1-\beta)\mathbf{D}.
```

This is the split implemented by `solvent_stress` (returns
$\boldsymbol{\tau}_s$), `total_stress` (adds $\boldsymbol{\tau}_s$ to a given
polymeric stress) and `oldroyd_b_residual` (delegates to `ucm_residual` on
$\boldsymbol{\tau}_p$, since $\boldsymbol{\tau}_s$ satisfies its own
constitutive relation algebraically and contributes no residual).

> **Equivalence with the single-equation form.** Substituting
> $\boldsymbol{\tau} = \boldsymbol{\tau}_p + 2\beta\mathbf{D}$ into the
> polymeric UCM equation recovers the classical Oldroyd-B single-equation
> form $\boldsymbol{\tau} + \mathrm{Wi}\,\overset{\nabla}{\boldsymbol{\tau}} =
> 2\mathbf{D} + 2\beta\,\mathrm{Wi}\,\overset{\nabla}{\mathbf{D}}$ for steady
> homogeneous flows (where $\overset{\nabla}{\mathbf{D}} =
> -(\mathbf{L}\mathbf{D}+\mathbf{D}\mathbf{L}^{T})$, by the same argument as
> Section 2). The split form is preferred here because it reuses the UCM
> operators directly and keeps the Newtonian limit explicit.

### 4.1 Steady simple shear (Couette flow)

For the same Couette profile $u(y)=\dot\gamma y$, $v=0$, the solvent stress
is

```math
\boldsymbol{\tau}_s = 2\beta\mathbf{D} =
\begin{bmatrix}0 & \beta\dot\gamma \\ \beta\dot\gamma & 0\end{bmatrix},
```

and the polymeric stress $\boldsymbol{\tau}_p$ is given by Section 3.1 with
$(1-\beta)$ as its viscosity factor. Adding the two gives the total
Oldroyd-B stress

```math
\boxed{
T_{yy} = 0,
\qquad
T_{xy} = (1-\beta)\,\dot\gamma + \beta\,\dot\gamma = \dot\gamma,
\qquad
T_{xx} = 2\,\mathrm{Wi}\,(1-\beta)\,\dot\gamma^{2}.
}
```

Two consistency checks follow directly from these formulas and are exercised
in `tests/test_model_oldroyd_b.py`:

- **$\beta = 0$** (no solvent): $\boldsymbol{\tau}_s = 0$ and the Oldroyd-B
  solution reduces to the pure UCM solution of Section 3.1.
- **$\beta = 1$** (no polymer): $(1-\beta) = 0$, so $\boldsymbol{\tau}_p = 0$
  and $T_{xy} = \dot\gamma$, $T_{xx}=T_{yy}=0$ - the purely Newtonian
  (solvent-only) simple-shear stress, with total shear stress equal to the
  shear rate because the dimensionless total viscosity is normalized to $1$.

Note that the *total* shear stress $T_{xy} = \dot\gamma$ is independent of
both $\mathrm{Wi}$ and $\beta$: redistributing viscosity between solvent and
polymer does not change the total shear stress in steady simple shear, only
the first normal-stress difference $N_1 = T_{xx} = 2\,\mathrm{Wi}\,(1-\beta)\,
\dot\gamma^2$, which vanishes in the Newtonian limit $\beta \to 1$ and in the
inertialess/zero-relaxation limit $\mathrm{Wi}\to 0$.

## 5. Reproducing

```powershell
python -m pytest tests\test_model_ucm.py tests\test_model_oldroyd_b.py -v
```

Both test modules check, for a grid of $(\mathrm{Wi}, \beta, \dot\gamma)$
values (including the Newtonian limits $\beta=0$ and $\beta=1$, and
array-valued shear rates):

- the closed-form boxed formulas above;
- that the closed-form stresses, combined with the simple-shear velocity
  gradient, drive `ucm_residual` / `oldroyd_b_residual` to zero (within
  `1e-12`);
- parameter validation (`Wi > 0`, `0 <= beta <= 1`).
