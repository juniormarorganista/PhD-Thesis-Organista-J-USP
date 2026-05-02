# Thesis mathematical overview

This document summarizes the core mathematical structure used in the README and
future presentation material.

## Incompressibility

```math
\nabla \cdot \mathbf{u} = 0.
```

For two-dimensional flows:

```math
\mathbf{u} = (u,v).
```

## Vorticity-streamfunction relation

```math
u = \frac{\partial \psi}{\partial y},
\qquad
v = -\frac{\partial \psi}{\partial x}.
```

The scalar vorticity is:

```math
\omega_z =
\frac{\partial v}{\partial x}
-
\frac{\partial u}{\partial y}.
```

The streamfunction is obtained from a Poisson-type equation:

```math
\nabla^2 \psi = -\omega_z.
```

## Generic constitutive structure

A compact representation of the viscoelastic stress evolution can be written as:

```math
\mathrm{Wi}\,\mathcal{D}\boldsymbol{\tau}
+
\mathcal{F}(\boldsymbol{\tau}, \nabla \mathbf{u})
=
2(1-\beta)\mathbf{D},
```

where:

```math
\mathbf{D}
=
\frac{1}{2}
\left(
\nabla\mathbf{u}
+
\nabla\mathbf{u}^{T}
\right).
```

The exact form of `\mathcal{F}` depends on the constitutive model:

- UCM;
- Oldroyd-B;
- Giesekus;
- LPTT.

## Manufactured solution idea

A manufactured exact field is prescribed:

```math
\mathbf{q}_e(x,y,t)
=
\left(
u_e, v_e, \psi_e, \omega_{z,e},
T_{xx,e}, T_{xy,e}, T_{yy,e}
\right).
```

Substituting it into the governing equations produces source terms:

```math
\mathcal{L}(\mathbf{q}_e) = \mathbf{f}_{\mathrm{MMS}}.
```

The numerical solver is then verified by measuring:

```math
\mathbf{e}_h
=
\mathbf{q}_h
-
\mathbf{q}_e.
```

## Observed order of accuracy

For grid spacing `h`, the observed convergence rate is estimated by:

```math
p
=
\frac{
\log\left(E_h/E_{h/2}\right)
}{
\log(2)
}.
```

This is the central diagnostic for code verification.
