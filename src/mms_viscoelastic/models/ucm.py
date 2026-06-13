"""Upper-Convected Maxwell (UCM) constitutive model.

This module implements the dimensionless UCM equation

    Wi * (D tau / Dt) + tau = 2 (1 - beta) D,

where ``tau`` is the polymeric extra-stress tensor, ``D`` is the
rate-of-deformation tensor, ``Wi`` is the Weissenberg number, ``beta`` is the
solvent-viscosity ratio and ``D tau / Dt`` is the upper-convected derivative

    D tau / Dt = d tau/dt + u . grad(tau) - L.tau - tau.L^T,

with ``L`` the velocity-gradient tensor ``L_ij = d u_i / d x_j``. See
``docs/constitutive-models.md`` for the full derivation and the closed-form
steady simple-shear (Couette) solution used in ``ucm_steady_simple_shear``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mms_viscoelastic.fields import FloatArray


@dataclass(frozen=True)
class VelocityGradient2D:
    """Components of the 2D velocity-gradient tensor ``L_ij = d u_i / d x_j``."""

    dudx: FloatArray
    dudy: FloatArray
    dvdx: FloatArray
    dvdy: FloatArray


@dataclass(frozen=True)
class StressTensor2D:
    """Symmetric components of a 2D extra-stress (or rate-of-deformation) tensor."""

    txx: FloatArray
    txy: FloatArray
    tyy: FloatArray


def rate_of_deformation_tensor(velocity_gradient: VelocityGradient2D) -> StressTensor2D:
    """Return the rate-of-deformation tensor ``D = (L + L^T) / 2``."""

    return StressTensor2D(
        txx=np.asarray(velocity_gradient.dudx, dtype=np.float64),
        txy=0.5
        * (
            np.asarray(velocity_gradient.dudy, dtype=np.float64)
            + np.asarray(velocity_gradient.dvdx, dtype=np.float64)
        ),
        tyy=np.asarray(velocity_gradient.dvdy, dtype=np.float64),
    )


def upper_convected_derivative(
    stress: StressTensor2D,
    velocity_gradient: VelocityGradient2D,
    *,
    material_derivative: StressTensor2D | None = None,
) -> StressTensor2D:
    """Return the upper-convected derivative ``D tau / Dt = d tau/dt + u.grad(tau) - L.tau - tau.L^T``.

    ``material_derivative`` should already contain the transport terms
    ``d tau/dt + u . grad(tau)``. It defaults to zero, which is exact for
    steady flows with a spatially uniform stress field (e.g. simple-shear /
    Couette flow), where only the ``- L.tau - tau.L^T`` contribution survives.
    """

    l11 = np.asarray(velocity_gradient.dudx, dtype=np.float64)
    l12 = np.asarray(velocity_gradient.dudy, dtype=np.float64)
    l21 = np.asarray(velocity_gradient.dvdx, dtype=np.float64)
    l22 = np.asarray(velocity_gradient.dvdy, dtype=np.float64)
    txx, txy, tyy = stress.txx, stress.txy, stress.tyy

    convective_xx = 2.0 * (l11 * txx + l12 * txy)
    convective_xy = l21 * txx + (l11 + l22) * txy + l12 * tyy
    convective_yy = 2.0 * (l21 * txy + l22 * tyy)

    if material_derivative is None:
        zero = np.zeros_like(np.asarray(txx, dtype=np.float64))
        material_derivative = StressTensor2D(txx=zero, txy=zero, tyy=zero)

    return StressTensor2D(
        txx=material_derivative.txx - convective_xx,
        txy=material_derivative.txy - convective_xy,
        tyy=material_derivative.tyy - convective_yy,
    )


def ucm_residual(
    stress: StressTensor2D,
    velocity_gradient: VelocityGradient2D,
    *,
    wi: float,
    beta: float,
    material_derivative: StressTensor2D | None = None,
) -> StressTensor2D:
    """Return the residual ``Wi * (D tau / Dt) + tau - 2 (1 - beta) D`` of the UCM equation.

    For an exact UCM solution this residual is (numerically) zero.
    """

    if wi <= 0:
        raise ValueError("wi must be positive.")
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must satisfy 0 <= beta <= 1.")

    convected = upper_convected_derivative(
        stress, velocity_gradient, material_derivative=material_derivative
    )
    deformation = rate_of_deformation_tensor(velocity_gradient)
    polymeric_viscosity = 1.0 - beta

    return StressTensor2D(
        txx=wi * convected.txx + stress.txx - 2.0 * polymeric_viscosity * deformation.txx,
        txy=wi * convected.txy + stress.txy - 2.0 * polymeric_viscosity * deformation.txy,
        tyy=wi * convected.tyy + stress.tyy - 2.0 * polymeric_viscosity * deformation.tyy,
    )


def ucm_steady_simple_shear(
    shear_rate: FloatArray | float, *, wi: float, beta: float
) -> StressTensor2D:
    """Return the closed-form steady UCM stress for simple-shear flow ``u = shear_rate * y``.

        txy = (1 - beta) * shear_rate
        txx = 2 * Wi * (1 - beta) * shear_rate ** 2
        tyy = 0

    See ``docs/constitutive-models.md`` for the derivation.
    """

    if wi <= 0:
        raise ValueError("wi must be positive.")
    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must satisfy 0 <= beta <= 1.")

    gamma = np.asarray(shear_rate, dtype=np.float64)
    polymeric_viscosity = 1.0 - beta

    return StressTensor2D(
        txx=2.0 * wi * polymeric_viscosity * gamma**2,
        txy=polymeric_viscosity * gamma,
        tyy=np.zeros_like(gamma),
    )


__all__ = [
    "StressTensor2D",
    "VelocityGradient2D",
    "rate_of_deformation_tensor",
    "ucm_residual",
    "ucm_steady_simple_shear",
    "upper_convected_derivative",
]
