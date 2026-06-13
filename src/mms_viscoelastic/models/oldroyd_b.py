"""Oldroyd-B constitutive model.

The Oldroyd-B extra-stress tensor splits into a polymeric contribution that
obeys the Upper-Convected Maxwell (UCM) equation and a Newtonian solvent
contribution:

    tau   = tau_p + tau_s,
    tau_s = 2 * beta * D,
    Wi * (D tau_p / Dt) + tau_p = 2 * (1 - beta) * D,

where ``D`` is the rate-of-deformation tensor, ``Wi`` is the Weissenberg
number and ``beta`` is the solvent-viscosity ratio. See
``docs/constitutive-models.md`` for the full derivation and the closed-form
steady simple-shear (Couette) solution used in
``oldroyd_b_steady_simple_shear``.
"""

from __future__ import annotations

import numpy as np

from mms_viscoelastic.fields import FloatArray
from mms_viscoelastic.models.ucm import (
    StressTensor2D,
    VelocityGradient2D,
    rate_of_deformation_tensor,
    ucm_residual,
    ucm_steady_simple_shear,
)


def solvent_stress(velocity_gradient: VelocityGradient2D, *, beta: float) -> StressTensor2D:
    """Return the Newtonian solvent contribution ``tau_s = 2 * beta * D``."""

    if not 0.0 <= beta <= 1.0:
        raise ValueError("beta must satisfy 0 <= beta <= 1.")

    deformation = rate_of_deformation_tensor(velocity_gradient)
    return StressTensor2D(
        txx=2.0 * beta * deformation.txx,
        txy=2.0 * beta * deformation.txy,
        tyy=2.0 * beta * deformation.tyy,
    )


def total_stress(
    polymeric_stress: StressTensor2D,
    velocity_gradient: VelocityGradient2D,
    *,
    beta: float,
) -> StressTensor2D:
    """Return the total Oldroyd-B extra-stress ``tau = tau_p + tau_s``."""

    tau_s = solvent_stress(velocity_gradient, beta=beta)
    return StressTensor2D(
        txx=polymeric_stress.txx + tau_s.txx,
        txy=polymeric_stress.txy + tau_s.txy,
        tyy=polymeric_stress.tyy + tau_s.tyy,
    )


def oldroyd_b_residual(
    polymeric_stress: StressTensor2D,
    velocity_gradient: VelocityGradient2D,
    *,
    wi: float,
    beta: float,
    material_derivative: StressTensor2D | None = None,
) -> StressTensor2D:
    """Return the residual of the polymeric UCM equation underlying Oldroyd-B.

    This is identical to :func:`mms_viscoelastic.models.ucm.ucm_residual`
    applied to the polymeric stress ``tau_p``; the Newtonian solvent
    contribution ``tau_s = 2 * beta * D`` satisfies its constitutive relation
    algebraically and contributes no residual.
    """

    return ucm_residual(
        polymeric_stress,
        velocity_gradient,
        wi=wi,
        beta=beta,
        material_derivative=material_derivative,
    )


def oldroyd_b_steady_simple_shear(
    shear_rate: FloatArray | float, *, wi: float, beta: float
) -> StressTensor2D:
    """Return the closed-form steady Oldroyd-B total stress for simple-shear flow.

        txy = shear_rate
        txx = 2 * Wi * (1 - beta) * shear_rate ** 2
        tyy = 0

    See ``docs/constitutive-models.md`` for the derivation.
    """

    gamma = np.asarray(shear_rate, dtype=np.float64)
    polymeric = ucm_steady_simple_shear(gamma, wi=wi, beta=beta)
    velocity_gradient = VelocityGradient2D(
        dudx=np.zeros_like(gamma),
        dudy=gamma,
        dvdx=np.zeros_like(gamma),
        dvdy=np.zeros_like(gamma),
    )
    return total_stress(polymeric, velocity_gradient, beta=beta)


__all__ = [
    "StressTensor2D",
    "VelocityGradient2D",
    "oldroyd_b_residual",
    "oldroyd_b_steady_simple_shear",
    "solvent_stress",
    "total_stress",
]
