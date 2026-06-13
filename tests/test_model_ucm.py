"""Tests for the Upper-Convected Maxwell (UCM) constitutive model."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.models.ucm import (
    VelocityGradient2D,
    rate_of_deformation_tensor,
    ucm_residual,
    ucm_steady_simple_shear,
)


def _simple_shear_velocity_gradient(shear_rate: float) -> VelocityGradient2D:
    gamma = np.asarray(shear_rate, dtype=np.float64)
    zero = np.zeros_like(gamma)
    return VelocityGradient2D(dudx=zero, dudy=gamma, dvdx=zero, dvdy=zero)


def test_rate_of_deformation_for_simple_shear() -> None:
    velocity_gradient = _simple_shear_velocity_gradient(2.0)

    deformation = rate_of_deformation_tensor(velocity_gradient)

    assert np.isclose(deformation.txx, 0.0)
    assert np.isclose(deformation.txy, 1.0)
    assert np.isclose(deformation.tyy, 0.0)


@pytest.mark.parametrize("wi", [0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.0, 0.5, 1.0, 3.0])
def test_steady_simple_shear_closed_form(wi: float, beta: float, shear_rate: float) -> None:
    stress = ucm_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    assert np.isclose(stress.txy, (1.0 - beta) * shear_rate)
    assert np.isclose(stress.txx, 2.0 * wi * (1.0 - beta) * shear_rate**2)
    assert np.isclose(stress.tyy, 0.0)


@pytest.mark.parametrize("wi", [0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.5, 1.0, 3.0])
def test_steady_simple_shear_satisfies_ucm_residual(
    wi: float, beta: float, shear_rate: float
) -> None:
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)
    stress = ucm_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    residual = ucm_residual(stress, velocity_gradient, wi=wi, beta=beta)

    assert np.isclose(residual.txx, 0.0, atol=1e-12)
    assert np.isclose(residual.txy, 0.0, atol=1e-12)
    assert np.isclose(residual.tyy, 0.0, atol=1e-12)


def test_steady_simple_shear_accepts_array_input() -> None:
    shear_rate = np.array([0.0, 1.0, 2.0, 4.0])

    stress = ucm_steady_simple_shear(shear_rate, wi=0.5, beta=0.4)

    assert stress.txx.shape == shear_rate.shape
    assert np.allclose(stress.txy, 0.6 * shear_rate)
    assert np.allclose(stress.txx, 2.0 * 0.5 * 0.6 * shear_rate**2)
    assert np.allclose(stress.tyy, 0.0)


@pytest.mark.parametrize("wi", [0.0, -1.0])
def test_ucm_steady_simple_shear_rejects_non_positive_wi(wi: float) -> None:
    with pytest.raises(ValueError):
        ucm_steady_simple_shear(1.0, wi=wi, beta=0.5)


@pytest.mark.parametrize("beta", [-0.1, 1.1])
def test_ucm_steady_simple_shear_rejects_beta_out_of_range(beta: float) -> None:
    with pytest.raises(ValueError):
        ucm_steady_simple_shear(1.0, wi=1.0, beta=beta)


def test_ucm_residual_rejects_invalid_parameters() -> None:
    velocity_gradient = _simple_shear_velocity_gradient(1.0)
    stress = ucm_steady_simple_shear(1.0, wi=1.0, beta=0.5)

    with pytest.raises(ValueError):
        ucm_residual(stress, velocity_gradient, wi=0.0, beta=0.5)

    with pytest.raises(ValueError):
        ucm_residual(stress, velocity_gradient, wi=1.0, beta=1.5)
