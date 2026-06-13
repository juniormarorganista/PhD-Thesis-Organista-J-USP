"""Tests for the Oldroyd-B constitutive model."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.models.oldroyd_b import (
    oldroyd_b_residual,
    oldroyd_b_steady_simple_shear,
    solvent_stress,
    total_stress,
)
from mms_viscoelastic.models.ucm import VelocityGradient2D, ucm_steady_simple_shear


def _simple_shear_velocity_gradient(shear_rate: float) -> VelocityGradient2D:
    gamma = np.asarray(shear_rate, dtype=np.float64)
    zero = np.zeros_like(gamma)
    return VelocityGradient2D(dudx=zero, dudy=gamma, dvdx=zero, dvdy=zero)


@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.0, 1.0, 2.5])
def test_solvent_stress_for_simple_shear(beta: float, shear_rate: float) -> None:
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)

    tau_s = solvent_stress(velocity_gradient, beta=beta)

    assert np.isclose(tau_s.txx, 0.0)
    assert np.isclose(tau_s.txy, beta * shear_rate)
    assert np.isclose(tau_s.tyy, 0.0)


@pytest.mark.parametrize("wi", [0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.0, 0.5, 1.0, 3.0])
def test_steady_simple_shear_closed_form(wi: float, beta: float, shear_rate: float) -> None:
    stress = oldroyd_b_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    assert np.isclose(stress.txy, shear_rate)
    assert np.isclose(stress.txx, 2.0 * wi * (1.0 - beta) * shear_rate**2)
    assert np.isclose(stress.tyy, 0.0)


@pytest.mark.parametrize("wi", [0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.5, 1.0, 3.0])
def test_steady_simple_shear_equals_polymeric_plus_solvent(
    wi: float, beta: float, shear_rate: float
) -> None:
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)
    polymeric = ucm_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    expected = total_stress(polymeric, velocity_gradient, beta=beta)
    actual = oldroyd_b_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    assert np.isclose(actual.txx, expected.txx)
    assert np.isclose(actual.txy, expected.txy)
    assert np.isclose(actual.tyy, expected.tyy)


@pytest.mark.parametrize("wi", [0.1, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("beta", [0.0, 0.3, 0.5, 1.0])
@pytest.mark.parametrize("shear_rate", [0.5, 1.0, 3.0])
def test_polymeric_stress_satisfies_oldroyd_b_residual(
    wi: float, beta: float, shear_rate: float
) -> None:
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)
    polymeric = ucm_steady_simple_shear(shear_rate, wi=wi, beta=beta)

    residual = oldroyd_b_residual(polymeric, velocity_gradient, wi=wi, beta=beta)

    assert np.isclose(residual.txx, 0.0, atol=1e-12)
    assert np.isclose(residual.txy, 0.0, atol=1e-12)
    assert np.isclose(residual.tyy, 0.0, atol=1e-12)


def test_beta_zero_reduces_to_pure_ucm() -> None:
    shear_rate = 1.5
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)

    oldroyd_b = oldroyd_b_steady_simple_shear(shear_rate, wi=0.7, beta=0.0)
    ucm = ucm_steady_simple_shear(shear_rate, wi=0.7, beta=0.0)
    tau_s = solvent_stress(velocity_gradient, beta=0.0)

    assert np.isclose(tau_s.txy, 0.0)
    assert np.isclose(oldroyd_b.txx, ucm.txx)
    assert np.isclose(oldroyd_b.txy, ucm.txy)
    assert np.isclose(oldroyd_b.tyy, ucm.tyy)


def test_beta_one_reduces_to_pure_newtonian_solvent() -> None:
    shear_rate = 2.0
    velocity_gradient = _simple_shear_velocity_gradient(shear_rate)

    oldroyd_b = oldroyd_b_steady_simple_shear(shear_rate, wi=1.0, beta=1.0)
    ucm = ucm_steady_simple_shear(shear_rate, wi=1.0, beta=1.0)
    tau_s = solvent_stress(velocity_gradient, beta=1.0)

    assert np.isclose(ucm.txx, 0.0)
    assert np.isclose(ucm.txy, 0.0)
    assert np.isclose(oldroyd_b.txx, 0.0)
    assert np.isclose(oldroyd_b.txy, tau_s.txy)
    assert np.isclose(oldroyd_b.txy, shear_rate)


def test_steady_simple_shear_accepts_array_input() -> None:
    shear_rate = np.array([0.0, 1.0, 2.0, 4.0])

    stress = oldroyd_b_steady_simple_shear(shear_rate, wi=0.5, beta=0.4)

    assert stress.txx.shape == shear_rate.shape
    assert np.allclose(stress.txy, shear_rate)
    assert np.allclose(stress.txx, 2.0 * 0.5 * 0.6 * shear_rate**2)
    assert np.allclose(stress.tyy, 0.0)


@pytest.mark.parametrize("beta", [-0.1, 1.1])
def test_solvent_stress_rejects_beta_out_of_range(beta: float) -> None:
    velocity_gradient = _simple_shear_velocity_gradient(1.0)

    with pytest.raises(ValueError):
        solvent_stress(velocity_gradient, beta=beta)
