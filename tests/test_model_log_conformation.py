"""Tests for the log-conformation reformulation (Etapa 15, ``conform.f``)."""

from __future__ import annotations

import numpy as np
import pytest

from mms_viscoelastic.models.log_conformation import (
    decompose_symmetric_2x2,
    kernel_conformation,
    log_conformation_to_stress,
    source_term_to_log_conformation,
    stress_to_log_conformation,
)

ETAP = 0.19


def _random_pd_conformation_stress(rng: np.random.Generator, n: int) -> tuple[np.ndarray, ...]:
    """Return ``(txx, txy, tyy)`` such that ``I + T/etap`` is positive definite."""

    eigs = rng.uniform(0.05, 8.0, (n, 2))
    angles = rng.uniform(0.0, 2.0 * np.pi, n)
    c, s = np.cos(angles), np.sin(angles)
    a11 = c * c * eigs[:, 0] + s * s * eigs[:, 1]
    a12 = c * s * (eigs[:, 1] - eigs[:, 0])
    a22 = s * s * eigs[:, 0] + c * c * eigs[:, 1]
    txx = (a11 - 1.0) * ETAP
    txy = a12 * ETAP
    tyy = (a22 - 1.0) * ETAP
    return txx, txy, tyy


def test_decompose_symmetric_2x2_diagonal_matrix_is_identity_rotation() -> None:
    a11 = np.array([3.0])
    a12 = np.array([0.0])
    a22 = np.array([5.0])

    spec = decompose_symmetric_2x2(a11, a12, a22)

    assert np.allclose(spec.eigval1, 3.0)
    assert np.allclose(spec.eigval2, 5.0)
    assert np.allclose(spec.cos, 1.0)
    assert np.allclose(spec.sin, 0.0)


def test_decompose_symmetric_2x2_matches_numpy_eigh() -> None:
    rng = np.random.default_rng(0)
    n = 200
    a11 = rng.uniform(-2.0, 5.0, n)
    a12 = rng.uniform(-3.0, 3.0, n)
    a22 = rng.uniform(-2.0, 5.0, n)

    spec = decompose_symmetric_2x2(a11, a12, a22)

    for i in range(0, n, 17):
        matrix = np.array([[a11[i], a12[i]], [a12[i], a22[i]]])
        expected = np.linalg.eigvalsh(matrix)
        actual = np.sort([spec.eigval1[i], spec.eigval2[i]])
        assert np.allclose(actual, expected, atol=1.0e-12)

        cos, sin = spec.cos[i], spec.sin[i]
        rotation = np.array([[cos, sin], [-sin, cos]])
        eigenvalues = np.diag([spec.eigval1[i], spec.eigval2[i]])
        reconstructed = rotation @ eigenvalues @ rotation.T
        assert np.allclose(reconstructed, matrix, atol=1.0e-12)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_stress_log_conformation_roundtrip(seed: int) -> None:
    rng = np.random.default_rng(seed)
    txx, txy, tyy = _random_pd_conformation_stress(rng, 500)

    psixx, psixy, psiyy = stress_to_log_conformation(txx, txy, tyy, etap=ETAP)
    txx2, txy2, tyy2 = log_conformation_to_stress(psixx, psixy, psiyy, etap=ETAP)

    assert np.allclose(txx2, txx, atol=1.0e-12)
    assert np.allclose(txy2, txy, atol=1.0e-12)
    assert np.allclose(tyy2, tyy, atol=1.0e-12)


def test_stress_to_log_conformation_zero_stress_is_zero_psi() -> None:
    zero = np.zeros(4)
    psixx, psixy, psiyy = stress_to_log_conformation(zero, zero, zero, etap=ETAP)

    assert np.allclose(psixx, 0.0)
    assert np.allclose(psixy, 0.0)
    assert np.allclose(psiyy, 0.0)


def test_source_term_to_log_conformation_scales_by_inverse_etap() -> None:
    tf_txx = np.array([2.0, -1.0])
    tf_txy = np.array([0.5, 3.0])
    tf_tyy = np.array([-2.0, 4.0])
    etapinv = 1.0 / ETAP

    tf_psixx, tf_psixy, tf_psiyy = source_term_to_log_conformation(
        tf_txx, tf_txy, tf_tyy, etapinv=etapinv
    )

    assert np.allclose(tf_psixx, tf_txx * etapinv)
    assert np.allclose(tf_psixy, tf_txy * etapinv)
    assert np.allclose(tf_psiyy, tf_tyy * etapinv)


def test_kernel_conformation_zero_velocity_gradient_gives_zero_rotation() -> None:
    rng = np.random.default_rng(3)
    txx, txy, tyy = _random_pd_conformation_stress(rng, 50)
    zero = np.zeros_like(txx)

    omegalog, bxx, bxy, byy, mxx, mxy, myy = kernel_conformation(
        txx, txy, tyy, zero, zero, zero, zero, etapinv=1.0 / ETAP, alpha_g=0.0
    )

    assert np.allclose(omegalog, 0.0, atol=1.0e-12)
    assert np.allclose(bxx, 0.0, atol=1.0e-12)
    assert np.allclose(bxy, 0.0, atol=1.0e-12)
    assert np.allclose(byy, 0.0, atol=1.0e-12)
    # Mxx/Mxy/Myy do not depend on the velocity gradient and need not vanish.
    assert np.all(np.isfinite(mxx))
    assert np.all(np.isfinite(mxy))
    assert np.all(np.isfinite(myy))


def test_kernel_conformation_matches_pointwise_linear_algebra() -> None:
    rng = np.random.default_rng(4)
    n = 64
    txx, txy, tyy = _random_pd_conformation_stress(rng, n)
    dudx = rng.uniform(-2.0, 2.0, n)
    dudy = rng.uniform(-2.0, 2.0, n)
    dvdx = rng.uniform(-2.0, 2.0, n)
    dvdy = rng.uniform(-2.0, 2.0, n)
    etapinv = 1.0 / ETAP
    alpha_g = 0.3

    omegalog, bxx, bxy, byy, mxx, mxy, myy = kernel_conformation(
        txx, txy, tyy, dudx, dudy, dvdx, dvdy, etapinv=etapinv, alpha_g=alpha_g
    )

    for i in range(n):
        matrix = np.array(
            [[txx[i] * etapinv + 1.0, txy[i] * etapinv], [txy[i] * etapinv, tyy[i] * etapinv + 1.0]]
        )
        eigvals, eigvecs = np.linalg.eigh(matrix)
        l1, l2 = eigvals
        q = eigvecs
        gradu = np.array([[dudx[i], dudy[i]], [dvdx[i], dvdy[i]]])
        m = q.T @ gradu @ q
        div = l2 - l1
        if abs(div) < 1.0e-14:
            div = 1.0e-14
        omegal = (m[0, 1] * l2 + m[1, 0] * l1) / div

        m_diag = np.diag([m[0, 0], m[1, 1]])
        lam = np.diag([l1, l2])
        jac = np.diag([1.0 / l1, 1.0 / l2])
        b = q @ m_diag @ lam @ jac @ q.T

        m_lambda = np.diag(
            [(1.0 - l1) * (1.0 + alpha_g * (l1 - 1.0)), (1.0 - l2) * (1.0 + alpha_g * (l2 - 1.0))]
        )
        mm = q @ m_lambda @ jac @ q.T

        omega_matrix = q @ np.array([[0.0, omegal], [-omegal, 0.0]]) @ q.T

        assert np.isclose(omegalog[i], omega_matrix[0, 1], atol=1.0e-10)
        assert np.isclose(bxx[i], b[0, 0], atol=1.0e-10)
        assert np.isclose(bxy[i], b[0, 1], atol=1.0e-10)
        assert np.isclose(byy[i], b[1, 1], atol=1.0e-10)
        assert np.isclose(mxx[i], mm[0, 0], atol=1.0e-10)
        assert np.isclose(mxy[i], mm[0, 1], atol=1.0e-10)
        assert np.isclose(myy[i], mm[1, 1], atol=1.0e-10)
