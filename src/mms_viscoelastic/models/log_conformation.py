"""Log-conformation reformulation of the polymeric stress.

Ports ``src/fortran_original/4Order_LPOG_Euler/conform.f``, which is not
discussed in the thesis text (zero occurrences of "conformation" in
``Thesis/tex/``) but is part of the original Fortran source and is ported
here for source-code fidelity (Etapa 15).

The conformation tensor is

    A = I + Txx_p / etap,

a symmetric positive-definite 2x2 matrix built from the polymeric stress
``(Txx, Txy, Tyy)`` (``etap = (1 - beta) / (Rey * Wi)``, see
:mod:`mms_viscoelastic.numerics.time_stepping`). The log-conformation
representation diagonalizes ``A = Q diag(eigval) Q^T`` and replaces it with

    Psi = Q diag(log(eigval)) Q^T,

which removes the loss-of-positive-definiteness failure mode of the raw
stress formulation at high Weissenberg number. ``Psi_to_T`` undoes the
transform via ``A = Q diag(exp(eigval)) Q^T``.

Eigendecomposition
-------------------
``conform.f``'s ``decomp_matrix`` is a generic Jacobi eigenvalue routine
(Numerical Recipes ``jacobi``) for a symmetric ``dimen x dimen`` matrix, but
``dimen = 2`` everywhere it is called. For a 2x2 symmetric matrix the
routine performs exactly one Jacobi rotation before the off-diagonal sum
``sm`` vanishes and it returns, so its output is reproduced here in closed
form rather than via an iterative routine or :func:`numpy.linalg.eigh`
(which would not reproduce the Fortran rotation's eigenvector sign/ordering
convention bit-for-bit). The closed form below follows the same branch
structure as ``decomp_matrix`` (the ``dabs(h)+g.eq.dabs(h)`` floating-point
absorption guard included), so it reduces to the identical rotation
``eigvec = [[c, s], [-s, c]]`` for every input.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mms_viscoelastic.fields import FloatArray


@dataclass(frozen=True)
class SpectralDecomposition2D:
    """Closed-form eigendecomposition of a symmetric 2x2 matrix ``[[a11, a12], [a12, a22]]``.

    ``eigval1``/``eigval2`` and the rotation ``(cos, sin)`` reproduce
    ``decomp_matrix`` (``conform.f``) for ``dimen = 2``: the eigenvectors are
    the columns of ``[[cos, sin], [-sin, cos]]``.
    """

    eigval1: FloatArray
    eigval2: FloatArray
    cos: FloatArray
    sin: FloatArray


def decompose_symmetric_2x2(a11: FloatArray, a12: FloatArray, a22: FloatArray) -> SpectralDecomposition2D:
    """Port of ``decomp_matrix`` (``conform.f``) specialized to ``dimen = 2``.

    Reproduces the single Jacobi rotation that diagonalizes a symmetric 2x2
    matrix, including the ``dabs(h)+g.eq.dabs(h)`` absorption guard used to
    pick the small-rotation branch when ``a12`` is (numerically) zero.
    """

    a11 = np.asarray(a11, dtype=np.float64)
    a12 = np.asarray(a12, dtype=np.float64)
    a22 = np.asarray(a22, dtype=np.float64)

    h = a22 - a11
    g = 100.0 * np.abs(a12)
    negligible = (np.abs(h) + g) == np.abs(h)

    t_small = np.divide(a12, h, out=np.zeros_like(a12), where=(h != 0.0))

    safe_a12 = np.where(a12 == 0.0, 1.0, a12)
    theta = 0.5 * h / safe_a12
    t_large = 1.0 / (np.abs(theta) + np.sqrt(1.0 + theta**2))
    t_large = np.where(theta < 0.0, -t_large, t_large)

    t = np.where(negligible, t_small, t_large)
    cos = 1.0 / np.sqrt(1.0 + t**2)
    sin = t * cos

    eigval1 = a11 - t * a12
    eigval2 = a22 + t * a12

    return SpectralDecomposition2D(eigval1=eigval1, eigval2=eigval2, cos=cos, sin=sin)


def stress_to_log_conformation(
    txx: FloatArray, txy: FloatArray, tyy: FloatArray, *, etap: float
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Port of ``T_to_Psi`` (``conform.f``): polymeric stress -> log-conformation.

    Builds the conformation matrix ``A = I + T_p / etap``, diagonalizes it
    and replaces the eigenvalues with their natural logarithm.
    """

    a11 = txx / etap + 1.0
    a12 = txy / etap
    a22 = tyy / etap + 1.0

    spec = decompose_symmetric_2x2(a11, a12, a22)
    log1 = np.log(spec.eigval1)
    log2 = np.log(spec.eigval2)
    c, s = spec.cos, spec.sin

    # Psi = Q diag(log1, log2) Q^T, Q = [[c, s], [-s, c]].
    psixx = c * c * log1 + s * s * log2
    psixy = c * s * (log2 - log1)
    psiyy = s * s * log1 + c * c * log2

    return psixx, psixy, psiyy


def log_conformation_to_stress(
    psixx: FloatArray, psixy: FloatArray, psiyy: FloatArray, *, etap: float
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Port of ``Psi_to_T`` (``conform.f``): log-conformation -> polymeric stress.

    Diagonalizes ``Psi``, exponentiates the eigenvalues to recover the
    conformation matrix ``A``, then maps back to stress via
    ``T_p = etap * (A - I)``.
    """

    spec = decompose_symmetric_2x2(psixx, psixy, psiyy)
    exp1 = np.exp(spec.eigval1)
    exp2 = np.exp(spec.eigval2)
    c, s = spec.cos, spec.sin

    a11 = c * c * exp1 + s * s * exp2
    a12 = c * s * (exp2 - exp1)
    a22 = s * s * exp1 + c * c * exp2

    txx = etap * (a11 - 1.0)
    txy = etap * a12
    tyy = etap * (a22 - 1.0)

    return txx, txy, tyy


def source_term_to_log_conformation(
    tf_txx: FloatArray, tf_txy: FloatArray, tf_tyy: FloatArray, *, etapinv: float
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Port of ``tf_conformation`` (``conform.f``): scale the MMS source term by ``1/etap``.

    The log-conformation transport equation for ``Psi`` is the ``Psi``-space
    image of the stress transport equation scaled by ``etapinv``, so its
    manufactured source term is simply the stress source term scaled the
    same way (no nonlinear transform is applied to the source term itself).
    """

    tf_psixx = tf_txx * etapinv
    tf_psixy = tf_txy * etapinv
    tf_psiyy = tf_tyy * etapinv
    return tf_psixx, tf_psixy, tf_psiyy


def _matmul2(
    a11: FloatArray,
    a12: FloatArray,
    a21: FloatArray,
    a22: FloatArray,
    b11: FloatArray,
    b12: FloatArray,
    b21: FloatArray,
    b22: FloatArray,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Elementwise product of two fields of 2x2 matrices (port of ``multi_matrix``)."""

    return (
        a11 * b11 + a12 * b21,
        a11 * b12 + a12 * b22,
        a21 * b11 + a22 * b21,
        a21 * b12 + a22 * b22,
    )


def kernel_conformation(
    txx: FloatArray,
    txy: FloatArray,
    tyy: FloatArray,
    dudx: FloatArray,
    dudy: FloatArray,
    dvdx: FloatArray,
    dvdy: FloatArray,
    *,
    etapinv: float,
    alpha_g: float,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Port of ``kernel_conformation`` (``conform.f``).

    Returns ``(omegalog, bxx, bxy, byy, mxx, mxy, myy)``, the rotation rate
    of the eigenbasis (``Omegalog``) and the symmetric tensors ``B`` and
    ``M`` that enter the ``Psi`` transport equation's right-hand side
    (``drv_Psi``, ``drvs.f``). Mirrors the original step-by-step matrix
    composition (rather than an algebraically simplified closed form) to
    keep the port mechanically checkable against the Fortran source.
    """

    a11 = txx * etapinv + 1.0
    a12 = txy * etapinv
    a22 = tyy * etapinv + 1.0
    spec = decompose_symmetric_2x2(a11, a12, a22)
    l1, l2, c, s = spec.eigval1, spec.eigval2, spec.cos, spec.sin

    q11, q12, q21, q22 = c, s, -s, c
    qt11, qt12, qt21, qt22 = c, -s, s, c

    zero = np.zeros_like(l1)
    invl1, invl2 = 1.0 / l1, 1.0 / l2

    # M = Q^T . gradu . Q
    mi11, mi12, mi21, mi22 = _matmul2(qt11, qt12, qt21, qt22, dudx, dudy, dvdx, dvdy)
    m11, m12, m21, m22 = _matmul2(mi11, mi12, mi21, mi22, q11, q12, q21, q22)

    div = l2 - l1
    div = np.where(np.abs(div) < 1e-14, 1e-14, div)
    omegal = (m12 * l2 + m21 * l1) / div

    # B = Q . diag(m11, m22) . diag(l1, l2) . diag(1/l1, 1/l2) . Q^T
    mi11, mi12, mi21, mi22 = _matmul2(q11, q12, q21, q22, m11, zero, zero, m22)
    m11b, m12b, m21b, m22b = _matmul2(mi11, mi12, mi21, mi22, l1, zero, zero, l2)
    mi11, mi12, mi21, mi22 = _matmul2(m11b, m12b, m21b, m22b, invl1, zero, zero, invl2)
    b11, b12, _, b22 = _matmul2(mi11, mi12, mi21, mi22, qt11, qt12, qt21, qt22)
    bxx, bxy, byy = b11, b12, b22

    # M(lambda) = (I - Lambda) . (I + alpha_g * (Lambda - I)), Lambda = diag(l1, l2).
    mlam1 = (1.0 - l1) * (1.0 + alpha_g * (l1 - 1.0))
    mlam2 = (1.0 - l2) * (1.0 + alpha_g * (l2 - 1.0))
    mi11, mi12, mi21, mi22 = _matmul2(q11, q12, q21, q22, mlam1, zero, zero, mlam2)
    m11c, m12c, m21c, m22c = _matmul2(mi11, mi12, mi21, mi22, invl1, zero, zero, invl2)
    mi11, mi12, _, mi22 = _matmul2(m11c, m12c, m21c, m22c, qt11, qt12, qt21, qt22)
    mxx, mxy, myy = mi11, mi12, mi22

    # Omegalog = Q . [[0, omegal], [-omegal, 0]] . Q^T, off-diagonal entry.
    mi11, mi12, mi21, mi22 = _matmul2(q11, q12, q21, q22, zero, omegal, -omegal, zero)
    _, m12d, _, _ = _matmul2(mi11, mi12, mi21, mi22, qt11, qt12, qt21, qt22)
    omegalog = m12d

    return omegalog, bxx, bxy, byy, mxx, mxy, myy


__all__ = [
    "SpectralDecomposition2D",
    "decompose_symmetric_2x2",
    "kernel_conformation",
    "log_conformation_to_stress",
    "source_term_to_log_conformation",
    "stress_to_log_conformation",
]
