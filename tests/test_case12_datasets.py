"""Tests for MMS Case 12 dataset generation."""

from __future__ import annotations

from mms_viscoelastic.ml.datasets import build_case12_dataset


def test_case12_dataset_shapes_and_columns() -> None:
    bundle = build_case12_dataset(nx=9, ny=7, validation_fraction=0.25)

    assert len(bundle.full) == 9 * 7
    assert len(bundle.train) + len(bundle.validation) == len(bundle.full)

    for column in ("x", "y", "time", "beta"):
        assert column in bundle.full.columns

    for column in ("u", "v", "psi", "wz", "txx", "txy", "tyy"):
        assert column in bundle.full.columns
        assert column in bundle.target_columns

    for column in ("x2", "y2", "xy", "sin_pi_x", "cos_pi_x"):
        assert column in bundle.full.columns

    assert "u" not in bundle.feature_columns
    assert "x" in bundle.feature_columns
    assert bundle.metadata["shape_convention"] == "(ny, nx)"
