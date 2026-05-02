from pathlib import Path

import numpy as np

from mms_viscoelastic.ml.datasets import (
    add_spatial_features,
    build_case12_dataset,
    deterministic_spatial_split,
    save_dataset_bundle,
)


def test_build_case12_dataset_has_expected_columns() -> None:
    table = build_case12_dataset(nx=5, ny=4, include_spatial_features=True)

    expected = {
        "x", "y", "time", "beta",
        "u", "v", "psi", "wz", "txx", "txy", "tyy",
        "x2", "y2", "xy", "sin_pi_x", "cos_pi_x", "sin_pi_y", "cos_pi_y",
    }
    assert table.shape[0] == 20
    assert expected.issubset(set(table.columns))
    assert np.isfinite(table.to_numpy()).all()


def test_add_spatial_features_does_not_mutate_input() -> None:
    table = build_case12_dataset(nx=3, ny=3, include_spatial_features=False)
    original_columns = list(table.columns)
    enriched = add_spatial_features(table)

    assert list(table.columns) == original_columns
    assert "x2" in enriched.columns
    assert "sin_pi_x" in enriched.columns


def test_deterministic_spatial_split_sizes() -> None:
    table = build_case12_dataset(nx=5, ny=5)
    train, validation = deterministic_spatial_split(table, validation_stride=5)

    assert len(table) == 25
    assert len(validation) == 5
    assert len(train) == 20


def test_save_dataset_bundle(tmp_path: Path) -> None:
    table = build_case12_dataset(nx=4, ny=4)
    paths = save_dataset_bundle(table, tmp_path, prefix="case12_test", validation_stride=4)

    assert set(paths) == {"full", "train", "validation", "metadata"}
    assert all(path.exists() and path.stat().st_size > 0 for path in paths.values())
