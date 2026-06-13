"""Tests for the parametric MMS Case 12 ML dataset builder."""

from __future__ import annotations

from mms_viscoelastic.ml.datasets import (
    PARAMETRIC_TARGET_COLUMNS,
    SOURCE_TERM_COLUMNS,
    build_case12_parametric_dataset,
)


def test_parametric_dataset_shapes_and_columns() -> None:
    bundle = build_case12_parametric_dataset(
        reynolds_values=(1.0, 10.0),
        wi_values=(1.0,),
        beta_values=(0.5,),
        time_values=(0.0, 0.5),
        nx=5,
        ny=5,
    )

    n_combinations = 2 * 1 * 1 * 2
    n_points_per_combination = 5 * 5

    assert bundle.metadata["n_combinations"] == n_combinations
    assert len(bundle.full) == n_combinations * n_points_per_combination

    for column in ("x", "y", "time", "beta", "reynolds", "wi", "case_id"):
        assert column in bundle.full.columns

    for column in PARAMETRIC_TARGET_COLUMNS:
        assert column in bundle.full.columns
        assert column in bundle.target_columns

    for column in SOURCE_TERM_COLUMNS:
        assert column in bundle.full.columns

    assert "case_id" not in bundle.feature_columns
    assert "u" not in bundle.feature_columns
    assert "x" in bundle.feature_columns
    assert "reynolds" in bundle.feature_columns
    assert "wi" in bundle.feature_columns


def test_parametric_dataset_case_ids_match_combinations() -> None:
    bundle = build_case12_parametric_dataset(
        reynolds_values=(1.0, 5.0),
        wi_values=(0.5, 1.0),
        beta_values=(0.5,),
        time_values=(0.0,),
        nx=4,
        ny=4,
    )

    n_combinations = 2 * 2 * 1 * 1
    assert sorted(bundle.full["case_id"].unique()) == list(range(n_combinations))


def test_parametric_dataset_split_covers_every_combination() -> None:
    bundle = build_case12_parametric_dataset(
        reynolds_values=(1.0, 10.0),
        wi_values=(1.0,),
        beta_values=(0.3, 0.7),
        time_values=(0.0,),
        nx=5,
        ny=5,
    )

    total = len(bundle.train) + len(bundle.validation) + len(bundle.test)
    assert total == len(bundle.full)

    for case_id in bundle.full["case_id"].unique():
        assert (bundle.train["case_id"] == case_id).any()
        assert (bundle.validation["case_id"] == case_id).any()
        assert (bundle.test["case_id"] == case_id).any()


def test_parametric_dataset_default_single_combination() -> None:
    bundle = build_case12_parametric_dataset(nx=5, ny=5)

    assert bundle.metadata["n_combinations"] == 1
    assert len(bundle.full) == 5 * 5
