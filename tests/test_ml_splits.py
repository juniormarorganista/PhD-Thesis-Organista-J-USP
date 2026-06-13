"""Tests for stratified train/validation/test splitting."""

from __future__ import annotations

import pandas as pd
import pytest

from mms_viscoelastic.ml.splits import stratified_train_validation_test_split


def _sample_dataframe(n_groups: int = 3, n_per_group: int = 20) -> pd.DataFrame:
    rows = [
        {"case_id": group, "value": index}
        for group in range(n_groups)
        for index in range(n_per_group)
    ]
    return pd.DataFrame(rows)


def test_missing_stratify_column_raises() -> None:
    df = _sample_dataframe()

    with pytest.raises(ValueError):
        stratified_train_validation_test_split(df, stratify_column="missing")


@pytest.mark.parametrize(
    "fractions",
    [
        (0.0, 0.5, 0.5),
        (0.5, -0.1, 0.6),
        (0.5, 0.5, 0.5),
        (0.3, 0.3, 0.3),
    ],
)
def test_invalid_fractions_raise(fractions: tuple[float, float, float]) -> None:
    df = _sample_dataframe()
    train_fraction, validation_fraction, test_fraction = fractions

    with pytest.raises(ValueError):
        stratified_train_validation_test_split(
            df,
            stratify_column="case_id",
            train_fraction=train_fraction,
            validation_fraction=validation_fraction,
            test_fraction=test_fraction,
        )


def test_split_covers_every_row_without_overlap() -> None:
    df = _sample_dataframe()

    split = stratified_train_validation_test_split(df, stratify_column="case_id")

    total = len(split.train) + len(split.validation) + len(split.test)
    assert total == len(df)

    train_index = set(zip(split.train["case_id"], split.train["value"]))
    validation_index = set(zip(split.validation["case_id"], split.validation["value"]))
    test_index = set(zip(split.test["case_id"], split.test["value"]))

    assert train_index.isdisjoint(validation_index)
    assert train_index.isdisjoint(test_index)
    assert validation_index.isdisjoint(test_index)


def test_split_is_stratified_per_group() -> None:
    df = _sample_dataframe(n_groups=3, n_per_group=20)

    split = stratified_train_validation_test_split(
        df,
        stratify_column="case_id",
        train_fraction=0.7,
        validation_fraction=0.15,
        test_fraction=0.15,
    )

    for case_id in df["case_id"].unique():
        n_train = (split.train["case_id"] == case_id).sum()
        n_validation = (split.validation["case_id"] == case_id).sum()
        n_test = (split.test["case_id"] == case_id).sum()

        assert n_train == 14
        assert n_validation == 3
        assert n_test == 3


def test_split_is_reproducible_with_fixed_random_state() -> None:
    df = _sample_dataframe()

    first = stratified_train_validation_test_split(df, stratify_column="case_id", random_state=7)
    second = stratified_train_validation_test_split(df, stratify_column="case_id", random_state=7)

    pd.testing.assert_frame_equal(first.train, second.train)
    pd.testing.assert_frame_equal(first.validation, second.validation)
    pd.testing.assert_frame_equal(first.test, second.test)
