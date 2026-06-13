"""Stratified train/validation/test splits for ML-ready datasets."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class StratifiedSplit:
    """Train/validation/test partition of a dataframe."""

    train: pd.DataFrame
    validation: pd.DataFrame
    test: pd.DataFrame


def stratified_train_validation_test_split(
    df: pd.DataFrame,
    *,
    stratify_column: str,
    train_fraction: float = 0.7,
    validation_fraction: float = 0.15,
    test_fraction: float = 0.15,
    random_state: int = 42,
) -> StratifiedSplit:
    """Split *df* into train/validation/test, stratified by ``stratify_column``.

    Each unique value of ``stratify_column`` is shuffled and split
    independently using the requested fractions, so every split contains
    rows from every group (e.g. every parameter combination of a parametric
    dataset).

    Parameters
    ----------
    df:
        Dataframe to split.
    stratify_column:
        Column whose unique values define the strata (e.g. ``"case_id"``).
    train_fraction, validation_fraction, test_fraction:
        Fractions of each stratum assigned to each split. Must be positive
        and sum to ``1.0``.
    random_state:
        Seed used to shuffle each stratum, for reproducibility.

    Returns
    -------
    StratifiedSplit
        The train, validation and test dataframes, each with a reset index.

    Raises
    ------
    ValueError
        If ``stratify_column`` is not a column of ``df``, if any fraction is
        not positive, or if the fractions do not sum to ``1.0``.
    """
    if stratify_column not in df.columns:
        raise ValueError(f"Column '{stratify_column}' not found in dataframe.")

    fractions = (train_fraction, validation_fraction, test_fraction)
    if any(fraction <= 0.0 for fraction in fractions):
        raise ValueError(
            "train_fraction, validation_fraction and test_fraction must be positive."
        )
    if not np.isclose(sum(fractions), 1.0):
        raise ValueError(
            "train_fraction + validation_fraction + test_fraction must equal 1.0."
        )

    train_parts: list[pd.DataFrame] = []
    validation_parts: list[pd.DataFrame] = []
    test_parts: list[pd.DataFrame] = []

    for _, group in df.groupby(stratify_column, sort=False):
        shuffled = group.sample(frac=1.0, random_state=random_state)
        n = len(shuffled)
        n_train = int(round(n * train_fraction))
        n_validation = int(round(n * validation_fraction))

        train_parts.append(shuffled.iloc[:n_train])
        validation_parts.append(shuffled.iloc[n_train : n_train + n_validation])
        test_parts.append(shuffled.iloc[n_train + n_validation :])

    return StratifiedSplit(
        train=pd.concat(train_parts).reset_index(drop=True),
        validation=pd.concat(validation_parts).reset_index(drop=True),
        test=pd.concat(test_parts).reset_index(drop=True),
    )


__all__ = ["StratifiedSplit", "stratified_train_validation_test_split"]
