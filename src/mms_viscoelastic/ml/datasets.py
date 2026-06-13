"""Dataset builders for manufactured-solution machine-learning workflows."""

from __future__ import annotations

import itertools
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.ml.features import TARGET_COLUMNS, field_to_feature_table
from mms_viscoelastic.ml.splits import StratifiedSplit, stratified_train_validation_test_split
from mms_viscoelastic.source_terms import (
    Case12SourceTermParameters,
    evaluate_case_12_source_terms,
)


DEFAULT_TARGET_COLUMNS = TARGET_COLUMNS

SOURCE_TERM_COLUMNS = ("tfwz", "tfTxx", "tfTxy", "tfTyy")

PARAMETRIC_TARGET_COLUMNS = TARGET_COLUMNS + SOURCE_TERM_COLUMNS


@dataclass(frozen=True)
class DatasetBundle:
    """Container for full/train/validation datasets and metadata."""

    full: pd.DataFrame
    train: pd.DataFrame
    validation: pd.DataFrame
    feature_columns: tuple[str, ...]
    target_columns: tuple[str, ...]
    metadata: dict[str, Any]


def infer_feature_columns(
    df: pd.DataFrame,
    *,
    target_columns: tuple[str, ...] = DEFAULT_TARGET_COLUMNS,
    exclude_columns: tuple[str, ...] = (),
) -> tuple[str, ...]:
    """Return ML input columns by excluding target and identifier columns."""

    excluded = set(target_columns) | set(exclude_columns)
    return tuple(column for column in df.columns if column not in excluded)


def split_train_validation(
    df: pd.DataFrame,
    *,
    validation_fraction: float = 0.2,
    random_state: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create a deterministic train/validation split."""

    if not 0.0 < validation_fraction < 1.0:
        raise ValueError("validation_fraction must be between 0 and 1.")

    validation = df.sample(frac=validation_fraction, random_state=random_state)
    train = df.drop(index=validation.index)

    return (
        train.reset_index(drop=True),
        validation.reset_index(drop=True),
    )


def build_case12_dataset(
    *,
    nx: int = 65,
    ny: int = 65,
    time: float = 0.0,
    beta: float = 0.5,
    at: float = 0.0,
    validation_fraction: float = 0.2,
    random_state: int = 42,
    include_engineered_features: bool = True,
) -> DatasetBundle:
    """Build an ML-ready dataset from the MMS Case 12 exact solution."""

    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    parameters = Case12Parameters(beta=beta, at=at)
    fields = evaluate_case_12_exact(grid, t=time, parameters=parameters)

    full = field_to_feature_table(
        grid,
        fields,
        time=time,
        beta=beta,
        include_engineered_features=include_engineered_features,
    )

    target_columns = tuple(
        column for column in DEFAULT_TARGET_COLUMNS if column in full.columns
    )
    feature_columns = infer_feature_columns(full, target_columns=target_columns)

    train, validation = split_train_validation(
        full,
        validation_fraction=validation_fraction,
        random_state=random_state,
    )

    metadata = {
        "case": "case_12",
        "nx": nx,
        "ny": ny,
        "n_points": len(full),
        "time": time,
        "beta": beta,
        "at": at,
        "validation_fraction": validation_fraction,
        "random_state": random_state,
        "feature_columns": list(feature_columns),
        "target_columns": list(target_columns),
        "shape_convention": "(ny, nx)",
    }

    return DatasetBundle(
        full=full,
        train=train,
        validation=validation,
        feature_columns=feature_columns,
        target_columns=target_columns,
        metadata=metadata,
    )


@dataclass(frozen=True)
class ParametricDatasetBundle:
    """Container for a parametric dataset and its train/validation/test splits."""

    full: pd.DataFrame
    train: pd.DataFrame
    validation: pd.DataFrame
    test: pd.DataFrame
    feature_columns: tuple[str, ...]
    target_columns: tuple[str, ...]
    metadata: dict[str, Any]


def build_case12_parametric_dataset(
    *,
    reynolds_values: Sequence[float] = (1.0,),
    wi_values: Sequence[float] = (1.0,),
    beta_values: Sequence[float] = (0.5,),
    time_values: Sequence[float] = (0.0,),
    nx: int = 17,
    ny: int = 17,
    at: float = 0.0,
    train_fraction: float = 0.7,
    validation_fraction: float = 0.15,
    test_fraction: float = 0.15,
    random_state: int = 42,
    include_engineered_features: bool = True,
) -> ParametricDatasetBundle:
    """Build a parametric ML dataset sweeping ``(Re, Wi, beta, t)`` for Case 12.

    For every combination of ``reynolds_values``, ``wi_values``,
    ``beta_values`` and ``time_values``, this evaluates both the exact
    manufactured fields (:func:`evaluate_case_12_exact`) and the manufactured
    source terms (:func:`evaluate_case_12_source_terms`) on a shared grid,
    and stacks the resulting pointwise tables into a single dataframe. Each
    combination is tagged with a ``case_id`` column, used to stratify the
    train/validation/test split so that every split contains points from
    every parameter combination.

    Parameters
    ----------
    reynolds_values, wi_values, beta_values, time_values:
        Parameter values to sweep. The dataset contains one block of grid
        points per element of the Cartesian product of these sequences.
    nx, ny:
        Grid resolution shared by every combination.
    at:
        Time-decay rate shared by every combination (see
        :class:`Case12Parameters`).
    train_fraction, validation_fraction, test_fraction:
        Fractions forwarded to
        :func:`mms_viscoelastic.ml.splits.stratified_train_validation_test_split`.
    random_state:
        Seed used for the stratified split.
    include_engineered_features:
        Whether to add polynomial/trigonometric coordinate features.

    Returns
    -------
    ParametricDatasetBundle
        The full dataset, its stratified train/validation/test splits, and
        metadata describing the sweep.
    """
    grid = Grid2D.from_bounds(nx=nx, ny=ny)

    combinations = list(
        itertools.product(reynolds_values, wi_values, beta_values, time_values)
    )

    tables: list[pd.DataFrame] = []
    for case_id, (reynolds, wi, beta, time) in enumerate(combinations):
        case_parameters = Case12Parameters(beta=beta, at=at)
        field = evaluate_case_12_exact(grid, t=time, parameters=case_parameters)

        source_parameters = Case12SourceTermParameters(
            reynolds=reynolds, wi=wi, beta=beta, at=at
        )
        source_terms = evaluate_case_12_source_terms(grid, t=time, parameters=source_parameters)

        merged = {**asdict(field), **asdict(source_terms)}

        table = field_to_feature_table(
            grid,
            merged,
            time=time,
            beta=beta,
            include_engineered_features=include_engineered_features,
        )
        table["reynolds"] = reynolds
        table["wi"] = wi
        table["case_id"] = case_id
        tables.append(table)

    full = pd.concat(tables, ignore_index=True)

    target_columns = tuple(
        column for column in PARAMETRIC_TARGET_COLUMNS if column in full.columns
    )
    feature_columns = infer_feature_columns(
        full, target_columns=target_columns, exclude_columns=("case_id",)
    )

    split = stratified_train_validation_test_split(
        full,
        stratify_column="case_id",
        train_fraction=train_fraction,
        validation_fraction=validation_fraction,
        test_fraction=test_fraction,
        random_state=random_state,
    )

    metadata = {
        "case": "case_12_parametric",
        "nx": nx,
        "ny": ny,
        "n_points": len(full),
        "n_combinations": len(combinations),
        "reynolds_values": list(reynolds_values),
        "wi_values": list(wi_values),
        "beta_values": list(beta_values),
        "time_values": list(time_values),
        "at": at,
        "train_fraction": train_fraction,
        "validation_fraction": validation_fraction,
        "test_fraction": test_fraction,
        "random_state": random_state,
        "feature_columns": list(feature_columns),
        "target_columns": list(target_columns),
        "shape_convention": "(ny, nx)",
    }

    return ParametricDatasetBundle(
        full=full,
        train=split.train,
        validation=split.validation,
        test=split.test,
        feature_columns=feature_columns,
        target_columns=target_columns,
        metadata=metadata,
    )


def save_dataset_bundle(bundle: DatasetBundle, output_dir: str | Path) -> dict[str, Path]:
    """Save dataset CSV files and metadata text file."""

    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)

    full_path = path / "case12_full.csv"
    train_path = path / "case12_train.csv"
    validation_path = path / "case12_validation.csv"
    metadata_path = path / "case12_metadata.txt"

    bundle.full.to_csv(full_path, index=False)
    bundle.train.to_csv(train_path, index=False)
    bundle.validation.to_csv(validation_path, index=False)

    lines = ["Case 12 ML dataset metadata", ""]
    for key, value in bundle.metadata.items():
        lines.append(f"{key}: {value}")
    lines.append("")
    lines.append(f"feature_columns: {list(bundle.feature_columns)}")
    lines.append(f"target_columns: {list(bundle.target_columns)}")

    metadata_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "full": full_path,
        "train": train_path,
        "validation": validation_path,
        "metadata": metadata_path,
    }


def save_parametric_dataset_bundle(
    bundle: ParametricDatasetBundle, output_dir: str | Path
) -> dict[str, Path]:
    """Save parametric dataset CSV files and metadata text file."""

    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)

    full_path = path / "case12_parametric_full.csv"
    train_path = path / "case12_parametric_train.csv"
    validation_path = path / "case12_parametric_validation.csv"
    test_path = path / "case12_parametric_test.csv"
    metadata_path = path / "case12_parametric_metadata.txt"

    bundle.full.to_csv(full_path, index=False)
    bundle.train.to_csv(train_path, index=False)
    bundle.validation.to_csv(validation_path, index=False)
    bundle.test.to_csv(test_path, index=False)

    lines = ["Case 12 parametric ML dataset metadata", ""]
    for key, value in bundle.metadata.items():
        lines.append(f"{key}: {value}")
    lines.append("")
    lines.append(f"feature_columns: {list(bundle.feature_columns)}")
    lines.append(f"target_columns: {list(bundle.target_columns)}")

    metadata_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "full": full_path,
        "train": train_path,
        "validation": validation_path,
        "test": test_path,
        "metadata": metadata_path,
    }


__all__ = [
    "DEFAULT_TARGET_COLUMNS",
    "PARAMETRIC_TARGET_COLUMNS",
    "SOURCE_TERM_COLUMNS",
    "DatasetBundle",
    "ParametricDatasetBundle",
    "StratifiedSplit",
    "build_case12_dataset",
    "build_case12_parametric_dataset",
    "infer_feature_columns",
    "save_dataset_bundle",
    "save_parametric_dataset_bundle",
    "split_train_validation",
    "stratified_train_validation_test_split",
]
