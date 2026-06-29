"""Machine-learning utilities for :mod:`mms_viscoelastic`.

This subpackage is an **exploratory extension**, not part of the thesis's
core scope. The thesis is a numerical-verification work (the Method of
Manufactured Solutions applied to high-order viscoelastic flow solvers);
it is not about machine learning. ``ml/`` only builds pointwise feature
tables and train/validation/test splits from the already-validated
manufactured fields and source terms, for optional future surrogate-model
experiments. No surrogate models are trained anywhere in this repository.
See ``README.md`` and ``.claude/PLANEJAMENTO.md`` (Etapa 11, descontinuada)
for the scope decision.
"""

from mms_viscoelastic.ml.features import (
    TARGET_COLUMNS,
    add_engineered_features,
    field_mapping,
    field_to_feature_table,
    grid_coordinates,
)
from mms_viscoelastic.ml.datasets import (
    DEFAULT_TARGET_COLUMNS,
    PARAMETRIC_TARGET_COLUMNS,
    SOURCE_TERM_COLUMNS,
    DatasetBundle,
    ParametricDatasetBundle,
    build_case12_dataset,
    build_case12_parametric_dataset,
    infer_feature_columns,
    save_dataset_bundle,
    save_parametric_dataset_bundle,
    split_train_validation,
)
from mms_viscoelastic.ml.splits import (
    StratifiedSplit,
    stratified_train_validation_test_split,
)

__all__ = [
    "DEFAULT_TARGET_COLUMNS",
    "PARAMETRIC_TARGET_COLUMNS",
    "SOURCE_TERM_COLUMNS",
    "DatasetBundle",
    "ParametricDatasetBundle",
    "StratifiedSplit",
    "TARGET_COLUMNS",
    "add_engineered_features",
    "build_case12_dataset",
    "build_case12_parametric_dataset",
    "field_mapping",
    "field_to_feature_table",
    "grid_coordinates",
    "infer_feature_columns",
    "save_dataset_bundle",
    "save_parametric_dataset_bundle",
    "split_train_validation",
    "stratified_train_validation_test_split",
]
