"""Machine-learning utilities built on top of verified numerical fields."""

from .features import field_to_feature_table
from .datasets import (
    add_spatial_features,
    build_case12_dataset,
    deterministic_spatial_split,
    save_dataset_bundle,
)

__all__ = [
    "add_spatial_features",
    "build_case12_dataset",
    "deterministic_spatial_split",
    "field_to_feature_table",
    "save_dataset_bundle",
]
