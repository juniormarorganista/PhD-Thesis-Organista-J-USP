"""Regression tests for ML package imports."""

from __future__ import annotations


def test_ml_package_exports_field_to_feature_table() -> None:
    from mms_viscoelastic.ml import field_to_feature_table

    assert callable(field_to_feature_table)


def test_dataset_builder_imports_without_circular_error() -> None:
    from mms_viscoelastic.ml.datasets import build_case12_dataset

    assert callable(build_case12_dataset)


def test_visualization_exports_backward_compatible_names() -> None:
    from mms_viscoelastic.visualization.fields import plot_field, plot_field_panel

    assert callable(plot_field)
    assert callable(plot_field_panel)
