"""Feature engineering helpers for future ML workflows."""

from __future__ import annotations

import pandas as pd

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D


def field_to_feature_table(grid: Grid2D, field: Field2D) -> pd.DataFrame:
    """Convert grid coordinates and available field components into a DataFrame."""
    data = {
        "x": grid.X.reshape(-1),
        "y": grid.Y.reshape(-1),
    }
    for component in field.available_components():
        value = getattr(field, component)
        if value is not None:
            data[component] = value.reshape(-1)
    return pd.DataFrame(data)
