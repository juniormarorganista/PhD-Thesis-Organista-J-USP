"""Containers for two-dimensional flow fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass
class Field2D:
    """Flow-field container prepared for numerical verification and ML datasets."""

    u: FloatArray | None = None
    v: FloatArray | None = None
    psi: FloatArray | None = None
    wz: FloatArray | None = None
    txx: FloatArray | None = None
    txy: FloatArray | None = None
    tyy: FloatArray | None = None

    def available_components(self) -> list[str]:
        """Return names of components currently stored in the object."""
        names = ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
        return [name for name in names if getattr(self, name) is not None]

    def as_feature_matrix(self, components: list[str] | None = None) -> FloatArray:
        """Flatten selected components into a matrix suitable for ML workflows.

        Returns
        -------
        numpy.ndarray
            Matrix with shape ``(n_points, n_components)``.
        """
        selected = components or self.available_components()
        if not selected:
            raise ValueError("No components available to build a feature matrix.")

        columns: list[FloatArray] = []
        expected_shape: tuple[int, ...] | None = None
        for name in selected:
            value = getattr(self, name)
            if value is None:
                raise ValueError(f"Component {name!r} is not available.")
            array = np.asarray(value, dtype=np.float64)
            if expected_shape is None:
                expected_shape = array.shape
            elif array.shape != expected_shape:
                raise ValueError("All selected components must have the same shape.")
            columns.append(array.reshape(-1))
        return np.column_stack(columns)
