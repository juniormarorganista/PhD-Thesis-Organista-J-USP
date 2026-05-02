"""Structured Cartesian grids."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Grid2D:
    """Uniform two-dimensional Cartesian grid.

    Arrays follow NumPy's matrix convention: field values have shape ``(ny, nx)``.
    Therefore, axis 0 is the y-direction and axis 1 is the x-direction.
    """

    x: FloatArray
    y: FloatArray

    @classmethod
    def from_bounds(
        cls,
        nx: int,
        ny: int,
        xlim: tuple[float, float] = (0.0, 1.0),
        ylim: tuple[float, float] = (0.0, 1.0),
    ) -> "Grid2D":
        if nx < 2 or ny < 2:
            raise ValueError("nx and ny must be at least 2.")
        x = np.linspace(xlim[0], xlim[1], nx, dtype=np.float64)
        y = np.linspace(ylim[0], ylim[1], ny, dtype=np.float64)
        return cls(x=x, y=y)

    @property
    def nx(self) -> int:
        return int(self.x.size)

    @property
    def ny(self) -> int:
        return int(self.y.size)

    @property
    def dx(self) -> float:
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        return float(self.y[1] - self.y[0])

    @property
    def X(self) -> FloatArray:
        return np.meshgrid(self.x, self.y, indexing="xy")[0]

    @property
    def Y(self) -> FloatArray:
        return np.meshgrid(self.x, self.y, indexing="xy")[1]
