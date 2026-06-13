"""Plotting helpers for two-dimensional manufactured fields.

The plotting routines preserve the NumPy convention used throughout the
project: arrays have shape ``(ny, nx)``. Axis 0 is ``y`` and axis 1 is ``x``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from mms_viscoelastic.grids import Grid2D


def _validate_component_array(grid: Grid2D, values: object) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    expected = (grid.ny, grid.nx)

    if array.shape != expected:
        raise ValueError(
            f"Expected array with shape {expected}, got {array.shape}. "
            "Remember: NumPy convention is (ny, nx)."
        )

    if not np.isfinite(array).all():
        raise ValueError("Field contains NaN or infinite values.")

    return array


def plot_scalar_field(
    grid: Grid2D,
    values: object,
    *,
    title: str,
    label: str | None = None,
    output_path: str | Path | None = None,
    dpi: int = 180,
) -> tuple[Figure, Axes]:
    """Plot one scalar field on the physical ``x-y`` domain."""

    array = _validate_component_array(grid, values)
    extent = (
        float(grid.x.min()),
        float(grid.x.max()),
        float(grid.y.min()),
        float(grid.y.max()),
    )

    fig, ax = plt.subplots(figsize=(6.0, 4.8), constrained_layout=True)
    image = ax.imshow(
        array,
        origin="lower",
        extent=extent,
        aspect="auto",
        interpolation="nearest",
    )
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(label or title)

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")

    return fig, ax


def plot_case12_panel(
    grid: Grid2D,
    field: object,
    *,
    components: Iterable[str] = ("u", "v", "psi", "wz", "txx", "txy", "tyy"),
    output_path: str | Path | None = None,
    title: str = "MMS case 12 - exact manufactured fields",
    dpi: int = 180,
) -> Figure:
    """Create a compact panel figure for selected Case 12 components."""

    selected = list(components)
    if not selected:
        raise ValueError("At least one component is required.")

    ncols = 3
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4.6 * ncols, 3.8 * nrows),
        squeeze=False,
    )
    axes_array = axes.ravel()

    extent = [
        float(grid.x.min()),
        float(grid.x.max()),
        float(grid.y.min()),
        float(grid.y.max()),
    ]

    for ax, name in zip(axes_array, selected, strict=False):
        values = getattr(field, name, None)
        if values is None:
            raise ValueError(f"Component {name!r} is not available.")
        array = _validate_component_array(grid, values)

        image = ax.imshow(
            array,
            origin="lower",
            extent=extent,
            aspect="auto",
            interpolation="nearest",
        )
        ax.set_title(name)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(image, ax=ax, shrink=0.82)

    for ax in axes_array[len(selected):]:
        ax.axis("off")

    fig.suptitle(title, y=1.02)
    fig.tight_layout()

    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")

    return fig


def save_case12_component_maps(
    grid: Grid2D,
    field: object,
    output_dir: str | Path,
    *,
    components: Iterable[str] = ("u", "v", "psi", "wz", "txx", "txy", "tyy"),
    dpi: int = 180,
) -> list[Path]:
    """Save one image per selected component and return created paths."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []
    for name in components:
        values = getattr(field, name, None)
        if values is None:
            raise ValueError(f"Component {name!r} is not available.")
        path = output / f"case12_{name}.png"
        fig, _ = plot_scalar_field(
            grid,
            values,
            title=f"Case 12 - {name}",
            label=name,
            output_path=path,
            dpi=dpi,
        )
        plt.close(fig)
        created.append(path)

    return created


def plot_field(
    grid: Grid2D,
    values: object,
    *,
    title: str,
    output_path: str | Path | None = None,
    label: str | None = None,
    dpi: int = 180,
) -> tuple[Figure, Axes]:
    """Backward-compatible alias for :func:`plot_scalar_field`."""

    return plot_scalar_field(
        grid,
        values,
        title=title,
        label=label,
        output_path=output_path,
        dpi=dpi,
    )


def plot_field_panel(
    grid: Grid2D,
    field: object,
    *,
    output_path: str | Path | None = None,
    title: str = "MMS case 12 - exact manufactured fields",
    dpi: int = 180,
) -> Figure:
    """Backward-compatible alias for :func:`plot_case12_panel`."""

    return plot_case12_panel(
        grid,
        field,
        output_path=output_path,
        title=title,
        dpi=dpi,
    )


__all__ = [
    "plot_case12_panel",
    "plot_field",
    "plot_field_panel",
    "plot_scalar_field",
    "save_case12_component_maps",
]
