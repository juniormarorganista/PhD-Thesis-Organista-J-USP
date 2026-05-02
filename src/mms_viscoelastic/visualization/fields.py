"""Plotting helpers for 2D manufactured fields.

The plotting routines preserve the project-wide NumPy convention:
arrays have shape ``(ny, nx)``, axis 0 is ``y`` and axis 1 is ``x``.
This avoids the common thesis-project bug: beautiful plots with swapped axes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib

# Safe backend for scripts, tests and GitHub Actions without a graphical display.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import NDArray

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D

FloatArray = NDArray[np.float64]


def _validate_component_array(grid: Grid2D, values: FloatArray) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
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
    values: FloatArray,
    *,
    title: str,
    label: str | None = None,
    output_path: str | Path | None = None,
    dpi: int = 180,
) -> tuple[Figure, Axes]:
    """Plot one scalar field on the physical x-y domain.

    Parameters
    ----------
    grid:
        Structured Cartesian grid.
    values:
        Two-dimensional array with shape ``(grid.ny, grid.nx)``.
    title:
        Figure title.
    label:
        Colorbar label.  If omitted, ``title`` is used.
    output_path:
        Optional PNG/PDF/SVG path.  Parent directories are created automatically.
    dpi:
        Resolution used when saving the figure.

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        Created figure and axes.  The caller may further customize or close them.
    """
    array = _validate_component_array(grid, values)
    extent = [float(grid.x.min()), float(grid.x.max()), float(grid.y.min()), float(grid.y.max())]

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
    field: Field2D,
    *,
    components: Iterable[str] = ("u", "v", "psi", "wz", "txx", "txy", "tyy"),
    output_path: str | Path | None = None,
    dpi: int = 180,
) -> Figure:
    """Create a compact panel figure for selected case-12 components."""
    selected = list(components)
    if not selected:
        raise ValueError("At least one component is required.")

    ncols = 3
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.6 * ncols, 3.8 * nrows))
    axes_array = np.atleast_1d(axes).reshape(-1)
    extent = [float(grid.x.min()), float(grid.x.max()), float(grid.y.min()), float(grid.y.max())]

    for ax, name in zip(axes_array, selected, strict=False):
        values = getattr(field, name, None)
        if values is None:
            raise ValueError(f"Component {name!r} is not available in Field2D.")
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

    fig.suptitle("MMS case 12 — exact manufactured fields", y=1.02)
    fig.tight_layout()

    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=dpi, bbox_inches="tight")

    return fig


def save_case12_component_maps(
    grid: Grid2D,
    field: Field2D,
    output_dir: str | Path,
    *,
    components: Iterable[str] = ("u", "v", "psi", "wz", "txx", "txy", "tyy"),
    dpi: int = 180,
) -> list[Path]:
    """Save one image per selected component and return the created paths."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []
    for name in components:
        values = getattr(field, name, None)
        if values is None:
            raise ValueError(f"Component {name!r} is not available in Field2D.")
        path = output / f"case12_{name}.png"
        fig, _ = plot_scalar_field(grid, values, title=f"case 12: {name}", label=name, output_path=path, dpi=dpi)
        plt.close(fig)
        created.append(path)
    return created
