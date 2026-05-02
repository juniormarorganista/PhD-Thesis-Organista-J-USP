"""Fase 4: visualizacao do MMS case 12 + datasets para Machine Learning.

Execute este script na RAIZ do repositorio reorganizado.

Ele NAO altera arquivos Fortran/LaTeX originais. Ele adiciona somente:
- utilitarios de visualizacao Python;
- utilitarios de dataset para ML;
- exemplo executavel;
- testes automatizados;
- documentacao tecnica.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path.cwd()


def write_text(path: str, content: str) -> None:
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content.strip() + "\n", encoding="utf-8", newline="\n")
    print(f"[OK] escrito: {path}")


def append_once(path: str, marker: str, block: str) -> None:
    target = ROOT / path
    if not target.exists():
        target.write_text("", encoding="utf-8", newline="\n")
    text = target.read_text(encoding="utf-8")
    if marker in text:
        print(f"[SKIP] bloco ja existe em: {path}")
        return
    if text and not text.endswith("\n"):
        text += "\n"
    text += "\n" + block.strip() + "\n"
    target.write_text(text, encoding="utf-8", newline="\n")
    print(f"[OK] atualizado: {path}")


def ensure_project_root() -> None:
    required = [
        ROOT / "src" / "mms_viscoelastic",
        ROOT / "src" / "mms_viscoelastic" / "manufactured_solutions" / "case_12.py",
        ROOT / "pyproject.toml",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit(
            "[ERRO] Execute este script na raiz do repositorio apos a Fase 3. "
            f"Itens ausentes: {missing}"
        )


VIS_INIT = r'''
"""Visualization utilities for numerical verification and thesis figures."""

from .fields import plot_case12_panel, plot_scalar_field, save_case12_component_maps

__all__ = [
    "plot_case12_panel",
    "plot_scalar_field",
    "save_case12_component_maps",
]
'''

VIS_FIELDS = r'''
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
'''

DATASETS = r'''
"""Dataset builders for future Machine Learning workflows.

The functions in this module intentionally avoid depending on scikit-learn.
They produce simple, reproducible tabular data that can later feed regression,
operator-learning or physics-informed workflows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.ml import field_to_feature_table


def add_spatial_features(table: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic spatial features useful for baseline ML models.

    Added columns:
    - ``x2``, ``y2`` and ``xy``;
    - low-order trigonometric features using ``pi*x`` and ``pi*y``.
    """
    result = table.copy()
    x = result["x"].to_numpy()
    y = result["y"].to_numpy()
    result["x2"] = x**2
    result["y2"] = y**2
    result["xy"] = x * y
    result["sin_pi_x"] = np.sin(np.pi * x)
    result["cos_pi_x"] = np.cos(np.pi * x)
    result["sin_pi_y"] = np.sin(np.pi * y)
    result["cos_pi_y"] = np.cos(np.pi * y)
    return result


def build_case12_dataset(
    *,
    nx: int = 65,
    ny: int = 65,
    time: float = 0.0,
    beta: float = 0.5,
    at: float = 0.0,
    include_spatial_features: bool = True,
) -> pd.DataFrame:
    """Build a tabular dataset from the exact MMS case-12 fields."""
    grid = Grid2D.from_bounds(nx=nx, ny=ny)
    params = Case12Parameters(beta=beta, at=at)
    field = evaluate_case_12_exact(grid, t=time, parameters=params)
    table = field_to_feature_table(grid, field)
    table.insert(2, "time", float(time))
    table.insert(3, "beta", float(beta))
    if include_spatial_features:
        table = add_spatial_features(table)
    return table


def deterministic_spatial_split(
    table: pd.DataFrame,
    *,
    validation_stride: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a table into train/validation sets using grid-index stride.

    This is deterministic and spatially distributed: every ``validation_stride``
    row goes to validation.  It is intentionally simple and reproducible.
    """
    if validation_stride < 2:
        raise ValueError("validation_stride must be at least 2.")
    index = np.arange(len(table))
    validation_mask = index % validation_stride == 0
    validation = table.loc[validation_mask].reset_index(drop=True)
    train = table.loc[~validation_mask].reset_index(drop=True)
    return train, validation


def save_dataset_bundle(
    table: pd.DataFrame,
    output_dir: str | Path,
    *,
    prefix: str = "case12",
    validation_stride: int = 5,
) -> dict[str, Path]:
    """Save full/train/validation CSV files and a small metadata file."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    train, validation = deterministic_spatial_split(table, validation_stride=validation_stride)

    full_path = output / f"{prefix}_full.csv"
    train_path = output / f"{prefix}_train.csv"
    validation_path = output / f"{prefix}_validation.csv"
    metadata_path = output / f"{prefix}_metadata.txt"

    table.to_csv(full_path, index=False)
    train.to_csv(train_path, index=False)
    validation.to_csv(validation_path, index=False)

    target_columns = ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
    feature_columns = [name for name in table.columns if name not in target_columns]
    metadata = [
        "MMS case 12 dataset bundle",
        f"rows_full={len(table)}",
        f"rows_train={len(train)}",
        f"rows_validation={len(validation)}",
        f"validation_stride={validation_stride}",
        "feature_columns=" + ",".join(feature_columns),
        "target_columns=" + ",".join(target_columns),
    ]
    metadata_path.write_text("\n".join(metadata) + "\n", encoding="utf-8", newline="\n")

    return {
        "full": full_path,
        "train": train_path,
        "validation": validation_path,
        "metadata": metadata_path,
    }
'''

ML_INIT = r'''
"""Machine-learning utilities built on top of verified numerical fields."""

from .datasets import (
    add_spatial_features,
    build_case12_dataset,
    deterministic_spatial_split,
    save_dataset_bundle,
)
from .features import field_to_feature_table

__all__ = [
    "add_spatial_features",
    "build_case12_dataset",
    "deterministic_spatial_split",
    "field_to_feature_table",
    "save_dataset_bundle",
]
'''

EXAMPLE = r'''
"""Generate MMS case-12 figures and ML-ready datasets.

Run from the repository root after installing the package in editable mode:

    python examples/03_case12_visualization_and_ml_dataset.py

Useful options:

    python examples/03_case12_visualization_and_ml_dataset.py --nx 65 --ny 65
    python examples/03_case12_visualization_and_ml_dataset.py --output-dir outputs/case12
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import Case12Parameters, evaluate_case_12_exact
from mms_viscoelastic.ml.datasets import build_case12_dataset, save_dataset_bundle
from mms_viscoelastic.visualization import plot_case12_panel, save_case12_component_maps


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate case-12 figures and ML datasets.")
    parser.add_argument("--nx", type=int, default=65, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=65, help="Number of grid points in y.")
    parser.add_argument("--time", type=float, default=0.0, help="Physical time.")
    parser.add_argument("--beta", type=float, default=0.5, help="Solvent-viscosity ratio.")
    parser.add_argument("--at", type=float, default=0.0, help="Temporal exponential coefficient.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs") / "case12",
        help="Directory where CSV files and figures will be written.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    parameters = Case12Parameters(beta=args.beta, at=args.at)
    field = evaluate_case_12_exact(grid, t=args.time, parameters=parameters)

    figure_dir = args.output_dir / "figures"
    dataset_dir = args.output_dir / "datasets"

    component_paths = save_case12_component_maps(grid, field, figure_dir)
    panel_path = figure_dir / "case12_panel.png"
    panel = plot_case12_panel(grid, field, output_path=panel_path)
    plt.close(panel)

    table = build_case12_dataset(
        nx=args.nx,
        ny=args.ny,
        time=args.time,
        beta=args.beta,
        at=args.at,
        include_spatial_features=True,
    )
    dataset_paths = save_dataset_bundle(table, dataset_dir, prefix="case12", validation_stride=5)

    print("MMS case 12 visualization + ML dataset")
    print(f"grid: nx={grid.nx}, ny={grid.ny}, points={grid.nx * grid.ny}")
    print(f"output directory: {args.output_dir}")
    print("\nFigures:")
    print(f"  {panel_path}")
    for path in component_paths:
        print(f"  {path}")
    print("\nDatasets:")
    for label, path in dataset_paths.items():
        print(f"  {label:>10s}: {path}")
    print("\nDataset columns:")
    print("  " + ", ".join(table.columns))


if __name__ == "__main__":
    main()
'''

TEST_VIS = r'''
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_12_exact
from mms_viscoelastic.visualization import plot_case12_panel, plot_scalar_field, save_case12_component_maps


def test_plot_scalar_field_saves_file(tmp_path: Path) -> None:
    grid = Grid2D.from_bounds(nx=9, ny=7)
    field = evaluate_case_12_exact(grid)
    output = tmp_path / "u.png"

    fig, ax = plot_scalar_field(grid, field.u, title="u", output_path=output)
    plt.close(fig)

    assert output.exists()
    assert output.stat().st_size > 0
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"


def test_plot_scalar_field_rejects_swapped_axes() -> None:
    grid = Grid2D.from_bounds(nx=9, ny=7)
    field = evaluate_case_12_exact(grid)

    with pytest.raises(ValueError, match="shape"):
        plot_scalar_field(grid, field.u.T, title="u transposed")


def test_case12_panel_and_component_maps(tmp_path: Path) -> None:
    grid = Grid2D.from_bounds(nx=11, ny=11)
    field = evaluate_case_12_exact(grid)

    panel_path = tmp_path / "panel.png"
    panel = plot_case12_panel(grid, field, output_path=panel_path)
    plt.close(panel)

    paths = save_case12_component_maps(grid, field, tmp_path / "components", components=["u", "v"])

    assert panel_path.exists()
    assert len(paths) == 2
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)
'''

TEST_DATASETS = r'''
from pathlib import Path

import numpy as np

from mms_viscoelastic.ml.datasets import (
    add_spatial_features,
    build_case12_dataset,
    deterministic_spatial_split,
    save_dataset_bundle,
)


def test_build_case12_dataset_has_expected_columns() -> None:
    table = build_case12_dataset(nx=5, ny=4, include_spatial_features=True)

    expected = {
        "x", "y", "time", "beta",
        "u", "v", "psi", "wz", "txx", "txy", "tyy",
        "x2", "y2", "xy", "sin_pi_x", "cos_pi_x", "sin_pi_y", "cos_pi_y",
    }
    assert table.shape[0] == 20
    assert expected.issubset(set(table.columns))
    assert np.isfinite(table.to_numpy()).all()


def test_add_spatial_features_does_not_mutate_input() -> None:
    table = build_case12_dataset(nx=3, ny=3, include_spatial_features=False)
    original_columns = list(table.columns)
    enriched = add_spatial_features(table)

    assert list(table.columns) == original_columns
    assert "x2" in enriched.columns
    assert "sin_pi_x" in enriched.columns


def test_deterministic_spatial_split_sizes() -> None:
    table = build_case12_dataset(nx=5, ny=5)
    train, validation = deterministic_spatial_split(table, validation_stride=5)

    assert len(table) == 25
    assert len(validation) == 5
    assert len(train) == 20


def test_save_dataset_bundle(tmp_path: Path) -> None:
    table = build_case12_dataset(nx=4, ny=4)
    paths = save_dataset_bundle(table, tmp_path, prefix="case12_test", validation_stride=4)

    assert set(paths) == {"full", "train", "validation", "metadata"}
    assert all(path.exists() and path.stat().st_size > 0 for path in paths.values())
'''

DOC = r'''
# Fase 4 — Visualização do case 12 e datasets para Machine Learning

Esta fase adiciona uma camada prática sobre o port Python do caso manufaturado 12.
O objetivo é transformar a solução exata em dois produtos úteis:

1. figuras de inspeção visual para a tese, apresentação e README;
2. datasets tabulares para futuros experimentos de Machine Learning.

## Arquivos criados

```text
src/mms_viscoelastic/visualization/__init__.py
src/mms_viscoelastic/visualization/fields.py
src/mms_viscoelastic/ml/datasets.py
examples/03_case12_visualization_and_ml_dataset.py
tests/test_case12_visualization.py
tests/test_case12_datasets.py
```

## Convenção de eixos

A função `plot_scalar_field` valida explicitamente que o array tem forma:

```text
(ny, nx)
```

Isso evita a troca acidental de eixos entre a convenção Fortran e a convenção
NumPy/Matplotlib. Esse ponto é crítico para comparar corretamente os resultados
numéricos com as figuras da tese.

## Como gerar figuras e datasets

Na raiz do repositório, após instalar o pacote em modo desenvolvimento:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py
```

Com grade mais fina:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --nx 129 --ny 129
```

Com pasta de saída personalizada:

```powershell
python examples\03_case12_visualization_and_ml_dataset.py --output-dir outputs\case12
```

## Saídas esperadas

```text
outputs/case12/figures/case12_panel.png
outputs/case12/figures/case12_u.png
outputs/case12/figures/case12_v.png
outputs/case12/figures/case12_psi.png
outputs/case12/figures/case12_wz.png
outputs/case12/figures/case12_txx.png
outputs/case12/figures/case12_txy.png
outputs/case12/figures/case12_tyy.png
outputs/case12/datasets/case12_full.csv
outputs/case12/datasets/case12_train.csv
outputs/case12/datasets/case12_validation.csv
outputs/case12/datasets/case12_metadata.txt
```

## Colunas do dataset

As colunas de alvo são:

```text
u, v, psi, wz, txx, txy, tyy
```

As colunas de entrada incluem:

```text
x, y, time, beta, x2, y2, xy, sin_pi_x, cos_pi_x, sin_pi_y, cos_pi_y
```

Esse formato já permite treinar modelos de regressão supervisionada para campos
manufaturados e, futuramente, comparar modelos clássicos de ML, PINNs e operadores
neurais.

## Limitação

Esta fase ainda não usa dados produzidos pelo Fortran. Ela trabalha com a solução
exata Python do `case_12`. A comparação Fortran vs Python deve vir em fase própria,
após a exportação padronizada dos resultados de referência.
'''

GITIGNORE_BLOCK = r'''
# phase-4-generated-outputs
outputs/
*.egg-info/
.pytest_cache/
.mypy_cache/
.ruff_cache/
__pycache__/
'''


def main() -> None:
    ensure_project_root()

    write_text("src/mms_viscoelastic/visualization/__init__.py", VIS_INIT)
    write_text("src/mms_viscoelastic/visualization/fields.py", VIS_FIELDS)
    write_text("src/mms_viscoelastic/ml/datasets.py", DATASETS)
    write_text("src/mms_viscoelastic/ml/__init__.py", ML_INIT)
    write_text("examples/03_case12_visualization_and_ml_dataset.py", EXAMPLE)
    write_text("tests/test_case12_visualization.py", TEST_VIS)
    write_text("tests/test_case12_datasets.py", TEST_DATASETS)
    write_text("docs/case-12-visualization-and-ml-dataset.md", DOC)
    append_once(".gitignore", "# phase-4-generated-outputs", GITIGNORE_BLOCK)

    print("\n[OK] Fase 4 criada. Rode agora:")
    print("  python -m pip install -e .")
    print("  python -m pytest -q")
    print("  python examples\\03_case12_visualization_and_ml_dataset.py")
    print("\nDepois confira:")
    print("  Get-ChildItem outputs\\case12 -Recurse")


if __name__ == "__main__":
    main()
