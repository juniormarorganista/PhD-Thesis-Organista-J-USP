# Reorganizacao inicial do repositorio da tese + fundacao da biblioteca Python
# Execute este script na RAIZ do repositorio clonado.
# Ele usa git mv para mover diretorios existentes sem alterar o conteudo dos arquivos originais.

$ErrorActionPreference = "Stop"

function Write-Info($msg) {
    Write-Host "[INFO] $msg" -ForegroundColor Cyan
}

function Write-Warn($msg) {
    Write-Host "[AVISO] $msg" -ForegroundColor Yellow
}

function Ensure-Dir($path) {
    if (-not (Test-Path $path)) {
        New-Item -ItemType Directory -Force -Path $path | Out-Null
    }
}

function Write-NewFile($path, $content) {
    $parent = Split-Path $path -Parent
    if ($parent) { Ensure-Dir $parent }
    if (Test-Path $path) {
        Write-Warn "Arquivo ja existe; nao vou sobrescrever: $path"
    }
    else {
        $content | Set-Content -Path $path -Encoding UTF8
        Write-Info "Criado: $path"
    }
}

Write-Info "Verificando se estamos dentro de um repositorio Git..."
git rev-parse --is-inside-work-tree | Out-Null

Write-Info "Configurando long paths no Git local deste repositorio..."
git config core.longpaths true

Write-Info "Criando branch de trabalho..."
$currentBranch = git branch --show-current
if ($currentBranch -ne "chore/reorganize-and-python-foundation") {
    git checkout -b chore/reorganize-and-python-foundation
}
else {
    Write-Info "Ja estamos na branch chore/reorganize-and-python-foundation"
}

Write-Info "Criando diretorios-base..."
Ensure-Dir "src"
Ensure-Dir "thesis"
Ensure-Dir "papers"
Ensure-Dir "docs"
Ensure-Dir "examples"
Ensure-Dir "tests"

Write-Info "Movendo diretorios originais com git mv, sem alterar conteudo..."
if ((Test-Path "Fortran_Code") -and (-not (Test-Path "src/fortran_original"))) {
    git mv "Fortran_Code" "src/fortran_original"
}
else { Write-Warn "Fortran_Code nao encontrado ou src/fortran_original ja existe; movimento ignorado." }

if ((Test-Path "Mathematica_Code") -and (-not (Test-Path "src/mathematica_original"))) {
    git mv "Mathematica_Code" "src/mathematica_original"
}
else { Write-Warn "Mathematica_Code nao encontrado ou src/mathematica_original ja existe; movimento ignorado." }

if ((Test-Path "Papers") -and (-not (Test-Path "papers/original"))) {
    git mv "Papers" "papers/original"
}
else { Write-Warn "Papers nao encontrado ou papers/original ja existe; movimento ignorado." }

if ((Test-Path "Thesis") -and (-not (Test-Path "thesis/original"))) {
    git mv "Thesis" "thesis/original"
}
else { Write-Warn "Thesis nao encontrado ou thesis/original ja existe; movimento ignorado." }

Write-Info "Criando fundacao da biblioteca Python mms_viscoelastic..."
Ensure-Dir "src/mms_viscoelastic"
Ensure-Dir "src/mms_viscoelastic/manufactured_solutions"
Ensure-Dir "src/mms_viscoelastic/numerics"
Ensure-Dir "src/mms_viscoelastic/models"
Ensure-Dir "src/mms_viscoelastic/ml"
Ensure-Dir "src/mms_viscoelastic/io"

Write-NewFile "pyproject.toml" @'
[build-system]
requires = ["setuptools>=69", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "mms-viscoelastic"
version = "0.1.0"
description = "Python foundation for manufactured-solution verification of incompressible viscoelastic flow solvers."
readme = "README.md"
requires-python = ">=3.10"
authors = [
    { name = "Juniormar Organista" }
]
keywords = [
    "method-of-manufactured-solutions",
    "viscoelastic-flows",
    "computational-fluid-dynamics",
    "finite-differences",
    "scientific-machine-learning"
]
classifiers = [
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Scientific/Engineering :: Physics"
]
dependencies = [
    "numpy>=1.26",
    "scipy>=1.11",
    "pandas>=2.0",
    "matplotlib>=3.8"
]

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "ruff>=0.6",
    "mypy>=1.8"
]
ml = [
    "scikit-learn>=1.4",
    "joblib>=1.3"
]
notebooks = [
    "jupyterlab>=4.0",
    "ipykernel>=6.0"
]

[tool.setuptools]
package-dir = {"" = "src"}

[tool.setuptools.packages.find]
where = ["src"]
include = ["mms_viscoelastic*"]

[tool.pytest.ini_options]
pythonpath = ["src"]
testpaths = ["tests"]
addopts = "-q"

[tool.ruff]
line-length = 100
target-version = "py310"

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "B"]
ignore = []
'@

Write-NewFile "requirements.txt" @'
numpy>=1.26
scipy>=1.11
pandas>=2.0
matplotlib>=3.8
'@

Write-NewFile "requirements-dev.txt" @'
-r requirements.txt
pytest>=8.0
ruff>=0.6
mypy>=1.8
scikit-learn>=1.4
joblib>=1.3
jupyterlab>=4.0
ipykernel>=6.0
'@

Write-NewFile "src/mms_viscoelastic/__init__.py" @'
"""Python foundation for MMS-based verification of viscoelastic flow solvers.

The package starts as a faithful, testable Python porting layer around the
original Fortran thesis code.  The original Fortran files are preserved in the
repository and should be treated as the reference implementation during the
porting process.
"""

from ._version import __version__
from .fields import Field2D
from .grids import Grid2D
from .parameters import FlowParameters

__all__ = [
    "__version__",
    "Field2D",
    "FlowParameters",
    "Grid2D",
]
'@

Write-NewFile "src/mms_viscoelastic/_version.py" @'
"""Package version."""

__version__ = "0.1.0"
'@

Write-NewFile "src/mms_viscoelastic/parameters.py" @'
"""Physical and numerical parameters used by the Python port."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FlowParameters:
    """Dimensionless parameters for incompressible viscoelastic simulations.

    Parameters
    ----------
    re:
        Reynolds number.
    wi:
        Weissenberg number.
    beta:
        Solvent-viscosity ratio.  In many formulations, 0 < beta <= 1.
    epsilon:
        PTT-type nonlinear parameter.  Keep zero for UCM/Oldroyd-B/Giesekus cases.
    xi:
        Slip/asymmetry parameter used by some PTT-type formulations.
    alpha_g:
        Giesekus mobility parameter.
    dt:
        Time step.
    t_final:
        Final simulation time.
    """

    re: float = 1.0
    wi: float = 1.0
    beta: float = 0.5
    epsilon: float = 0.0
    xi: float = 0.0
    alpha_g: float = 0.0
    dt: float = 1.0e-6
    t_final: float = 0.05

    def validate(self) -> None:
        """Raise ``ValueError`` if parameters are outside basic admissible ranges."""
        if self.re <= 0:
            raise ValueError("re must be positive.")
        if self.wi <= 0:
            raise ValueError("wi must be positive.")
        if not (0 < self.beta <= 1):
            raise ValueError("beta must satisfy 0 < beta <= 1.")
        if self.dt <= 0:
            raise ValueError("dt must be positive.")
        if self.t_final <= 0:
            raise ValueError("t_final must be positive.")
        if self.alpha_g < 0:
            raise ValueError("alpha_g must be non-negative.")
        if self.epsilon < 0:
            raise ValueError("epsilon must be non-negative.")
'@

Write-NewFile "src/mms_viscoelastic/grids.py" @'
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
'@

Write-NewFile "src/mms_viscoelastic/fields.py" @'
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
'@

Write-NewFile "src/mms_viscoelastic/numerics/__init__.py" @'
"""Numerical operators and verification utilities."""

from .derivatives import derivative_x, derivative_y, laplacian
from .errors import estimate_order, l2_error, linf_error, relative_l2_error

__all__ = [
    "derivative_x",
    "derivative_y",
    "estimate_order",
    "l2_error",
    "laplacian",
    "linf_error",
    "relative_l2_error",
]
'@

Write-NewFile "src/mms_viscoelastic/numerics/derivatives.py" @'
"""Finite-difference derivative operators.

This module starts with robust second-order operators.  The compact high-order
operators from the original Fortran implementation should be ported here in a
later phase and validated against manufactured solutions.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def _as_float_array(values: FloatArray) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError("Expected a two-dimensional array with shape (ny, nx).")
    return array


def derivative_x(values: FloatArray, dx: float) -> FloatArray:
    """Return the numerical derivative with respect to x.

    Axis convention: x corresponds to axis 1.
    """
    if dx <= 0:
        raise ValueError("dx must be positive.")
    f = _as_float_array(values)
    out = np.empty_like(f)
    out[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dx)
    out[:, 0] = (-3.0 * f[:, 0] + 4.0 * f[:, 1] - f[:, 2]) / (2.0 * dx)
    out[:, -1] = (3.0 * f[:, -1] - 4.0 * f[:, -2] + f[:, -3]) / (2.0 * dx)
    return out


def derivative_y(values: FloatArray, dy: float) -> FloatArray:
    """Return the numerical derivative with respect to y.

    Axis convention: y corresponds to axis 0.
    """
    if dy <= 0:
        raise ValueError("dy must be positive.")
    f = _as_float_array(values)
    out = np.empty_like(f)
    out[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * dy)
    out[0, :] = (-3.0 * f[0, :] + 4.0 * f[1, :] - f[2, :]) / (2.0 * dy)
    out[-1, :] = (3.0 * f[-1, :] - 4.0 * f[-2, :] + f[-3, :]) / (2.0 * dy)
    return out


def laplacian(values: FloatArray, dx: float, dy: float) -> FloatArray:
    """Return the second-order finite-difference Laplacian."""
    if dx <= 0 or dy <= 0:
        raise ValueError("dx and dy must be positive.")
    f = _as_float_array(values)
    out = np.zeros_like(f)
    out[1:-1, 1:-1] = (
        (f[1:-1, 2:] - 2.0 * f[1:-1, 1:-1] + f[1:-1, :-2]) / dx**2
        + (f[2:, 1:-1] - 2.0 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / dy**2
    )
    return out
'@

Write-NewFile "src/mms_viscoelastic/numerics/errors.py" @'
"""Error norms and convergence-rate utilities."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


def l2_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the discrete L2 error between numerical and exact arrays."""
    numerical = np.asarray(numerical, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    if numerical.shape != exact.shape:
        raise ValueError("numerical and exact arrays must have the same shape.")
    return float(np.sqrt(np.mean((numerical - exact) ** 2)))


def relative_l2_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the relative discrete L2 error."""
    denominator = l2_error(np.zeros_like(exact, dtype=np.float64), exact)
    if denominator == 0:
        raise ValueError("Cannot compute relative error because exact field has zero L2 norm.")
    return l2_error(numerical, exact) / denominator


def linf_error(numerical: FloatArray, exact: FloatArray) -> float:
    """Return the maximum absolute error."""
    numerical = np.asarray(numerical, dtype=np.float64)
    exact = np.asarray(exact, dtype=np.float64)
    if numerical.shape != exact.shape:
        raise ValueError("numerical and exact arrays must have the same shape.")
    return float(np.max(np.abs(numerical - exact)))


def estimate_order(h_values: FloatArray, errors: FloatArray) -> float:
    """Estimate convergence order from mesh sizes and errors using log-log fit."""
    h = np.asarray(h_values, dtype=np.float64)
    err = np.asarray(errors, dtype=np.float64)
    if h.ndim != 1 or err.ndim != 1 or h.size != err.size:
        raise ValueError("h_values and errors must be one-dimensional arrays with same size.")
    if h.size < 2:
        raise ValueError("At least two mesh levels are required.")
    if np.any(h <= 0) or np.any(err <= 0):
        raise ValueError("Mesh sizes and errors must be positive.")
    slope, _ = np.polyfit(np.log(h), np.log(err), deg=1)
    return float(slope)
'@

Write-NewFile "src/mms_viscoelastic/manufactured_solutions/__init__.py" @'
"""Manufactured solutions for verification studies."""

from .toy import evaluate_toy_solution

__all__ = ["evaluate_toy_solution"]
'@

Write-NewFile "src/mms_viscoelastic/manufactured_solutions/toy.py" @'
"""Small analytical manufactured field used for smoke tests.

This is not yet one of the thesis Fortran cases.  It exists only to validate the
Python package structure, grid conventions and derivative/error utilities before
the faithful port of ``case_01.f`` ... ``case_14.f``.
"""

from __future__ import annotations

import numpy as np

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D


def evaluate_toy_solution(grid: Grid2D) -> Field2D:
    """Evaluate a smooth divergence-free toy field on ``grid``."""
    x = grid.X
    y = grid.Y
    psi = np.sin(np.pi * x) * np.sin(np.pi * y)
    u = np.pi * np.sin(np.pi * x) * np.cos(np.pi * y)
    v = -np.pi * np.cos(np.pi * x) * np.sin(np.pi * y)
    wz = 2.0 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    txx = 1.0 + 0.1 * np.sin(np.pi * x)
    txy = 0.1 * np.sin(np.pi * x) * np.sin(np.pi * y)
    tyy = 1.0 + 0.1 * np.cos(np.pi * y)
    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)
'@

Write-NewFile "src/mms_viscoelastic/models/__init__.py" @'
"""Constitutive-model helpers."""

from .constitutive import ConstitutiveModel

__all__ = ["ConstitutiveModel"]
'@

Write-NewFile "src/mms_viscoelastic/models/constitutive.py" @'
"""Constitutive model identifiers used during the Fortran-to-Python port."""

from __future__ import annotations

from enum import StrEnum


class ConstitutiveModel(StrEnum):
    """Supported model labels."""

    NEWTONIAN = "newtonian"
    UCM = "ucm"
    OLDROYD_B = "oldroyd_b"
    GIESEKUS = "giesekus"
    LPTT = "lptt"
'@

Write-NewFile "src/mms_viscoelastic/ml/__init__.py" @'
"""Machine-learning utilities built on top of verified numerical fields."""

from .features import field_to_feature_table

__all__ = ["field_to_feature_table"]
'@

Write-NewFile "src/mms_viscoelastic/ml/features.py" @'
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
'@

Write-NewFile "src/mms_viscoelastic/io/__init__.py" @'
"""Input/output helpers."""

from .reference import read_reference_csv

__all__ = ["read_reference_csv"]
'@

Write-NewFile "src/mms_viscoelastic/io/reference.py" @'
"""Reference-data loading utilities."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_reference_csv(path: str | Path) -> pd.DataFrame:
    """Read a CSV file containing reference numerical data."""
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    return pd.read_csv(csv_path)
'@

Write-NewFile "tests/test_import.py" @'
def test_package_imports() -> None:
    import mms_viscoelastic

    assert mms_viscoelastic.__version__ == "0.1.0"
'@

Write-NewFile "tests/test_grid_and_fields.py" @'
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_toy_solution


def test_grid_shape_and_field_features() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_toy_solution(grid)

    assert grid.X.shape == (4, 5)
    assert grid.Y.shape == (4, 5)
    assert field.as_feature_matrix().shape == (20, 7)
'@

Write-NewFile "tests/test_numerics.py" @'
import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.numerics import derivative_x, derivative_y, estimate_order, l2_error, linf_error


def test_derivatives_for_linear_field() -> None:
    grid = Grid2D.from_bounds(nx=8, ny=6)
    f = 2.0 * grid.X - 3.0 * grid.Y + 1.0

    dfdx = derivative_x(f, grid.dx)
    dfdy = derivative_y(f, grid.dy)

    assert np.allclose(dfdx, 2.0)
    assert np.allclose(dfdy, -3.0)


def test_error_norms_and_order() -> None:
    exact = np.array([1.0, 2.0, 3.0])
    numerical = np.array([1.0, 2.0, 4.0])

    assert l2_error(numerical, exact) > 0.0
    assert linf_error(numerical, exact) == 1.0

    h = np.array([0.5, 0.25, 0.125])
    err = h**2
    assert np.isclose(estimate_order(h, err), 2.0)
'@

Write-NewFile "examples/01_smoke_python_package.py" @'
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_toy_solution
from mms_viscoelastic.ml import field_to_feature_table
from mms_viscoelastic.numerics import derivative_x, l2_error


def main() -> None:
    grid = Grid2D.from_bounds(nx=32, ny=32)
    field = evaluate_toy_solution(grid)
    table = field_to_feature_table(grid, field)

    dudx = derivative_x(field.u, grid.dx)
    zero_error = l2_error(field.u, field.u)

    print(table.head())
    print(f"dudx shape: {dudx.shape}")
    print(f"zero error: {zero_error:.3e}")


if __name__ == "__main__":
    main()
'@

Write-NewFile "docs/repository-map.md" @'
# Mapa inicial do repositório

Este repositório reúne os arquivos originais da tese, códigos Fortran, notebooks
Mathematica, artigo científico e a nova fundação Python.

## Regra de preservação

Os códigos Fortran e os arquivos LaTeX originais não devem ser alterados nesta
fase. Eles foram apenas reorganizados com `git mv`, preservando o histórico do Git.

## Estrutura após a reorganização

```text
src/
  fortran_original/        # código Fortran original da tese
  mathematica_original/    # notebooks Mathematica originais
  mms_viscoelastic/        # nova biblioteca Python em construção

thesis/
  original/                # tese LaTeX original

papers/
  original/                # artigo e figuras originais

docs/                      # documentação técnica do repositório
examples/                  # exemplos executáveis da biblioteca Python
tests/                     # testes automatizados da biblioteca Python
```

## Objetivo da fundação Python

A pasta `src/mms_viscoelastic` foi criada para evoluir para uma biblioteca Python.
A transcrição Fortran -> Python deve ser incremental, testável e orientada por
comparação contra a implementação original.
'@

Write-NewFile "docs/python-porting-plan.md" @'
# Plano de transcrição Fortran -> Python

## Princípio

A versão Python deve nascer como biblioteca científica, não como um conjunto solto
de scripts. O código Fortran original permanece como referência.

## Ordem sugerida de transcrição

1. Parâmetros físicos e numéricos.
2. Malhas e convenções de indexação.
3. Campos bidimensionais e saída de dados.
4. Soluções manufaturadas dos arquivos `case_01.f` até `case_14.f`.
5. Operadores diferenciais de segunda ordem.
6. Operadores compactos de alta ordem.
7. Normas de erro e estimativa de ordem de convergência.
8. Solver de Poisson.
9. Integração temporal Euler/RK.
10. Modelos UCM, Oldroyd-B, Giesekus e LPTT.
11. Comparação automatizada Fortran vs Python.
12. Geração de datasets para Machine Learning.

## Preparação para Machine Learning

A função `field_to_feature_table` já converte malhas e campos em `pandas.DataFrame`,
criando uma base natural para regressão, classificação, surrogates, PINNs/FNOs ou
outros métodos de Scientific Machine Learning.

## Cuidados críticos

- Manter a convenção NumPy: campos com shape `(ny, nx)`.
- Documentar qualquer diferença em relação aos índices Fortran.
- Testar cada módulo isoladamente.
- Não validar o Python apenas visualmente; comparar normas e ordens de convergência.
'@

Write-Info "Validando estrutura Python..."
python -m compileall src/mms_viscoelastic

Write-Info "Tentando executar testes, se pytest estiver instalado..."
try {
    python -m pytest -q
}
catch {
    Write-Warn "pytest nao executou. Instale dependencias com: python -m pip install -r requirements-dev.txt"
}

Write-Info "Status Git apos reorganizacao:"
git status --short

Write-Info "Fase 1 concluida. Revise o status e faca commit quando estiver tudo certo."
