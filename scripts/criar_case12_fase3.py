"""Fase 3: cria a primeira transcricao Python do caso MMS 12.

Execute este script na RAIZ do repositorio reorganizado.
Ele NAO altera os arquivos Fortran/LaTeX originais.  Ele cria/atualiza apenas
arquivos da biblioteca Python, testes, exemplo e documentacao tecnica.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path.cwd()


def write_text(path: str, content: str) -> None:
    target = ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content.strip() + "\n", encoding="utf-8", newline="\n")
    print(f"[OK] escrito: {path}")


def ensure_project_root() -> None:
    required = [
        ROOT / "src" / "mms_viscoelastic",
        ROOT / "src" / "fortran_original",
        ROOT / "pyproject.toml",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        raise SystemExit(
            "[ERRO] Execute este script na raiz do repositorio reorganizado. "
            f"Itens ausentes: {missing}"
        )


CASE_12_MODULE = r'''
"""Manufactured solution: Fortran ``case_12.f`` exact fields.

This module ports the exact manufactured fields from the original Fortran
subroutine ``exact_solutions_case_12(t)``.  It is intentionally focused on the
exact solution first; the very long source terms from ``source_term_case_12(t)``
should be ported in a later phase after this field-level translation is tested.

Array convention
----------------
Fortran stores fields as ``field(i, j)``.  NumPy stores the same physical field
as arrays with shape ``(ny, nx)``, where axis 0 is ``y`` and axis 1 is ``x``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from mms_viscoelastic.fields import Field2D
from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.parameters import FlowParameters

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Case12Coefficients:
    """Trigonometric stress coefficients used by MMS case 12.

    The defaults reproduce the values defined in the original ``par.nn`` file:
    all ``a`` and ``b`` coefficients are equal to one.
    """

    a1xx: float = 1.0
    a2xx: float = 1.0
    a3xx: float = 1.0
    b1xx: float = 1.0
    b2xx: float = 1.0
    b3xx: float = 1.0
    a1xy: float = 1.0
    a2xy: float = 1.0
    a3xy: float = 1.0
    b1xy: float = 1.0
    b2xy: float = 1.0
    b3xy: float = 1.0
    a1yy: float = 1.0
    a2yy: float = 1.0
    a3yy: float = 1.0
    b1yy: float = 1.0
    b2yy: float = 1.0
    b3yy: float = 1.0


@dataclass(frozen=True)
class Case12Parameters:
    """Parameters specific to the exact MMS case 12 fields."""

    beta: float = 0.5
    at: float = 0.0
    coefficients: Case12Coefficients = Case12Coefficients()

    @classmethod
    def from_flow_parameters(
        cls,
        flow: FlowParameters,
        *,
        at: float = 0.0,
        coefficients: Case12Coefficients | None = None,
    ) -> "Case12Parameters":
        """Build case parameters from the generic ``FlowParameters`` object."""
        return cls(beta=flow.beta, at=at, coefficients=coefficients or Case12Coefficients())

    def validate(self) -> None:
        """Validate the admissible range of the MMS parameters."""
        if not (0.0 < self.beta <= 1.0):
            raise ValueError("beta must satisfy 0 < beta <= 1.")


def _stress_factor_x(x: FloatArray, a1: float, a2: float, a3: float) -> FloatArray:
    return a1 + a3 * np.cos(np.pi * x) + a2 * np.sin(np.pi * x)


def _stress_factor_y(y: FloatArray, b1: float, b2: float, b3: float) -> FloatArray:
    return b1 + b3 * np.cos(np.pi * y) + b2 * np.sin(np.pi * y)


def evaluate_case_12_exact(
    grid: Grid2D,
    t: float = 0.0,
    parameters: Case12Parameters | None = None,
) -> Field2D:
    """Evaluate the exact manufactured solution of Fortran ``case_12.f``.

    Parameters
    ----------
    grid:
        Uniform Cartesian grid.  The returned arrays have shape ``(grid.ny, grid.nx)``.
    t:
        Physical time.
    parameters:
        Case-12 parameters.  If omitted, defaults mirror the thesis code for
        ``beta`` and trigonometric stress coefficients.

    Returns
    -------
    Field2D
        Exact fields ``u``, ``v``, ``wz``, ``psi``, ``txx``, ``txy`` and ``tyy``.
    """
    params = parameters or Case12Parameters()
    params.validate()

    x = grid.X
    y = grid.Y
    pi = np.pi
    decay = np.exp(params.at * t)
    c = params.coefficients

    sin_px = np.sin(pi * x)
    cos_px = np.cos(pi * x)
    sin_py = np.sin(pi * y)
    cos_py = np.cos(pi * y)

    x_poly = (-1.0 + x) * x

    u = x_poly * sin_px * (pi * y * cos_py + sin_py) / decay

    v = -(
        y
        * (pi * x_poly * cos_px + (-1.0 + 2.0 * x) * sin_px)
        * sin_py
    ) / decay

    wz = (
        2.0
        * (
            pi * x_poly * cos_py * sin_px
            + y
            * (
                pi * (-1.0 + 2.0 * x) * cos_px
                + (1.0 - pi**2 * x_poly) * sin_px
            )
            * sin_py
        )
    ) / decay

    psi = x_poly * y * sin_px * sin_py / decay

    txx = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1xx, c.a2xx, c.a3xx)
        * _stress_factor_y(y, c.b1xx, c.b2xx, c.b3xx)
        / decay
    )
    txy = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1xy, c.a2xy, c.a3xy)
        * _stress_factor_y(y, c.b1xy, c.b2xy, c.b3xy)
        / decay
    )
    tyy = (
        (1.0 - params.beta)
        * _stress_factor_x(x, c.a1yy, c.a2yy, c.a3yy)
        * _stress_factor_y(y, c.b1yy, c.b2yy, c.b3yy)
        / decay
    )

    return Field2D(u=u, v=v, psi=psi, wz=wz, txx=txx, txy=txy, tyy=tyy)


def case_12_feature_names() -> list[str]:
    """Return the canonical order of fields produced by case 12."""
    return ["u", "v", "psi", "wz", "txx", "txy", "tyy"]
'''

INIT_MODULE = r'''
"""Manufactured solutions for verification studies."""

from .case_12 import (
    Case12Coefficients,
    Case12Parameters,
    case_12_feature_names,
    evaluate_case_12_exact,
)
from .toy import evaluate_toy_solution

__all__ = [
    "Case12Coefficients",
    "Case12Parameters",
    "case_12_feature_names",
    "evaluate_case_12_exact",
    "evaluate_toy_solution",
]
'''

TEST_CASE_12 = r'''
import numpy as np

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import (
    Case12Parameters,
    evaluate_case_12_exact,
)
from mms_viscoelastic.ml import field_to_feature_table


def test_case_12_shapes_and_components() -> None:
    grid = Grid2D.from_bounds(nx=33, ny=33)
    field = evaluate_case_12_exact(grid, t=0.0)

    assert field.u.shape == (33, 33)
    assert field.v.shape == (33, 33)
    assert field.wz.shape == (33, 33)
    assert field.psi.shape == (33, 33)
    assert field.txx.shape == (33, 33)
    assert field.txy.shape == (33, 33)
    assert field.tyy.shape == (33, 33)
    assert field.available_components() == ["u", "v", "psi", "wz", "txx", "txy", "tyy"]


def test_case_12_known_boundary_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters(beta=0.5, at=0.0))

    # Point (x, y) = (0, 0), with all stress coefficients equal to one.
    assert np.isclose(field.u[0, 0], 0.0)
    assert np.isclose(field.v[0, 0], 0.0)
    assert np.isclose(field.psi[0, 0], 0.0)
    assert np.isclose(field.wz[0, 0], 0.0)
    assert np.isclose(field.txx[0, 0], 2.0)
    assert np.isclose(field.txy[0, 0], 2.0)
    assert np.isclose(field.tyy[0, 0], 2.0)


def test_case_12_known_center_values() -> None:
    grid = Grid2D.from_bounds(nx=3, ny=3)
    field = evaluate_case_12_exact(grid, t=0.0, parameters=Case12Parameters(beta=0.5, at=0.0))

    center = (1, 1)  # y-index, x-index -> (x, y) = (0.5, 0.5)
    assert np.isclose(field.u[center], -0.25)
    assert np.isclose(field.v[center], 0.0, atol=1.0e-15)
    assert np.isclose(field.psi[center], -0.125)
    assert np.isclose(field.wz[center], 1.0 + np.pi**2 / 4.0)
    assert np.isclose(field.txx[center], 2.0)
    assert np.isclose(field.txy[center], 2.0)
    assert np.isclose(field.tyy[center], 2.0)


def test_case_12_dataframe_for_ml_pipeline() -> None:
    grid = Grid2D.from_bounds(nx=5, ny=4)
    field = evaluate_case_12_exact(grid, t=0.0)
    table = field_to_feature_table(grid, field)

    assert table.shape == (20, 9)
    assert list(table.columns) == ["x", "y", "u", "v", "psi", "wz", "txx", "txy", "tyy"]
    assert np.isfinite(table.to_numpy()).all()
'''

EXAMPLE_CASE_12 = r'''
"""Example: evaluate MMS case 12 and prepare a tabular ML dataset.

Run from the repository root after installing the package in editable mode:

    python examples/02_case12_manufactured_solution.py

Optional CSV export:

    python examples/02_case12_manufactured_solution.py --export-csv outputs/case12_dataset.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

from mms_viscoelastic.grids import Grid2D
from mms_viscoelastic.manufactured_solutions import evaluate_case_12_exact
from mms_viscoelastic.ml import field_to_feature_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate MMS case 12 exact fields.")
    parser.add_argument("--nx", type=int, default=33, help="Number of grid points in x.")
    parser.add_argument("--ny", type=int, default=33, help="Number of grid points in y.")
    parser.add_argument("--time", type=float, default=0.0, help="Physical time.")
    parser.add_argument(
        "--export-csv",
        type=Path,
        default=None,
        help="Optional path to export the generated dataset as CSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = Grid2D.from_bounds(nx=args.nx, ny=args.ny)
    field = evaluate_case_12_exact(grid, t=args.time)
    table = field_to_feature_table(grid, field)

    print("MMS case 12 exact field")
    print(f"grid: nx={grid.nx}, ny={grid.ny}, dx={grid.dx:.6e}, dy={grid.dy:.6e}")
    print(table.head())
    print("\ncomponent ranges:")
    for component in field.available_components():
        values = getattr(field, component)
        print(f"  {component:>3s}: min={values.min(): .6e}, max={values.max(): .6e}")

    if args.export_csv is not None:
        args.export_csv.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.export_csv, index=False)
        print(f"\nCSV exported to: {args.export_csv}")


if __name__ == "__main__":
    main()
'''

DOC_CASE_12 = r'''
# Port Python do caso manufaturado 12

Este documento registra a primeira transcrição efetiva do código Fortran para a
biblioteca Python `mms_viscoelastic`.

## Arquivo de origem

O módulo Python foi baseado na subrotina Fortran:

```text
src/fortran_original/4Order_LPOG_Euler/case_12.f
```

A rotina transcrita nesta fase é:

```fortran
exact_solutions_case_12(t)
```

## Arquivo Python criado

```text
src/mms_viscoelastic/manufactured_solutions/case_12.py
```

## Campos implementados

A função `evaluate_case_12_exact` retorna um objeto `Field2D` com:

- `u`: velocidade na direção x;
- `v`: velocidade na direção y;
- `wz`: vorticidade;
- `psi`: função de corrente;
- `txx`: componente normal xx do tensor de tensões;
- `txy`: componente cisalhante xy do tensor de tensões;
- `tyy`: componente normal yy do tensor de tensões.

## Convenção de indexação

No Fortran, os campos são acessados como:

```fortran
campo(i,j)
```

No Python/NumPy, a convenção adotada é:

```python
campo[j, i]
```

ou seja, os arrays têm forma:

```text
(ny, nx)
```

Isso é essencial para evitar troca de eixos nos gráficos e nos datasets.

## Preparação para Machine Learning

O exemplo:

```text
examples/02_case12_manufactured_solution.py
```

permite gerar uma tabela com colunas:

```text
x, y, u, v, psi, wz, txx, txy, tyy
```

Essa tabela é a primeira base limpa para futuros estudos com Machine Learning,
como regressão de campos, operadores neurais, surrogates ou modelos híbridos.

## Limitação desta fase

A subrotina `source_term_case_12(t)` ainda não foi transcrita. Ela contém uma
expressão algébrica muito extensa e deve ser portada em fase própria, com testes
separados e comparação cuidadosa com os dados gerados pelo Fortran.

## Testes criados

```text
tests/test_case_12_manufactured.py
```

Os testes verificam:

1. formato dos campos;
2. valores conhecidos em pontos específicos;
3. conversão do campo para `pandas.DataFrame`;
4. compatibilidade inicial com pipeline de Machine Learning.
'''


def main() -> None:
    ensure_project_root()
    write_text("src/mms_viscoelastic/manufactured_solutions/case_12.py", CASE_12_MODULE)
    write_text("src/mms_viscoelastic/manufactured_solutions/__init__.py", INIT_MODULE)
    write_text("tests/test_case_12_manufactured.py", TEST_CASE_12)
    write_text("examples/02_case12_manufactured_solution.py", EXAMPLE_CASE_12)
    write_text("docs/case-12-python-port.md", DOC_CASE_12)
    print("\n[OK] Fase 3 criada. Rode agora:")
    print("  python -m pytest -q")
    print("  python examples\\02_case12_manufactured_solution.py")
    print("  python examples\\02_case12_manufactured_solution.py --export-csv outputs\\case12_dataset.csv")


if __name__ == "__main__":
    main()
