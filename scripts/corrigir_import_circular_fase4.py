"""Corrige o import circular criado na Fase 4.

Execute na raiz do repositorio, depois de rodar criar_visualizacao_dataset_fase4.py.

Problema corrigido:
- src/mms_viscoelastic/ml/datasets.py importava field_to_feature_table via
  `from mms_viscoelastic.ml import ...`, enquanto `ml/__init__.py` importava
  datasets. Isso gerava inicializacao circular.

Solucao:
- datasets.py passa a importar diretamente de mms_viscoelastic.ml.features;
- ml/__init__.py exporta primeiro field_to_feature_table e depois os builders.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path.cwd()
PKG = ROOT / "src" / "mms_viscoelastic"
ML_DIR = PKG / "ml"


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content.strip() + "\n", encoding="utf-8", newline="\n")
    print(f"[OK] escrito: {path.relative_to(ROOT)}")


def ensure_root() -> None:
    required = [
        ROOT / "pyproject.toml",
        PKG / "__init__.py",
        ML_DIR / "datasets.py",
        ML_DIR / "__init__.py",
    ]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        raise SystemExit(
            "[ERRO] Execute este script na raiz do repositorio apos a Fase 4. "
            f"Arquivos ausentes: {missing}"
        )


def ensure_features_module() -> None:
    features_path = ML_DIR / "features.py"
    if features_path.exists():
        print(f"[OK] mantido: {features_path.relative_to(ROOT)}")
        return

    write_text(
        features_path,
        r'''
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
''',
    )


def patch_datasets_import() -> None:
    datasets_path = ML_DIR / "datasets.py"
    text = datasets_path.read_text(encoding="utf-8")
    old = "from mms_viscoelastic.ml import field_to_feature_table"
    new = "from mms_viscoelastic.ml.features import field_to_feature_table"
    if old in text:
        text = text.replace(old, new)
        datasets_path.write_text(text, encoding="utf-8", newline="\n")
        print(f"[OK] corrigido import em: {datasets_path.relative_to(ROOT)}")
    elif new in text:
        print(f"[SKIP] import ja estava correto em: {datasets_path.relative_to(ROOT)}")
    else:
        raise SystemExit(
            "[ERRO] Nao encontrei o import esperado em datasets.py. "
            "Abra o arquivo e verifique manualmente."
        )


def rewrite_ml_init() -> None:
    init_path = ML_DIR / "__init__.py"
    write_text(
        init_path,
        r'''
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
''',
    )


def main() -> None:
    ensure_root()
    ensure_features_module()
    patch_datasets_import()
    rewrite_ml_init()
    print("\n[OK] Import circular corrigido.")
    print("Agora rode:")
    print("  python -m pip install -e .")
    print("  python -m pytest -q")
    print("  python examples\\03_case12_visualization_and_ml_dataset.py")


if __name__ == "__main__":
    main()
