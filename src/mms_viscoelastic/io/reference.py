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
