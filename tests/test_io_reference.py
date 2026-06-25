from pathlib import Path

import pytest

from mms_viscoelastic.io.reference import read_reference_csv

REPO_ROOT = Path(__file__).resolve().parents[1]
SAMPLE_CSV = REPO_ROOT / "validation" / "reference_data" / "case12_fields_3x3_default.csv"


def test_read_reference_csv_returns_expected_columns() -> None:
    df = read_reference_csv(SAMPLE_CSV)

    assert list(df.columns) == ["x", "y", "u", "v", "wz", "psi", "txx", "txy", "tyy"]
    assert len(df) == 9


def test_read_reference_csv_accepts_str_path() -> None:
    df = read_reference_csv(str(SAMPLE_CSV))

    assert len(df) == 9


def test_read_reference_csv_missing_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "does_not_exist.csv"

    with pytest.raises(FileNotFoundError):
        read_reference_csv(missing)
