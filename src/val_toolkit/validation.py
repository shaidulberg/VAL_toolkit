from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def require_columns(df: pd.DataFrame, columns: Iterable[str], table_name: str = "table") -> None:
    """Raise a readable error if required columns are missing."""
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns in {table_name}: {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def ensure_output_dir(path: str | Path) -> Path:
    """Create output directory if needed and return it as a Path."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def read_csv_checked(path: str | Path, required_columns: Iterable[str] | None = None) -> pd.DataFrame:
    """Read a CSV and optionally check required columns."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input CSV does not exist: {path}")
    df = pd.read_csv(path)
    if required_columns is not None:
        require_columns(df, required_columns, table_name=str(path))
    return df
