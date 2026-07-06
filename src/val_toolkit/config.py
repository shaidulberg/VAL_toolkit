from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def read_yaml(path: str | Path) -> dict[str, Any]:
    """Read a YAML config file and return a dictionary."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file does not exist: {path}")
    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if config is None:
        raise ValueError(f"Config file is empty: {path}")
    if not isinstance(config, dict):
        raise TypeError(f"Config file must contain a YAML mapping at top level: {path}")
    return config


def resolve_path(path: str | Path, base_dir: str | Path | None = None) -> Path:
    """Resolve a path relative to base_dir unless it is already absolute."""
    path = Path(path)
    if path.is_absolute():
        return path
    if base_dir is None:
        base_dir = Path.cwd()
    return Path(base_dir).resolve() / path
