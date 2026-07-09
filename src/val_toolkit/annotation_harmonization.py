from __future__ import annotations

from typing import Any, Mapping, Sequence

import pandas as pd

CANONICAL_ACS_CELL_TYPES: tuple[str, ...] = ("B", "CD4 T", "CD8 T", "DC", "Monocytes", "NK")

DEFAULT_LABEL_SYNONYMS: dict[str, str] = {
    # B lineage
    "b": "B",
    "b cell": "B",
    "b cells": "B",
    "b_cell": "B",
    "b_cells": "B",
    "b-cell": "B",
    "memory b": "B",
    "memory b cell": "B",
    "naive b": "B",
    "naive b cell": "B",
    "plasma": "B",
    "plasma cell": "B",
    "plasmablast": "B",
    # CD4 T lineage
    "cd4": "CD4 T",
    "cd4 t": "CD4 T",
    "cd4 t cell": "CD4 T",
    "cd4 t cells": "CD4 T",
    "cd4-positive t cell": "CD4 T",
    "helper t": "CD4 T",
    "t helper": "CD4 T",
    "treg": "CD4 T",
    "regulatory t": "CD4 T",
    # CD8 T lineage
    "cd8": "CD8 T",
    "cd8 t": "CD8 T",
    "cd8 t cell": "CD8 T",
    "cd8 t cells": "CD8 T",
    "cd8-positive t cell": "CD8 T",
    "cytotoxic t": "CD8 T",
    "cytotoxic t cell": "CD8 T",
    # Dendritic cell lineage
    "dc": "DC",
    "dendritic": "DC",
    "dendritic cell": "DC",
    "dendritic cells": "DC",
    "cdc": "DC",
    "cdc1": "DC",
    "cdc2": "DC",
    "conventional dendritic cell": "DC",
    "plasmacytoid dendritic cell": "DC",
    "pdc": "DC",
    "lamp3 dc": "DC",
    "mregdc": "DC",
    # Monocyte lineage
    "mono": "Monocytes",
    "monocyte": "Monocytes",
    "monocytes": "Monocytes",
    "classical monocyte": "Monocytes",
    "non-classical monocyte": "Monocytes",
    "intermediate monocyte": "Monocytes",
    # NK lineage
    "nk": "NK",
    "nk cell": "NK",
    "nk cells": "NK",
    "natural killer": "NK",
    "natural killer cell": "NK",
    "natural killer cells": "NK",
}


def normalize_label_text(label: Any) -> str:
    """Normalize label text for conservative matching."""
    if pd.isna(label):
        return ""
    return (
        str(label)
        .strip()
        .lower()
        .replace("_", " ")
        .replace("-", " ")
        .replace("/", " ")
    )


def flatten_label_mapping(label_mapping: Mapping[str, Any] | None) -> dict[str, str]:
    """
    Convert config label mapping into raw-label -> canonical-label mapping.

    Supports either of these YAML styles:
      label_mapping:
        B cells: B
        CD4 T cells: CD4 T

    or:
      label_mapping:
        B: [B cell, B cells, Naive B]
        CD4 T: [CD4 T cell, Treg]
    """
    flat: dict[str, str] = {}
    if not label_mapping:
        return flat

    canonical = set(CANONICAL_ACS_CELL_TYPES) | {"Other"}
    for key, value in label_mapping.items():
        key_str = str(key)
        if isinstance(value, (list, tuple, set)):
            if key_str not in canonical:
                raise ValueError(
                    "List-style label_mapping must use canonical cell types as keys. "
                    f"Got key={key_str!r}."
                )
            for raw in value:
                flat[str(raw)] = key_str
        else:
            flat[key_str] = str(value)
    return flat


def harmonize_label(
    label: Any,
    label_mapping: Mapping[str, Any] | None = None,
    allow_substring_match: bool = True,
) -> str:
    """Map a raw annotation label to B/CD4 T/CD8 T/DC/Monocytes/NK/Other."""
    if pd.isna(label):
        return "Other"

    text = str(label).strip()
    if text == "":
        return "Other"

    flat = flatten_label_mapping(label_mapping)
    if text in flat:
        return flat[text]

    normalized_flat = {normalize_label_text(k): v for k, v in flat.items()}
    normalized = normalize_label_text(text)
    if normalized in normalized_flat:
        return normalized_flat[normalized]

    if normalized in DEFAULT_LABEL_SYNONYMS:
        return DEFAULT_LABEL_SYNONYMS[normalized]

    if allow_substring_match:
        # Use longer keys first to avoid short tokens like "b" matching too broadly.
        for key in sorted(DEFAULT_LABEL_SYNONYMS, key=len, reverse=True):
            if len(key) >= 3 and key in normalized:
                return DEFAULT_LABEL_SYNONYMS[key]

    return "Other"


def validate_major_cell_types(cell_types: Sequence[str]) -> list[str]:
    invalid = [ct for ct in cell_types if ct not in CANONICAL_ACS_CELL_TYPES]
    if invalid:
        raise ValueError(
            "Unsupported major cell type(s): "
            + ", ".join(invalid)
            + f". Supported values are: {', '.join(CANONICAL_ACS_CELL_TYPES)}"
        )
    return list(cell_types)
