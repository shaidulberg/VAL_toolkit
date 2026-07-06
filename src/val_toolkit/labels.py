from __future__ import annotations

CELL_TYPE_DISPLAY = {
    "B": "B",
    "CD4": "CD4 T",
    "CD4 T": "CD4 T",
    "CD8": "CD8 T",
    "CD8 T": "CD8 T",
    "CD4/8": "CD4/8 T",
    "CD4/8 T": "CD4/8 T",
    "Mono": "Monocytes",
    "Monocytes": "Monocytes",
    "NK": "NK",
    "DC": "DC",
    "Other": "Other",
}


def display_cell_type(label: str) -> str:
    """Convert internal cell-type labels to manuscript-friendly display labels."""
    return CELL_TYPE_DISPLAY.get(str(label), str(label))


def display_method(label: str) -> str:
    """Convert method labels to manuscript-friendly display labels."""
    label = str(label)
    lower = label.lower()
    if lower.startswith("singler."):
        return "SingleR." + label.split(".", 1)[1]
    if lower == "singler":
        return "SingleR"
    if lower == "celltypist":
        return "CellTypist"
    if lower == "scimilarity":
        return "SCimilarity"
    return label
