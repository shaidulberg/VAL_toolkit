#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import anndata as ad
except ImportError as exc:
    raise ImportError(
        "This script requires anndata. Install with `pip install anndata` or recreate the conda env from environment.yml."
    ) from exc

from val_toolkit.val_panels import DEFAULT_MARKER_SETS

rng = np.random.default_rng(7)

cell_types = ["B", "CD4 T", "CD8 T", "DC", "Monocytes", "NK"]
patients = [f"P{i:02d}" for i in range(1, 17)]
responses = {patient: ("R" if i < 8 else "NR") for i, patient in enumerate(patients)}

# Raw labels intentionally mimic names that can appear across different annotation tools.
raw_label_options = {
    "B": ["B cell", "B cells", "Memory B", "Naive B"],
    "CD4 T": ["CD4 T cell", "CD4 T cells", "Helper T", "Treg"],
    "CD8 T": ["CD8 T cell", "CD8 T cells", "Cytotoxic T cell"],
    "DC": ["Dendritic cell", "cDC", "pDC", "LAMP3 DC"],
    "Monocytes": ["Monocyte", "Classical monocyte", "Non-classical monocyte"],
    "NK": ["NK cell", "Natural killer cell", "NK cells"],
}
method_names = ["celltypist", "singler_hpca", "singler_monaco", "scimilarity"]
method_accuracy = {
    "celltypist": 0.90,
    "singler_hpca": 0.84,
    "singler_monaco": 0.82,
    "scimilarity": 0.86,
}

all_markers = []
for genes in DEFAULT_MARKER_SETS.values():
    all_markers.extend(genes)
all_markers = list(dict.fromkeys(all_markers))
background_genes = [f"GENE{i:03d}" for i in range(1, 121)]
genes = all_markers + background_genes

gene_index = {gene: i for i, gene in enumerate(genes)}
rows = []
obs_records = []

cell_id_counter = 0
for patient in patients:
    response = responses[patient]
    for cell_type in cell_types:
        base_n = 36
        # Create a small response-associated enrichment for the example dataset.
        if response == "R" and cell_type == "B":
            n_cells = base_n + 16
        elif response == "R" and cell_type == "CD8 T":
            n_cells = base_n + 8
        else:
            n_cells = base_n
        for _ in range(n_cells):
            expr = rng.poisson(0.05, size=len(genes)).astype(float)
            for marker in DEFAULT_MARKER_SETS[cell_type]:
                if marker in gene_index:
                    expr[gene_index[marker]] += rng.poisson(3.5)
            # Add a little cross-lineage noise.
            noise_type = rng.choice(cell_types)
            for marker in rng.choice(
                DEFAULT_MARKER_SETS[noise_type], size=min(2, len(DEFAULT_MARKER_SETS[noise_type])), replace=False
            ):
                if marker in gene_index:
                    expr[gene_index[marker]] += rng.poisson(0.6)

            obs_record = {
                "cell_id": f"cell_{cell_id_counter:05d}",
                "patient_id": patient,
                "timepoint": "pre",
                "response": response,
                "truth_cell_type_for_example_only": cell_type,
            }

            # Add fake multi-method annotations so the ACS workflow can be tested without
            # installing CellTypist, SingleR, or SCimilarity. These columns are examples
            # of user-supplied annotation-method outputs, not real package outputs.
            for method in method_names:
                is_correct = rng.random() < method_accuracy[method]
                if is_correct:
                    called_type = cell_type
                    confidence = rng.uniform(0.72, 0.98)
                else:
                    called_type = rng.choice([ct for ct in cell_types if ct != cell_type])
                    confidence = rng.uniform(0.35, 0.75)
                obs_record[f"{method}_label"] = rng.choice(raw_label_options[called_type])
                obs_record[f"{method}_confidence"] = round(float(confidence), 4)

            rows.append(expr)
            obs_records.append(obs_record)
            cell_id_counter += 1

X = np.vstack(rows)
obs = pd.DataFrame(obs_records).set_index("cell_id")
var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
adata = ad.AnnData(X=sparse.csr_matrix(X), obs=obs, var=var)

out = Path("example_data/val_ranked_bins_example.h5ad")
out.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out)
print(f"Saved example h5ad: {out.resolve()}")
print(f"Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print("Example obs annotation columns added:")
for method in method_names:
    print(f"  - {method}_label / {method}_confidence")
