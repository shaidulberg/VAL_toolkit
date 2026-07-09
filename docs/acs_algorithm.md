# ACS algorithm specification

This document defines the ACS-ranked figure workflow implemented in this repository.

## Terms

**Annotation method**  
A cell-type annotation source, such as CellTypist, SingleR with a particular reference, SCimilarity, or a user-provided classifier.

**Raw label**  
The label emitted by an annotation method. Examples: `Memory B`, `CD4 T cell`, `Classical monocyte`, `NK cell`, `cDC`.

**Major immune lineage**  
The harmonized cell-type label used in the ACS figures:

```text
B
CD4 T
CD8 T
DC
Monocytes
NK
Other
```

**VAL**  
For a given cell, VAL is the number of annotation methods assigning the consensus major immune lineage.

**ACS ranking**  
Cells are ranked within each major lineage by a two-level confidence hierarchy:

1. VAL descending.
2. Mean confidence among methods supporting the consensus label descending.

This ranking mirrors the manuscript description of weighted VALs used to generate Aggregated Confidence Scores (ACS).

---

## Inputs

The ACS workflow requires:

1. An `.h5ad` file.
2. `adata.obs` columns for sample/patient and response metadata.
3. Annotation labels and optional confidence scores from one or more methods, either stored in `adata.obs` or supplied as an external CSV.

---

## Step-by-step algorithm

### 1. Harmonize annotation labels

For each method and each cell, the raw method label is mapped to one of:

```text
B, CD4 T, CD8 T, DC, Monocytes, NK, Other
```

The mapping uses built-in conservative synonyms and an optional user-provided `label_mapping` block in the YAML config.

### 2. Select consensus cell type

For each cell, the consensus label is selected from the harmonized method labels.

Rules:

1. Ignore `Other` when at least one non-`Other` label is available.
2. Choose the label with the most method votes.
3. If vote counts tie, choose the label with the highest mean confidence among supporting methods.
4. If all methods are `Other`, the consensus label is `Other`.

### 3. Compute VAL

For each cell:

```text
VAL = number of methods whose harmonized label equals the consensus label
```

The output also includes:

```text
VAL_fraction = VAL / number_of_annotation_methods
```

### 4. Compute mean consensus confidence

For each cell:

```text
mean_consensus_confidence = mean(confidence scores among methods supporting the consensus label)
```

If a method has no confidence column, confidence is treated as 1.0 for that method.

### 5. Rank cells by ACS within each cell type

For each major lineage separately, cells are sorted by:

```text
VAL descending
mean_consensus_confidence descending
cell_id ascending as a deterministic final tie-breaker
```

### 6. Assign non-overlapping ACS bins

Within each major lineage, cells are divided into non-overlapping 10% bins by ACS rank:

```text
0-10, 10-20, 20-30, ..., 90-100
```

### 7. Compute bin abundance per sample

For each sample, cell type, and ACS bin, the workflow counts cells and calculates a proportion.

The denominator is controlled by config:

```yaml
denominator: all_cells
```

or:

```yaml
denominator: cell_type_cells
```

The default is `all_cells`, matching an abundance-style interpretation.

### 8. Test response association

For each cell type and ACS bin, the workflow compares bin abundance between response groups using two-sided Mann-Whitney U tests.

FDR correction can be global across all cell-type-by-bin tests or within each cell type:

```yaml
statistics:
  fdr_scope: global
```

or:

```yaml
statistics:
  fdr_scope: within_cell_type
```

### 9. Plot ACS figures A-F

The plotted y-axis is:

```text
-log10(FDR q-value)
```

The x-axis is the ACS-ranked bin endpoint:

```text
10, 20, 30, ..., 100
```

The vertical dashed line marks the bin with the smallest q-value for that cell type.

---

## Core output files

```text
acs_cell_table.csv
acs_ranked_cell_bins.csv
acs_ranked_bin_abundance_by_sample.csv
acs_ranked_bin_association_results.csv
acs_figures_A_to_F.png/pdf/svg
```

`acs_cell_table.csv` is the central table for inspecting method votes, consensus labels, VAL, and ACS ranking.
