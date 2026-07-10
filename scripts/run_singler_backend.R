#!/usr/bin/env Rscript

# SingleR backend bridge for VAL_toolkit.
# Input is a MatrixMarket cell-by-gene matrix plus cell and gene ID files.
# Output is a compact CSV: cell_id, label, confidence, delta_next.

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Expected argument beginning with --, got: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg, call. = FALSE)
  }
}

normalize_ref_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

load_reference <- function(reference_name) {
  key <- normalize_ref_name(reference_name)
  if (key %in% c("hpca", "humanprimarycellatlas")) {
    return(celldex::HumanPrimaryCellAtlasData())
  }
  if (key %in% c("monaco", "monacoimmune")) {
    return(celldex::MonacoImmuneData())
  }
  if (key %in% c("dice", "databaseimmunecellexpression", "databaseimmunecellexpressiondata")) {
    return(celldex::DatabaseImmuneCellExpressionData())
  }
  if (key %in% c("blueprintencode", "encode", "blueprint")) {
    return(celldex::BlueprintEncodeData())
  }
  if (key %in% c("hema", "novershtern", "novershternhematopoietic")) {
    return(celldex::NovershternHematopoieticData())
  }
  stop(
    "Unknown SingleR reference '", reference_name, "'. Supported aliases: ",
    "HPCA, Monaco, DICE, BlueprintEncode/ENCODE, Novershtern/HEMA.",
    call. = FALSE
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

required <- c("matrix", "cell_ids", "gene_ids", "output_csv", "reference", "label_field", "prune_labels", "score_column")
missing <- required[!required %in% names(args)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "), call. = FALSE)
}

require_pkg("SingleR")
require_pkg("celldex")
require_pkg("Matrix")

matrix_path <- args[["matrix"]]
cell_ids_path <- args[["cell_ids"]]
gene_ids_path <- args[["gene_ids"]]
output_csv <- args[["output_csv"]]
reference_name <- args[["reference"]]
label_field <- args[["label_field"]]
prune_labels <- tolower(args[["prune_labels"]]) %in% c("true", "1", "yes")
score_column <- args[["score_column"]]
assay_type_ref <- if ("assay_type_ref" %in% names(args)) args[["assay_type_ref"]] else "logcounts"

message("[SingleR:R] Reading query matrix")
query_cells_by_genes <- Matrix::readMM(matrix_path)
cell_ids <- readLines(cell_ids_path, warn = FALSE)
gene_ids <- readLines(gene_ids_path, warn = FALSE)

if (nrow(query_cells_by_genes) != length(cell_ids)) {
  stop("Matrix row count does not match number of cell IDs.", call. = FALSE)
}
if (ncol(query_cells_by_genes) != length(gene_ids)) {
  stop("Matrix column count does not match number of gene IDs.", call. = FALSE)
}

rownames(query_cells_by_genes) <- make.unique(cell_ids)
colnames(query_cells_by_genes) <- make.unique(gene_ids)
query_genes_by_cells <- Matrix::t(query_cells_by_genes)

message("[SingleR:R] Loading reference: ", reference_name)
ref <- load_reference(reference_name)
if (!label_field %in% colnames(SummarizedExperiment::colData(ref))) {
  stop(
    "Requested label_field '", label_field, "' was not found in reference colData. Available fields: ",
    paste(colnames(SummarizedExperiment::colData(ref)), collapse = ", "),
    call. = FALSE
  )
}
labels <- SummarizedExperiment::colData(ref)[[label_field]]

message("[SingleR:R] Running SingleR")
pred <- SingleR::SingleR(
  test = query_genes_by_cells,
  ref = ref,
  labels = labels,
  assay.type.ref = assay_type_ref
)

assigned_labels <- as.character(pred$labels)
if (prune_labels && "pruned.labels" %in% colnames(pred)) {
  pruned <- as.character(pred$pruned.labels)
  assigned_labels[!is.na(pruned) & pruned != ""] <- pruned[!is.na(pruned) & pruned != ""]
}

scores <- as.data.frame(pred$scores)
assigned_score <- rep(NA_real_, length(assigned_labels))
for (i in seq_along(assigned_labels)) {
  label <- assigned_labels[[i]]
  if (!is.na(label) && label %in% colnames(scores)) {
    assigned_score[[i]] <- as.numeric(scores[i, label])
  } else {
    assigned_score[[i]] <- suppressWarnings(max(as.numeric(scores[i, ]), na.rm = TRUE))
  }
}

delta_next <- rep(NA_real_, length(assigned_labels))
if ("delta.next" %in% colnames(pred)) {
  delta_next <- as.numeric(pred$delta.next)
}

confidence <- assigned_score
if (score_column == "delta.next") {
  confidence <- delta_next
} else if (score_column != "assigned_label_score") {
  stop("Unsupported score_column: ", score_column, call. = FALSE)
}

out <- data.frame(
  cell_id = cell_ids,
  label = assigned_labels,
  confidence = confidence,
  delta_next = delta_next,
  stringsAsFactors = FALSE
)

utils::write.csv(out, output_csv, row.names = FALSE, quote = TRUE)
message("[SingleR:R] Saved: ", output_csv)
