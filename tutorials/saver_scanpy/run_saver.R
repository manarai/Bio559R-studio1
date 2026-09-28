#!/usr/bin/env Rscript
# Run optional SAVER expression recovery on post-QC UMI RNA counts.
# Input: Matrix Market count matrix (genes x cells) plus genes/barcodes text files.
# Output: an RDS saver object with estimate, se, and run metadata.
#
# Usage:
#   Rscript run_saver.R <counts.mtx> <genes.tsv> <barcodes.tsv> <output.rds> [ncores]
#
# Do not pass log-normalized data or ATAC peaks to SAVER.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop(
    "Usage: Rscript run_saver.R <counts.mtx> <genes.tsv> <barcodes.tsv> <output.rds> [ncores]",
    call. = FALSE
  )
}

counts_path <- args[[1]]
genes_path <- args[[2]]
barcodes_path <- args[[3]]
output_path <- args[[4]]
ncores <- if (length(args) == 5) as.integer(args[[5]]) else 1L
if (is.na(ncores) || ncores < 1L) stop("ncores must be a positive integer.", call. = FALSE)

for (path in c(counts_path, genes_path, barcodes_path)) {
  if (!file.exists(path)) stop(sprintf("Input not found: %s", path), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Matrix)
  library(SAVER)
})

counts <- readMM(counts_path) # required orientation: genes x cells
genes <- read.delim(genes_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)[[1]]

if (nrow(counts) != length(genes) || ncol(counts) != length(barcodes)) {
  stop("Matrix dimensions do not match the supplied genes/barcodes files.", call. = FALSE)
}
if (anyDuplicated(genes)) {
  stop("Gene names must be unique before running SAVER.", call. = FALSE)
}

rownames(counts) <- genes
colnames(counts) <- barcodes
message(sprintf("Running SAVER on %d genes x %d cells using %d core(s).", nrow(counts), ncol(counts), ncores))
message("Input must be post-QC, unnormalized UMI counts. This may require substantial memory and time.")

result <- SAVER::saver(counts, ncores = ncores)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    saver = result,
    input = list(
      counts_path = normalizePath(counts_path),
      genes_path = normalizePath(genes_path),
      barcodes_path = normalizePath(barcodes_path),
      ncores = ncores,
      completed_at = as.character(Sys.time()),
      package_version = as.character(utils::packageVersion("SAVER"))
    )
  ),
  output_path
)
message(sprintf("Wrote SAVER result to %s", output_path))
