#!/usr/bin/env Rscript
# F2 — nnSVG: Nearest Neighbor-based Spatially Variable Gene detection
# =====================================================================
# Citation: Weber et al. Nat Commun 14, 4059 (2023). doi:10.1038/s41467-023-39748-0
# Benchmark: Salas et al. Nat Methods 22, 813-823 (2025): 95.3% TPR
#
# Input: numpy files pre-exported by F2_export_for_nnsvg.py
#   counts_matrix.npy  — float32 (n_cells, n_genes)
#   gene_names.csv     — one gene per line
#   cell_names.csv     — one cell per line
#   spatial_coords.npy — float32 (n_cells, 2)

suppressPackageStartupMessages({
  library(optparse)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(jsonlite)
  library(nnSVG)
})

opt_list <- list(
  make_option("--datadir",    type="character", help="Directory with numpy export files"),
  make_option("--n_neighbors",type="integer",   default=10),
  make_option("--n_threads",  type="integer",   default=4),
  make_option("--outdir",     type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("=======================================================================\n")
cat("F2 — nnSVG Spatially Variable Gene Detection\n")
cat("=======================================================================\n")
cat(sprintf("  datadir     : %s\n", opt$datadir))
cat(sprintf("  n_neighbors : %d\n", opt$n_neighbors))
cat(sprintf("  n_threads   : %d\n", opt$n_threads))

t0 <- Sys.time()

# ─── Load pre-exported numpy files ───────────────────────────────────────────
cat("\n[1/4] Loading numpy exports...\n")

# Matrix Market file is (genes, cells) sparse
counts_sparse <- readMM(file.path(opt$datadir, "counts_matrix.mtx"))
gene_names    <- readLines(file.path(opt$datadir, "gene_names.csv"))
cell_names    <- readLines(file.path(opt$datadir, "cell_names.csv"))
spatial_df    <- read.csv(file.path(opt$datadir, "spatial_coords.csv"))
spatial_xy    <- as.matrix(spatial_df[, c("x","y")])

counts_sparse <- as(counts_sparse, "dgCMatrix")
rownames(counts_sparse) <- gene_names
colnames(counts_sparse) <- cell_names

cat(sprintf("  Counts: %d genes × %d cells\n", nrow(counts_sparse), ncol(counts_sparse)))
cat(sprintf("  Coords: %d cells × 2\n",        nrow(spatial_xy)))

# ─── Build SpatialExperiment ─────────────────────────────────────────────────
cat("\n[2/4] Building SpatialExperiment...\n")

# Log-normalize for nnSVG (requires logcounts assay)
lib_size <- Matrix::colSums(counts_sparse)
lib_size[lib_size == 0] <- 1
logcounts_mat <- log1p(t(t(counts_sparse) / lib_size) * 1e4)

spe <- SpatialExperiment(
  assays        = list(counts = counts_sparse, logcounts = logcounts_mat),
  spatialCoords = spatial_xy,
  colData       = DataFrame(cell_id = cell_names)
)
cat(sprintf("  SPE: %d genes × %d cells\n", nrow(spe), ncol(spe)))

# ─── nnSVG ───────────────────────────────────────────────────────────────────
cat(sprintf("\n[3/4] Running nnSVG (k=%d, threads=%d)...\n",
            opt$n_neighbors, opt$n_threads))

spe <- nnSVG(
  spe,
  assay_name  = "logcounts",
  n_neighbors = opt$n_neighbors,
  n_threads   = opt$n_threads,
  verbose     = FALSE
)

# nnSVG rowData output columns: sigma.sq, tau.sq, prop_sv, phi, mean, var,
#   LR_stat, rank, pval, padj
rd <- as.data.frame(rowData(spe))
rd$gene <- rownames(spe)

keep_cols <- intersect(
  c("gene", "rank", "LR_stat", "prop_sv", "pval", "padj", "mean", "var"),
  colnames(rd)
)
results_df  <- rd[, keep_cols]
results_df  <- results_df[order(results_df$rank), ]
results_df$significant <- results_df$padj < 0.05
n_sig <- sum(results_df$significant, na.rm=TRUE)

cat(sprintf("  nnSVG SVGs (padj<0.05): %d/%d genes\n", n_sig, nrow(results_df)))

# ─── Export ──────────────────────────────────────────────────────────────────
cat("\n[4/4] Saving outputs...\n")
out_csv <- file.path(opt$outdir, "nnsvg_svg_scores.csv")
write.csv(results_df, out_csv, row.names=FALSE)
cat(sprintf("  Saved: %s\n", out_csv))

cat("  Top 10 SVGs:\n")
for (i in seq_len(min(10, nrow(results_df)))) {
  r <- results_df[i,]
  cat(sprintf("    %-15s rank=%-4d LR=%.2f prop_sv=%.4f padj=%.2e\n",
              r$gene, r$rank, r$LR_stat, r$prop_sv, r$padj))
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
write_json(
  list(n_genes=nrow(results_df), n_cells=ncol(spe),
       n_neighbors=opt$n_neighbors, n_svg_significant=n_sig,
       execution_seconds=round(elapsed, 1)),
  file.path(opt$outdir, "nnsvg_summary.json"),
  pretty=TRUE, auto_unbox=TRUE
)
cat(sprintf("\n✓ nnSVG complete in %.1fs\n", elapsed))
