#!/usr/bin/env Rscript
# F1 — BANKSY Spatial Domain Identification
# ==========================================
# Cross-validation tool for Novae domains (best in Salas 2025 Xenium benchmark).
#
# Loads counts + coords from H5 files (pre-exported by F1_export_for_banksy.py)
# and runs the BANKSY pipeline:
#   1. computeBanksy with k_geom=6 (recommended for Xenium)
#   2. runBanksyPCA at lambda=0.2 (compromise of expression vs spatial)
#   3. clusterBanksy with leiden algo (resolution-tuned to target n_clusters)
#
# Citation: Singhal V., Chou N., Lee J. et al. BANKSY unifies cell typing and
# tissue domain segmentation for scalable spatial omics data analysis.
# Nat Genet 56, 431-441 (2024). doi:10.1038/s41588-024-01664-3

suppressPackageStartupMessages({
  library(optparse)
  library(Banksy)
  library(SpatialExperiment)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(jsonlite)
  library(rhdf5)
})

# ─── CLI ────────────────────────────────────────────────────────────────────
opt_list <- list(
  make_option("--counts_h5", type="character",
              help="HDF5 file with counts matrix (cells x genes)"),
  make_option("--coords_csv", type="character",
              help="CSV with cell_id, x, y"),
  make_option("--lambda", type="numeric", default=0.2,
              help="BANKSY mixing param [%default]"),
  make_option("--k_geom", type="integer", default=6,
              help="k-nearest geometric neighbors [%default]"),
  make_option("--n_clusters", type="integer", default=10),
  make_option("--algo", type="character", default="leiden"),
  make_option("--n_pcs", type="integer", default=20),
  make_option("--seed", type="integer", default=42),
  make_option("--outdir", type="character", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("=", strrep("=", 68), "\n", sep="")
cat("F1 — BANKSY Spatial Domain Identification\n")
cat("=", strrep("=", 68), "\n", sep="")
cat(sprintf("  counts_h5 : %s\n", opt$counts_h5))
cat(sprintf("  coords_csv: %s\n", opt$coords_csv))
stopifnot(file.exists(opt$counts_h5), file.exists(opt$coords_csv))
cat(sprintf("  lambda    : %s\n", opt$lambda))
cat(sprintf("  k_geom    : %s\n", opt$k_geom))
cat(sprintf("  n_clusters: %s\n", opt$n_clusters))

t0 <- Sys.time()

# ─── Load counts (cells x genes -> transpose to genes x cells) ─────────────
cat("\n[1/5] Loading counts from HDF5 ...\n")
# rhdf5::h5read in older versions requires file path explicitly
counts_h5_path <- normalizePath(opt$counts_h5)
cat(sprintf("  Reading: %s\n", counts_h5_path))
counts_data <- h5read(counts_h5_path, "X")
gene_names <- h5read(counts_h5_path, "var_names")
cell_names <- h5read(counts_h5_path, "obs_names")
gene_names <- as.character(gene_names)
cell_names <- as.character(cell_names)

# counts_data is (n_genes, n_cells) when stored row-major from python
if (is.null(dim(counts_data))) {
  stop("counts_data must be a 2D matrix")
}
counts <- as(counts_data, "dgCMatrix")
# Verify orientation: rows=genes, cols=cells
if (nrow(counts) == length(cell_names) && ncol(counts) == length(gene_names)) {
  counts <- t(counts)
}
rownames(counts) <- gene_names
colnames(counts) <- cell_names
cat(sprintf("  counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# ─── Coordinates ───────────────────────────────────────────────────────────
coord_df <- read.csv(normalizePath(opt$coords_csv))
stopifnot(all(c("cell_id", "x", "y") %in% colnames(coord_df)))
rownames(coord_df) <- coord_df$cell_id
coords <- as.matrix(coord_df[cell_names, c("x", "y")])
cat(sprintf("  coords: %s\n", paste(dim(coords), collapse=" x ")))

# ─── Build SpatialExperiment ───────────────────────────────────────────────
cat("\n[2/5] Building SpatialExperiment...\n")
spe <- SpatialExperiment(
  assays = list(counts = counts),
  spatialCoords = coords,
  colData = DataFrame(cell_id = cell_names)
)
# Banksy expects a normalized assay
libsize <- colSums(counts)
medsize <- median(libsize)
logcounts_mat <- log1p(t(t(counts) / libsize) * medsize)
assay(spe, "logcounts") <- logcounts_mat
cat("  spe ready.\n")

# ─── computeBanksy ─────────────────────────────────────────────────────────
cat(sprintf("\n[3/5] computeBanksy(k_geom=%d, lambda=%.2f)...\n",
            opt$k_geom, opt$lambda))
spe <- computeBanksy(
  spe,
  assay_name = "logcounts",
  k_geom = opt$k_geom,
  compute_agf = FALSE,
  verbose = TRUE,
  seed = opt$seed
)
cat("  done.\n")

# ─── PCA ───────────────────────────────────────────────────────────────────
cat(sprintf("\n[4/5] runBanksyPCA(npcs=%d, lambda=%.2f)...\n",
            opt$n_pcs, opt$lambda))
spe <- runBanksyPCA(
  spe,
  use_agf = FALSE,
  lambda = opt$lambda,
  npcs = opt$n_pcs,
  seed = opt$seed
)
cat("  done.\n")

# ─── Cluster ───────────────────────────────────────────────────────────────
cat(sprintf("\n[5/5] clusterBanksy(algo=%s, target n_clusters≈%d)...\n",
            opt$algo, opt$n_clusters))
res_attempts <- c(0.3, 0.5, 0.8, 1.0, 1.5)
best_spe <- NULL
best_diff <- Inf
for (res in res_attempts) {
  spe_try <- clusterBanksy(
    spe, use_agf = FALSE, lambda = opt$lambda,
    algo = opt$algo,
    resolution = res,
    k_neighbors = 50, npcs = opt$n_pcs, seed = opt$seed
  )
  cluster_col <- tail(clusterNames(spe_try), 1)
  k <- length(unique(colData(spe_try)[[cluster_col]]))
  cat(sprintf("    resolution=%.2f -> %d clusters\n", res, k))
  d <- abs(k - opt$n_clusters)
  if (d < best_diff) { best_diff <- d; best_spe <- spe_try }
  if (k >= opt$n_clusters) break
}
spe <- best_spe
cluster_col <- tail(clusterNames(spe), 1)
cat(sprintf("  Selected: %s\n", cluster_col))

# ─── Export ────────────────────────────────────────────────────────────────
out_csv <- file.path(opt$outdir, "banksy_domains.csv")
domains <- colData(spe)[[cluster_col]]
out_df <- data.frame(
  cell_id = colnames(spe),
  banksy_domain = paste0("B", domains)
)
write.csv(out_df, out_csv, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", out_csv))
cat(sprintf("  Banksy domains found: %d\n", length(unique(domains))))

# Summary JSON
summary_json <- file.path(opt$outdir, "F1_banksy_summary.json")
elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
write_json(list(
  algo = opt$algo,
  lambda = opt$lambda,
  k_geom = opt$k_geom,
  n_clusters_target = opt$n_clusters,
  n_clusters_found = length(unique(domains)),
  n_cells = ncol(spe),
  execution_seconds = round(elapsed, 1)
), summary_json, pretty=TRUE, auto_unbox=TRUE)

cat(sprintf("\n✓ F1 (Banksy) complete in %.1fs\n", elapsed))
