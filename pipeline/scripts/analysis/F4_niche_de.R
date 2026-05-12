#!/usr/bin/env Rscript
# F4 — NicheDE: Niche-specific Differential Expression
# =======================================================
#
# Identifies genes that are differentially expressed in a cell type
# ONLY when that cell type is located within a specific spatial domain/niche.
#
# Method: NicheDE (Mason JC et al. Genome Biol 25, 107 (2024))
# Inputs:
#   - counts_matrix.mtx   : genes × cells raw counts
#   - spatial_coords.csv  : x, y per cell
#   - cell_types.csv      : cell_id, cell_type (L1), cell_type_fine (L2)
#   - cell_domains.csv    : cell_id, banksy_domain (14 domains from Phase B)
#
# Design: use L1 cell types (7 types) × 14 Banksy domains
#         Focus: which genes in e.g. Macrophages change based on domain context?

suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(nicheDE)
  library(jsonlite)
  library(dplyr)
})

opt_list <- list(
  make_option("--datadir",   type="character", help="Directory with exported files"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--cell_type_col", type="character", default="cell_type",
              help="Column to use: cell_type (L1) or cell_type_fine (L2) [%default]"),
  make_option("--sigma",     type="numeric",   default=100,
              help="Bandwidth for effective niche (µm) [%default]"),
  make_option("--num_cores", type="integer",   default=4,
              help="Parallel cores [%default]"),
  make_option("--C",         type="integer",   default=150,
              help="Min cells per niche [%default]"),
  make_option("--M",         type="integer",   default=10,
              help="Min mean expression [%default]"),
  make_option("--max_cells", type="integer",   default=50000,
              help="Max cells to use (stratified subsample; 0=all) [%default]")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

cat("=======================================================================\n")
cat("F4 — NicheDE: Niche-specific Differential Expression\n")
cat("=======================================================================\n")
cat(sprintf("  datadir       : %s\n", opt$datadir))
cat(sprintf("  cell_type_col : %s\n", opt$cell_type_col))
cat(sprintf("  sigma         : %g µm\n", opt$sigma))
cat(sprintf("  num_cores     : %d\n", opt$num_cores))

t0 <- Sys.time()

# ─── Load data ────────────────────────────────────────────────────────────────
cat("\n[1/5] Loading pre-exported data...\n")

counts_sparse <- readMM(file.path(opt$datadir, "counts_matrix.mtx"))
gene_names    <- readLines(file.path(opt$datadir, "gene_names.csv"))
cell_names    <- readLines(file.path(opt$datadir, "cell_names.csv"))
coords_df     <- read.csv(file.path(opt$datadir, "spatial_coords.csv"))
ct_df         <- read.csv(file.path(opt$datadir, "cell_types.csv"))
domain_df     <- read.csv(file.path(opt$datadir, "cell_domains.csv"))

# Matrix Market is genes × cells; NicheDE needs cells × genes → transpose
counts_sparse <- as(t(counts_sparse), "dgCMatrix")  # cells × genes
rownames(counts_sparse) <- cell_names
colnames(counts_sparse) <- gene_names

cat(sprintf("  Counts: %d cells × %d genes\n", nrow(counts_sparse), ncol(counts_sparse)))
cat(sprintf("  Coords: %d cells\n", nrow(coords_df)))

# ─── Subsample (SOCK parallel requires data small enough to serialize) ────────
if (opt$max_cells > 0 && nrow(counts_sparse) > opt$max_cells) {
  cat(sprintf("\n  Subsampling %d → %d cells (random, seed=42)...\n",
              nrow(counts_sparse), opt$max_cells))
  set.seed(42)
  keep          <- sort(sample(seq_len(nrow(counts_sparse)), opt$max_cells))
  counts_sparse <- counts_sparse[keep, ]
  cell_names    <- cell_names[keep]
  rownames(counts_sparse) <- cell_names
  cat(sprintf("  Subsampled: %d cells × %d genes\n",
              nrow(counts_sparse), ncol(counts_sparse)))
}

# ─── Coordinate matrix ────────────────────────────────────────────────────────
cat("\n[2/5] Preparing inputs...\n")

# Align everything to counts ROW order (cells × genes, so rows are cells)
# Force cell_id to character in all data frames to match rownames (strings)
cell_order         <- rownames(counts_sparse)
coords_df$cell_id  <- as.character(coords_df$cell_id)
ct_df$cell_id      <- as.character(ct_df$cell_id)
domain_df$cell_id  <- as.character(domain_df$cell_id)

coords_aligned <- coords_df[match(cell_order, coords_df$cell_id), c("x","y")]
ct_aligned     <- ct_df[match(cell_order, ct_df$cell_id), ]
domain_aligned <- domain_df[match(cell_order, domain_df$cell_id), ]

# Verify no NAs in coordinates
if (any(is.na(coords_aligned))) {
  n_na <- sum(is.na(coords_aligned$x))
  stop(sprintf("Coordinate alignment failed: %d cells have NA coords", n_na))
}

coord_mat <- as.matrix(coords_aligned)
rownames(coord_mat) <- cell_order
# Force numeric and verify
coord_mat[,1] <- as.numeric(coord_mat[,1])
coord_mat[,2] <- as.numeric(coord_mat[,2])
na_count <- sum(is.na(coord_mat))
if (na_count > 0) stop(sprintf("%d NA values in coord_mat", na_count))
cat(sprintf("  coord_mat: %d × %d\n", nrow(coord_mat), ncol(coord_mat)))

# ─── Deconvolution matrix (binary cell type indicator) ───────────────────────
# NicheDE expects deconv_mat: cells × cell_types (proportion or binary)
ct_labels  <- as.character(ct_aligned[[opt$cell_type_col]])
ct_labels[is.na(ct_labels)] <- "Unknown"

# Exclude "Unknown" and "Unassigned" for clean analysis
valid_mask  <- !ct_labels %in% c("Unknown", "Unassigned")
cat(sprintf("  Cells with valid cell type: %d/%d\n",
            sum(valid_mask), length(valid_mask)))

# Subset to valid cells (rows in cells × genes matrix)
counts_valid <- counts_sparse[valid_mask, ]
coord_valid  <- coord_mat[valid_mask, ]
ct_valid     <- ct_labels[valid_mask]
domain_valid <- as.character(domain_aligned$banksy_domain[valid_mask])
domain_valid[is.na(domain_valid)] <- "Unknown"

# Pre-scale coordinates to Visium-scale so NicheDE's internal Rfast::dista
# scale computation (scale = sigma / min_dist) doesn't produce NA.
# Rfast::dista with trans=TRUE returns query×ref matrix; apply() over rows
# returns NA for some rows → median(mdist)=NA → scale=NA → all coords→NA.
# Fix: normalize so approx_min_dist ≈ sigma, giving scale ≈ 1.
{
  sample_idx   <- sample(nrow(coord_valid), min(500, nrow(coord_valid)))
  sample_D     <- as.matrix(dist(coord_valid[sample_idx, ]))
  approx_min   <- mean(apply(sample_D, 2, function(x) sort(x[x > 0])[1]))
  pre_scale    <- opt$sigma / approx_min
  coord_valid  <- coord_valid * pre_scale
  cat(sprintf("  Coord pre-scale: %.4fx (approx_min_dist=%.2f µm → target sigma=%g)\n",
              pre_scale, approx_min, opt$sigma))
}

# Build binary deconv_mat (cells × cell_types)
ct_levels   <- sort(unique(ct_valid))
deconv_mat  <- matrix(0, nrow=length(ct_valid), ncol=length(ct_levels),
                      dimnames=list(rownames(counts_valid), ct_levels))
for (i in seq_along(ct_valid)) {
  deconv_mat[i, ct_valid[i]] <- 1
}
cat(sprintf("  deconv_mat: %d cells × %d types\n", nrow(deconv_mat), ncol(deconv_mat)))
cat(sprintf("  Cell types: %s\n", paste(ct_levels, collapse=", ")))

# Library matrix: cell_types × genes reference expression (mean counts per type)
# CreateLibraryMatrix expects:
#   data      : cells × genes matrix
#   cell_type : 2-column matrix (col1=cell_name, col2=cell_type_label)
ct_mat <- cbind(rownames(counts_valid), ct_valid)
library_mat <- CreateLibraryMatrix(
  data      = counts_valid,
  cell_type = ct_mat
)
cat(sprintf("  library_mat: %d types × %d genes\n", nrow(library_mat), ncol(library_mat)))
# Reorder library_mat rows to match deconv_mat column order
library_mat <- library_mat[colnames(deconv_mat), , drop=FALSE]
cat(sprintf("  library_mat reordered to: %s\n", paste(rownames(library_mat), collapse=", ")))

# ─── Create NicheDE object ────────────────────────────────────────────────────
cat(sprintf("\n[3/5] Creating NicheDE object (sigma=%g)...\n", opt$sigma))

nde_obj <- CreateNicheDEObject(
  counts_mat    = counts_valid,
  coordinate_mat = coord_valid,
  library_mat   = library_mat,
  deconv_mat    = deconv_mat,
  sigma         = opt$sigma,
  Int           = TRUE
)

# NicheDE bug: for ≥10000 cells, Rfast::dista(trans=TRUE) returns a query×ref
# matrix; apply(D,1,...)[3] gives NA for some rows → median(mdist)=NA →
# scale=NA → object@coord becomes all-NA. Detect and repair here.
if (any(is.na(nde_obj@coord))) {
  cat("  [FIX] NicheDE coord NA detected (Rfast::dista bug) — patching with pre-scaled coords\n")
  nde_obj@coord <- coord_valid
  cat(sprintf("  [FIX] coord patched: %d × %d, range x=[%.1f,%.1f] y=[%.1f,%.1f]\n",
              nrow(nde_obj@coord), ncol(nde_obj@coord),
              min(nde_obj@coord[,1]), max(nde_obj@coord[,1]),
              min(nde_obj@coord[,2]), max(nde_obj@coord[,2])))
} else {
  cat(sprintf("  coord OK: range x=[%.1f,%.1f] y=[%.1f,%.1f]\n",
              min(nde_obj@coord[,1]), max(nde_obj@coord[,1]),
              min(nde_obj@coord[,2]), max(nde_obj@coord[,2])))
}

# Calculate effective niche — use LargeScale version for 268k cells
# batch_size controls memory: 5000 processes ~5000 cells at a time
cat("  Calculating effective niche (LargeScale, batch_size=5000)...\n")
nde_obj <- CalculateEffectiveNicheLargeScale(nde_obj, batch_size=5000, cutoff=0.05)
cat("  Effective niche calculated.\n")

# ─── Run NicheDE ──────────────────────────────────────────────────────────────
cat(sprintf("\n[4/5] Running NicheDE (cores=%d, C=%d, M=%d)...\n",
            opt$num_cores, opt$C, opt$M))

# NicheDE parallel bug: SOCK workers crash when serializing 620 MB dense
# counts_chunk. Fix: register doMC (FORK-based, no data serialization)
# before niche_DE runs, so foreach uses mclapply under the hood.
# doMC overrides the doParallel registration that niche_DE would do.
tryCatch({
  library(doMC)
  doMC::registerDoMC(cores = opt$num_cores)
  cat(sprintf("  [FIX] Registered doMC FORK backend (%d cores) — no SOCK serialization\n",
              opt$num_cores))
}, error = function(e) {
  cat(sprintf("  [FIX] doMC unavailable (%s) — NicheDE will use SOCK (may fail on large data)\n",
              e$message))
})

outfile <- file.path(opt$outdir, "nichede_raw_output.csv")
nde_obj <- niche_DE(
  object    = nde_obj,
  num_cores = opt$num_cores,
  outfile   = outfile,
  C         = opt$C,
  M         = opt$M,
  gamma     = 0.8,
  print     = TRUE,
  Int       = TRUE,
  batch     = TRUE,
  self_EN   = FALSE,
  G         = 1
)

cat("  NicheDE complete.\n")

# ─── Extract and save results ─────────────────────────────────────────────────
cat("\n[5/5] Extracting and saving results...\n")

# get_niche_DE_genes returns significant niche-DE gene list
# niche_DE_markers returns a summary table
tryCatch({
  markers <- niche_DE_markers(nde_obj)
  write.csv(markers, file.path(opt$outdir, "nichede_markers.csv"), row.names=FALSE)
  cat(sprintf("  Markers saved: %d entries\n", nrow(markers)))
}, error = function(e) {
  cat(sprintf("  Warning: niche_DE_markers failed: %s\n", e$message))
})

# Extract p-values and LFC per cell type × niche combination
tryCatch({
  pvals <- get_niche_DE_pval_raw(nde_obj)
  # pvals is a list of matrices (genes × cell_types) indexed by niche
  # Flatten into a tidy data frame
  result_rows <- list()
  for (niche_name in names(pvals)) {
    p_mat <- pvals[[niche_name]]
    if (is.null(p_mat)) next
    for (ct_name in colnames(p_mat)) {
      p_vec <- p_mat[, ct_name]
      genes_sig <- names(p_vec)[!is.na(p_vec) & p_vec < 0.05]
      if (length(genes_sig) > 0) {
        result_rows[[length(result_rows)+1]] <- data.frame(
          niche     = niche_name,
          cell_type = ct_name,
          gene      = genes_sig,
          pval      = p_vec[genes_sig],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(result_rows) > 0) {
    results_df <- do.call(rbind, result_rows)
    # BH correction across all tests
    results_df$padj <- p.adjust(results_df$pval, method="BH")
    results_df <- results_df[order(results_df$padj), ]
    write.csv(results_df, file.path(opt$outdir, "nichede_significant_genes.csv"),
              row.names=FALSE)
    cat(sprintf("  Significant niche-DE genes (p<0.05): %d\n", nrow(results_df)))
    cat(sprintf("  Significant niche-DE genes (padj<0.05): %d\n",
                sum(results_df$padj < 0.05, na.rm=TRUE)))
  } else {
    cat("  No significant niche-DE genes found at p<0.05\n")
    write.csv(data.frame(), file.path(opt$outdir, "nichede_significant_genes.csv"),
              row.names=FALSE)
  }
}, error = function(e) {
  cat(sprintf("  Warning: p-value extraction failed: %s\n", e$message))
  # Fall back to raw output CSV if it exists
  if (file.exists(outfile)) {
    raw <- read.csv(outfile)
    write.csv(raw, file.path(opt$outdir, "nichede_significant_genes.csv"), row.names=FALSE)
    cat(sprintf("  Saved raw output: %d rows\n", nrow(raw)))
  }
})

# Save cell type × domain summary
summary_df <- data.frame(
  cell_type  = ct_levels,
  n_cells    = as.integer(table(ct_valid)[ct_levels])
)
domain_ct_tab <- table(domain_valid, ct_valid)
write.csv(as.data.frame.matrix(domain_ct_tab),
          file.path(opt$outdir, "domain_celltype_crosstab.csv"))
cat(sprintf("  Domain × cell type crosstab saved.\n"))

elapsed <- as.numeric(difftime(Sys.time(), t0, units="secs"))
write_json(list(
  n_cells       = sum(valid_mask),
  n_genes       = ncol(counts_valid),
  n_cell_types  = length(ct_levels),
  sigma         = opt$sigma,
  C             = opt$C,
  M             = opt$M,
  execution_seconds = round(elapsed, 1)
), file.path(opt$outdir, "nichede_summary.json"), pretty=TRUE, auto_unbox=TRUE)

cat(sprintf("\n✓ NicheDE complete in %.1fs\n", elapsed))
