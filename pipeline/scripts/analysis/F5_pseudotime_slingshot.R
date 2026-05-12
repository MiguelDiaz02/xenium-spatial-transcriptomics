#!/usr/bin/env Rscript
# F5_pseudotime_slingshot.R  ── DEPRECATED 2026-05-06 ──
#
# DO NOT USE. Replaced by F5b_macrophage_pseudotime.py + F5c_macrophage_slingshot.R.
#
# Two methodological problems identified in audit (2026-05-06):
#   (1) Slingshot was given start.clus = Tumor_proliferating with end.clus =
#       (T cells, AT1, endothelial). Tumor cells do not differentiate into
#       these targets — the inferred "trajectory" is a spatial gradient
#       relabelled as developmental pseudotime.
#   (2) tradeSeq fitGAM was applied to Novae embedding dimensions scaled and
#       rounded to "counts". tradeSeq's negative-binomial GAM assumes count
#       data with a specific mean-variance relationship; embedding dims do
#       not satisfy this and the resulting p-values are statistically
#       invalid.
#
# The replacement (F5b + F5c) restricts pseudotime to a biologically valid
# axis (Macrophage M1↔M2 polarization) and applies tradeSeq to real Xenium
# gene counts with explicit Bonferroni correction.
#
# Original input: cells_for_slingshot.csv (Novae 64d embeddings + metadata)
# Original output: F5_pseudotime/

suppressPackageStartupMessages({
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(scater)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(viridis)
  library(jsonlite)
})

# ── Paths ──────────────────────────────────────────────────────────────────────
.this_file <- normalizePath(sub("^--file=", "",
                                grep("^--file=", commandArgs(trailingOnly = FALSE),
                                     value = TRUE)[1]),
                            mustWork = FALSE)
source(file.path(dirname(.this_file), "..", "utils", "paths.R"))
BASE  <- xenium_dataset_root()
INPUT <- file.path(BASE, "results/03_phase3_spatial/F5_pseudotime/cells_for_slingshot.csv")
EXPR  <- file.path(BASE, "results/03_phase3_spatial/F5_pseudotime")  # same dir
OUT   <- file.path(BASE, "results/03_phase3_spatial/F5_pseudotime")
FIG   <- file.path(BASE, "results/figures/phase5_pseudotime")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(format(Sys.time(), "[%H:%M:%S]"), msg, "\n", flush=TRUE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
log_msg("Loading cell data...")
cells <- read.csv(INPUT, row.names=1)
log_msg(paste("Loaded", nrow(cells), "cells,", ncol(cells), "columns"))

# Novae embedding cols
embed_cols <- grep("^novae_dim_", names(cells), value=TRUE)
log_msg(paste("Using", length(embed_cols), "Novae embedding dims"))

# ── 2. Focus subset: TME-relevant types ───────────────────────────────────────
# Pseudotime along tumor → immune → stromal gradient
# Exclude: Unassigned, rare (<100 cells), pure structural types
KEEP_TYPES <- c(
  "Tumor_proliferating", "Tumor_resting",
  "Macrophage_M1", "Macrophage_M2", "Monocyte_classical",
  "CD8_T_cytotoxic", "CD8_T_exhausted", "CD4_T_helper", "Treg",
  "NK_cytotoxic", "cDC1", "cDC2",
  "Endothelial_blood", "Endothelial_lymphatic",
  "Epithelial_general", "AT1", "Ciliated"
)

cells_filt <- cells[cells$cell_type_L2 %in% KEEP_TYPES, ]
log_msg(paste("After type filtering:", nrow(cells_filt), "cells"))

# Subsample for tractability (Slingshot scales O(n^2))
set.seed(42)
MAX_CELLS <- 30000
if (nrow(cells_filt) > MAX_CELLS) {
  # Stratified subsample by cell type
  idx <- cells_filt %>%
    mutate(.row = row_number()) %>%
    group_by(cell_type_L2) %>%
    slice_sample(prop = MAX_CELLS / nrow(cells_filt), replace=FALSE) %>%
    pull(.row)
  cells_sub <- cells_filt[idx, ]
} else {
  cells_sub <- cells_filt
}
log_msg(paste("Working with", nrow(cells_sub), "cells (stratified subsample)"))
log_msg(paste("Type distribution:"))
print(table(cells_sub$cell_type_L2))

# ── 3. Build SCE with Novae embeddings as reducedDim ─────────────────────────
log_msg("Building SingleCellExperiment...")
embed_mat <- as.matrix(cells_sub[, embed_cols])
rownames(embed_mat) <- rownames(cells_sub)

sce <- SingleCellExperiment(
  colData = cells_sub[, c("cell_type_L1","cell_type_L2","cell_type_L3","x","y")]
)
reducedDim(sce, "NovaeEmbed") <- embed_mat

# UMAP on Novae embeddings for visualization
log_msg("Computing UMAP for visualization...")
set.seed(42)
sce <- runUMAP(sce, dimred="NovaeEmbed", name="UMAP", n_neighbors=30, min_dist=0.3)

# ── 4. Slingshot trajectory ───────────────────────────────────────────────────
# Start cluster: Tumor (proliferating = most undifferentiated)
# End clusters: Terminal immune/stromal states
log_msg("Running Slingshot...")
START_CLUSTER <- "Tumor_proliferating"
END_CLUSTERS  <- c("CD8_T_exhausted", "Treg", "Endothelial_blood", "AT1", "Ciliated")

set.seed(42)
sce <- slingshot(
  sce,
  clusterLabels = colData(sce)$cell_type_L2,
  reducedDim    = "NovaeEmbed",
  start.clus    = START_CLUSTER,
  end.clus      = END_CLUSTERS,
  approx_points = 150
)

# Extract pseudotime
pt_cols  <- grep("slingPseudotime", names(colData(sce)), value=TRUE)
n_lineages <- length(pt_cols)
log_msg(paste("Slingshot found", n_lineages, "lineages:", paste(pt_cols, collapse=", ")))

# Save pseudotime table
pt_df <- as.data.frame(colData(sce)[, c("cell_type_L2","x","y", pt_cols)])
write.csv(pt_df, file.path(OUT, "pseudotime_per_cell.csv"))
log_msg("Saved pseudotime_per_cell.csv")

# ── 5. Spatial gradient: distance to tumor centroid ───────────────────────────
log_msg("Computing spatial distance to tumor centroid...")
tumor_cells <- cells_sub[cells_sub$cell_type_L2 %in% c("Tumor_proliferating","Tumor_resting"), ]
tumor_centroid <- c(mean(tumor_cells$x), mean(tumor_cells$y))

pt_df$dist_to_tumor <- sqrt(
  (cells_sub$x - tumor_centroid[1])^2 +
  (cells_sub$y - tumor_centroid[2])^2
)

write.csv(pt_df, file.path(OUT, "pseudotime_with_spatial.csv"))
log_msg("Saved pseudotime_with_spatial.csv")

# ── 6. Figures ────────────────────────────────────────────────────────────────
log_msg("Generating figures...")

# Color palette for cell types
type_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(KEEP_TYPES)),
  KEEP_TYPES
)

# Fig 1: UMAP colored by cell type
umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP1","UMAP2")
umap_df$cell_type_L2 <- colData(sce)$cell_type_L2
# Add first lineage pseudotime if available
if (length(pt_cols) > 0) {
  umap_df$pseudotime <- colData(sce)[[pt_cols[1]]]
}

p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color=cell_type_L2)) +
  geom_point(size=0.3, alpha=0.5) +
  scale_color_manual(values=type_colors, na.value="grey80") +
  theme_minimal(base_size=12) +
  guides(color=guide_legend(override.aes=list(size=3), ncol=2)) +
  labs(title="UMAP — Novae embeddings (TME subset)", color="Cell type L2")

ggsave(file.path(FIG, "F5_umap_celltypes.png"), p1, width=12, height=8, dpi=200)
log_msg("Saved F5_umap_celltypes.png")

# Fig 2: UMAP colored by pseudotime (lineage 1)
if ("pseudotime" %in% names(umap_df)) {
  p2 <- ggplot(umap_df, aes(UMAP1, UMAP2, color=pseudotime)) +
    geom_point(size=0.3, alpha=0.5) +
    scale_color_viridis(option="plasma", na.value="grey90") +
    theme_minimal(base_size=12) +
    labs(title="Pseudotime — Lineage 1 (Tumor → terminal states)", color="Pseudotime")

  ggsave(file.path(FIG, "F5_umap_pseudotime_L1.png"), p2, width=10, height=8, dpi=200)
  log_msg("Saved F5_umap_pseudotime_L1.png")
}

# Fig 3: Spatial map of pseudotime
p3 <- ggplot(pt_df, aes(x, y, color=.data[[pt_cols[1]]])) +
  geom_point(size=0.2, alpha=0.4) +
  scale_color_viridis(option="plasma", na.value="grey90") +
  coord_equal() +
  theme_void(base_size=12) +
  labs(title="Spatial pseudotime (Lineage 1)", color="Pseudotime")

ggsave(file.path(FIG, "F5_spatial_pseudotime.png"), p3, width=10, height=9, dpi=200)
log_msg("Saved F5_spatial_pseudotime.png")

# Fig 4: Pseudotime vs distance to tumor centroid
if (length(pt_cols) > 0) {
  p4 <- ggplot(pt_df, aes(dist_to_tumor, .data[[pt_cols[1]]], color=cell_type_L2)) +
    geom_point(size=0.3, alpha=0.3) +
    geom_smooth(method="loess", color="black", se=TRUE, size=1) +
    scale_color_manual(values=type_colors, na.value="grey80") +
    theme_minimal(base_size=12) +
    labs(
      x="Distance to tumor centroid (µm)",
      y="Pseudotime (Lineage 1)",
      title="Pseudotime gradient: tumor core → periphery",
      color="Cell type L2"
    ) +
    guides(color=guide_legend(override.aes=list(size=3), ncol=2))

  ggsave(file.path(FIG, "F5_pseudotime_vs_spatial_gradient.png"), p4, width=12, height=8, dpi=200)
  log_msg("Saved F5_pseudotime_vs_spatial_gradient.png")
}

# ── 7. tradeSeq — GAM on top genes ────────────────────────────────────────────
# We use Novae embedding dims as proxy "expression" for GAM fitting
# This avoids needing the full count matrix in R and is valid since Novae
# dims capture gene expression variation
log_msg("Fitting tradeSeq GAMs on Novae embedding dims (proxy expression)...")

# Use first lineage only for GAM (most cells assigned)
pt_lineage1 <- colData(sce)[[pt_cols[1]]]
keep_gam    <- !is.na(pt_lineage1)
log_msg(paste("Cells with Lineage 1 pseudotime:", sum(keep_gam)))

if (sum(keep_gam) > 500) {
  # Build count-like matrix from Novae dims (scaled to [0,100])
  embed_gam <- t(embed_mat[keep_gam, ])
  # Normalize dims to positive integers (GAM needs counts-like input)
  embed_min  <- apply(embed_gam, 1, min)
  embed_gam_pos <- embed_gam - embed_min + 1
  embed_gam_int <- round(embed_gam_pos * 10)  # scale to reasonable counts

  cellWeights <- matrix(1, nrow=sum(keep_gam), ncol=1)
  pseudotime  <- matrix(pt_lineage1[keep_gam], nrow=sum(keep_gam), ncol=1)

  log_msg("Fitting GAMs (this takes ~5-15 min)...")
  set.seed(42)
  tryCatch({
    sce_gam <- fitGAM(
      counts      = embed_gam_int,
      pseudotime  = pseudotime,
      cellWeights = cellWeights,
      nknots      = 6,
      verbose     = FALSE
    )

    # Association test: which dims change with pseudotime
    assoc_res <- associationTest(sce_gam)
    assoc_res$dim <- rownames(assoc_res)
    write.csv(assoc_res[order(assoc_res$pvalue),],
              file.path(OUT, "tradeseq_association_dims.csv"))
    log_msg("Saved tradeseq_association_dims.csv")

    # Plot top 4 variable dims along pseudotime
    top_dims <- rownames(assoc_res)[order(assoc_res$pvalue)][1:4]

    for (dim_name in top_dims) {
      png(file.path(FIG, paste0("F5_gam_", dim_name, ".png")),
          width=800, height=600, res=120)
      plotSmoothers(sce_gam, dim_name, border=TRUE)
      title(main=paste("tradeSeq GAM —", dim_name, "along Lineage 1"))
      dev.off()
    }
    log_msg("Saved GAM smoother plots")

  }, error=function(e) {
    log_msg(paste("tradeSeq GAM error (non-fatal):", conditionMessage(e)))
  })
} else {
  log_msg("Not enough cells for GAM — skipping tradeSeq")
}

# ── 8. Summary JSON ───────────────────────────────────────────────────────────
log_msg("Writing summary...")
type_counts <- as.list(table(cells_sub$cell_type_L2))

summary_obj <- list(
  date           = format(Sys.time(), "%Y-%m-%d"),
  n_cells_input  = nrow(cells),
  n_cells_subset = nrow(cells_sub),
  keep_types     = KEEP_TYPES,
  n_lineages     = n_lineages,
  lineage_names  = pt_cols,
  start_cluster  = START_CLUSTER,
  end_clusters   = END_CLUSTERS,
  tumor_centroid = list(x=tumor_centroid[1], y=tumor_centroid[2]),
  type_counts    = type_counts
)
writeLines(toJSON(summary_obj, pretty=TRUE, auto_unbox=TRUE),
           file.path(OUT, "F5_summary.json"))

log_msg("=== Fase E COMPLETE ===")
log_msg(paste("Outputs in:", OUT))
log_msg(paste("Figures in:", FIG))
