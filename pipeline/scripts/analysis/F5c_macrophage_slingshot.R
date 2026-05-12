#!/usr/bin/env Rscript
# F5c — Macrophage M1↔M2 Slingshot trajectory (audit-fixed)
#
# Replaces F5_pseudotime_slingshot.R for the macrophage sub-trajectory.
# Differences from F5:
#   1. Trajectory restricted to a biologically valid axis (M1↔M2 polarization)
#   2. tradeSeq fit on REAL gene counts (not embedding dims)
#   3. Bonferroni correction explicitly applied
#
# Inputs (from F5b):
#   - macrophage_counts_for_R.csv       (n_cells × 289 genes raw counts)
#   - macrophage_metadata_for_R.csv     (cell_type_L2, leiden, dpt, umap)
#   - macrophage_pca_for_R.csv          (30 PCs)

suppressPackageStartupMessages({
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
.this_file <- normalizePath(sub("^--file=", "",
                                grep("^--file=", commandArgs(trailingOnly = FALSE),
                                     value = TRUE)[1]),
                            mustWork = FALSE)
source(file.path(dirname(.this_file), "..", "utils", "paths.R"))
BASE  <- xenium_dataset_root()
INDIR <- file.path(BASE, "results/03_phase3_spatial/F5b_macrophage_pseudotime")
OUT   <- INDIR
FIG   <- file.path(BASE, "results/figures/phase5_pseudotime_v2")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

log_msg <- function(m) cat(format(Sys.time(), "[%H:%M:%S]"), m, "\n")

# ── 1. Load ───────────────────────────────────────────────────────────────────
log_msg("Loading exports from F5b ...")
counts <- read.csv(file.path(INDIR, "macrophage_counts_for_R.csv"), row.names=1, check.names=FALSE)
meta   <- read.csv(file.path(INDIR, "macrophage_metadata_for_R.csv"), row.names=1)
pca    <- read.csv(file.path(INDIR, "macrophage_pca_for_R.csv"), row.names=1)

stopifnot(all(rownames(counts) == rownames(meta)))
stopifnot(all(rownames(pca) == rownames(meta)))
log_msg(sprintf("  counts: %d cells × %d genes", nrow(counts), ncol(counts)))
log_msg(sprintf("  PCA dims: %d", ncol(pca)))

# ── 2. Build SCE ──────────────────────────────────────────────────────────────
log_msg("Building SingleCellExperiment ...")
sce <- SingleCellExperiment(
  assays = list(counts = t(as.matrix(counts))),
  colData = meta
)
reducedDim(sce, "PCA") <- as.matrix(pca)

# ── 3. Slingshot ──────────────────────────────────────────────────────────────
log_msg("Running Slingshot — start = Macrophage_M1 ...")
set.seed(42)
sce <- slingshot(
  sce,
  clusterLabels = colData(sce)$cell_type_L2,
  reducedDim    = "PCA",
  start.clus    = "Macrophage_M1",
  approx_points = 150
)

pt_cols <- grep("slingPseudotime", names(colData(sce)), value=TRUE)
log_msg(sprintf("  Slingshot lineages: %d", length(pt_cols)))

# Save pseudotime
pt_df <- as.data.frame(colData(sce)[, c("cell_type_L2", "macro_leiden", "dpt_pseudotime",
                                         "umap_1", "umap_2", pt_cols)])
write.csv(pt_df, file.path(OUT, "macrophage_slingshot_pseudotime.csv"))
log_msg("  Saved macrophage_slingshot_pseudotime.csv")

# ── 4. Concordance with DPT ───────────────────────────────────────────────────
log_msg("Computing DPT vs Slingshot concordance ...")
sl_pt <- pt_df[[pt_cols[1]]]
dpt   <- pt_df$dpt_pseudotime
keep  <- !is.na(sl_pt) & !is.na(dpt)
spearman_r <- cor(sl_pt[keep], dpt[keep], method="spearman")
pearson_r  <- cor(sl_pt[keep], dpt[keep], method="pearson")
log_msg(sprintf("  Spearman r (Slingshot vs DPT): %.3f", spearman_r))
log_msg(sprintf("  Pearson r:                     %.3f", pearson_r))

# Concordance figure
df_conc <- data.frame(slingshot=sl_pt[keep], dpt=dpt[keep],
                      celltype=pt_df$cell_type_L2[keep])
p_conc <- ggplot(df_conc, aes(slingshot, dpt, color=celltype)) +
  geom_point(size=0.3, alpha=0.4) +
  scale_color_manual(values=c(Macrophage_M1="#1f77b4", Macrophage_M2="#d62728")) +
  geom_smooth(method="loess", color="black", se=FALSE, size=0.7) +
  theme_minimal(base_size=11) +
  labs(x="Slingshot pseudotime (Lineage 1)", y="DPT pseudotime",
       title=sprintf("Pseudotime concordance — Spearman r = %.3f", spearman_r),
       color="Cell type")
ggsave(file.path(FIG, "F5c_concordance_slingshot_vs_dpt.png"), p_conc,
       width=7, height=5, dpi=200)
log_msg("  Saved concordance figure")

# ── 5. tradeSeq with REAL gene counts (audit fix) ────────────────────────────
# tradeSeq is statistically appropriate ONLY for count data. We use the actual
# 289-gene Xenium counts here (not embedding dimensions as F5 did).
log_msg("Fitting tradeSeq GAM on real gene counts ...")
keep_gam <- !is.na(sl_pt)
sl_pt_filt <- sl_pt[keep_gam]
counts_gam <- t(as.matrix(counts[keep_gam, ]))

# Filter low-count genes (tradeSeq defaults)
gene_filter <- rowSums(counts_gam) > 50
counts_gam_f <- counts_gam[gene_filter, ]
log_msg(sprintf("  Genes after count filter (sum>50): %d / %d",
                nrow(counts_gam_f), nrow(counts_gam)))

cellWeights <- matrix(1, nrow=ncol(counts_gam_f), ncol=1)
pseudotime  <- matrix(sl_pt_filt, nrow=ncol(counts_gam_f), ncol=1)

set.seed(42)
result <- tryCatch({
  sce_gam <- fitGAM(
    counts = counts_gam_f,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    nknots = 6,
    verbose = FALSE
  )

  assoc_res <- associationTest(sce_gam)
  assoc_res$gene <- rownames(assoc_res)

  # Bonferroni correction (audit fix — was missing in F5)
  n_tests <- nrow(assoc_res)
  assoc_res$pval_bonferroni <- pmin(assoc_res$pvalue * n_tests, 1)

  # BH for comparison
  assoc_res$pval_BH_FDR <- p.adjust(assoc_res$pvalue, method="BH")

  assoc_sig <- subset(assoc_res, pval_BH_FDR < 0.05)
  assoc_sig_bonf <- subset(assoc_res, pval_bonferroni < 0.05)
  log_msg(sprintf("  Genes pseudotime-associated (BH FDR < 0.05): %d", nrow(assoc_sig)))
  log_msg(sprintf("  Genes pseudotime-associated (Bonferroni  < 0.05): %d", nrow(assoc_sig_bonf)))

  write.csv(assoc_res[order(assoc_res$pvalue), ],
            file.path(OUT, "macrophage_tradeseq_gene_assoc.csv"))
  log_msg("  Saved tradeseq gene association table")

  # Top 6 genes — smoother plots
  top_genes <- rownames(assoc_res)[order(assoc_res$pvalue)][1:6]
  for (g in top_genes) {
    png(file.path(FIG, paste0("F5c_gam_", g, ".png")),
        width=800, height=600, res=120)
    plotSmoothers(sce_gam, counts_gam_f, gene=g, border=TRUE)
    title(main=paste("tradeSeq GAM —", g, "across M1↔M2 trajectory"))
    dev.off()
  }
  log_msg("  Saved top-6 gene smoother plots")

  list(success = TRUE, n_sig_bh = nrow(assoc_sig), n_sig_bonf = nrow(assoc_sig_bonf),
       top_genes = top_genes)
}, error = function(e) {
  log_msg(paste("tradeSeq error:", conditionMessage(e)))
  list(success = FALSE, error = conditionMessage(e))
})

# ── 6. Summary ────────────────────────────────────────────────────────────────
summary_obj <- list(
  date              = format(Sys.time(), "%Y-%m-%d"),
  n_cells           = nrow(meta),
  n_M1              = sum(meta$cell_type_L2 == "Macrophage_M1"),
  n_M2              = sum(meta$cell_type_L2 == "Macrophage_M2"),
  n_lineages        = length(pt_cols),
  start_cluster     = "Macrophage_M1",
  spearman_dpt_vs_slingshot = round(spearman_r, 3),
  pearson_dpt_vs_slingshot  = round(pearson_r, 3),
  tradeseq_results  = result
)
writeLines(toJSON(summary_obj, pretty=TRUE, auto_unbox=TRUE),
           file.path(OUT, "F5c_summary.json"))

log_msg("=== F5c COMPLETE ===")
log_msg(paste("Outputs:", OUT))
log_msg(paste("Figures:", FIG))
