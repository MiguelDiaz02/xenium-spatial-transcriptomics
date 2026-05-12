#!/usr/bin/env Rscript
# M1 — Pseudobulk DGE via muscat::aggregateData + edgeR/DESeq2/limma
# ===================================================================
# R counterpart to pseudobulk_decoupler.py. muscat is the methodologically
# rigorous reference implementation for spatial-/single-cell pseudobulk DE
# (Crowell et al., Nat Commun 2020 — formal pseudo-replication treatment).
#
# Run AFTER concat_samples.py has produced cohort.h5ad. This script reads
# the h5ad via zellkonverter, builds a SingleCellExperiment, then runs
# muscat::pbDS with three methods: edgeR-LRT, DESeq2-Wald, limma-voom.
#
# Outputs:
#   TBDs/cohort/results/pseudobulk_muscat/<organ>/<celltype>/
#       <contrast>__edgeR.tsv
#       <contrast>__DESeq2.tsv
#       <contrast>__limma.tsv
#
# Same columns as decoupler-py output (gene, log2FoldChange, lfcSE, stat,
# pvalue, padj, baseMean) so pseudobulk_consensus.py can join them.
#
# Args (via commandArgs(trailingOnly=TRUE), in order):
#   --input  <cohort.h5ad>
#   --cohort <cohort_TBDs.yaml-json or cohort_meta.json>
#   --outdir <output root>
#   --celltype_col cell_type_L2
#   --organ lung|liver|all
#
# Environment: conda env xenium_R_analysis (R ≥ 4.4, muscat ≥ 1.20,
# zellkonverter ≥ 1.16, SingleCellExperiment, edgeR, DESeq2, limma).

suppressPackageStartupMessages({
  library(jsonlite)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(muscat)
  library(edgeR)
  library(DESeq2)
  library(limma)
})

# ──────────────────────────────────────────────────────────────────────
# CLI parsing (lightweight; full optparse would add a dep)
# ──────────────────────────────────────────────────────────────────────
parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input = NA_character_,
    cohort = NA_character_,
    outdir = NA_character_,
    celltype_col = "cell_type_L2",
    organ = "all"
  )
  i <- 1
  while (i <= length(argv)) {
    k <- sub("^--", "", argv[i])
    v <- argv[i + 1]
    if (!is.null(defaults[[k]])) defaults[[k]] <- v
    i <- i + 2
  }
  defaults
}

# ──────────────────────────────────────────────────────────────────────
# Build SCE from h5ad
# ──────────────────────────────────────────────────────────────────────
load_cohort_sce <- function(h5ad_path) {
  message(sprintf("[muscat] reading %s", h5ad_path))
  sce <- readH5AD(h5ad_path, raw = FALSE, layers = "counts")

  # muscat expects assay named 'counts'
  if (!"counts" %in% assayNames(sce)) {
    if ("X" %in% assayNames(sce)) {
      message("[muscat] no 'counts' assay → using .X (assume log1p inverted to counts upstream)")
      assayNames(sce)[assayNames(sce) == "X"] <- "counts"
    } else {
      stop("h5ad lacks both 'counts' and 'X' assays")
    }
  }
  sce
}

# ──────────────────────────────────────────────────────────────────────
# Run one organ × contrast × method
# ──────────────────────────────────────────────────────────────────────
run_method <- function(pb, contrast_name, ref_lvl, test_lvl, method) {
  res <- tryCatch(
    pbDS(
      pb,
      method = method,
      design = ~ condition,
      coef = NULL,
      contrast = limma::makeContrasts(
        contrasts = sprintf("condition%s - condition%s", test_lvl, ref_lvl),
        levels = levels(droplevels(pb$condition))
      ),
      verbose = FALSE
    ),
    error = function(e) { message(sprintf("  [%s] %s", method, e$message)); NULL }
  )
  if (is.null(res)) return(NULL)

  out <- do.call(rbind, lapply(names(res$table[[1]]), function(ct) {
    tbl <- res$table[[1]][[ct]]
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    data.frame(
      gene = tbl$gene,
      log2FoldChange = tbl$logFC,
      lfcSE = if ("logCPM" %in% names(tbl)) NA_real_ else NA_real_,
      stat = if ("F" %in% names(tbl)) tbl$F else NA_real_,
      pvalue = tbl$p_val,
      padj = tbl$p_adj.loc,
      baseMean = NA_real_,
      celltype = ct,
      contrast = contrast_name,
      method = method,
      stringsAsFactors = FALSE
    )
  }))
  out
}

write_results <- function(df, outdir, organ) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  for (ct in unique(df$celltype)) {
    for (m in unique(df$method)) {
      sub <- df[df$celltype == ct & df$method == m, , drop = FALSE]
      if (nrow(sub) == 0) next
      d <- file.path(outdir, organ, gsub("[ /]", "_", ct))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      out <- file.path(d, sprintf("%s__%s.tsv", unique(sub$contrast), m))
      write.table(sub, out, sep = "\t", row.names = FALSE, quote = FALSE)
      message(sprintf("  wrote %s (%d rows)", out, nrow(sub)))
    }
  }
}

# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
main <- function() {
  args <- parse_args()
  stopifnot(file.exists(args$input))
  stopifnot(file.exists(args$cohort))

  cohort <- fromJSON(args$cohort, simplifyVector = FALSE)
  organs <- if (args$organ == "all") names(cohort$organs) else args$organ

  sce <- load_cohort_sce(args$input)
  message(sprintf("[muscat] SCE: %d cells × %d genes", ncol(sce), nrow(sce)))

  # muscat expects three columns: cluster_id (cell type), sample_id (donor),
  # group_id (condition).
  sce$cluster_id <- factor(colData(sce)[[args$celltype_col]])
  sce$sample_id <- factor(sce$subject_id)
  sce$group_id  <- factor(sce$condition, levels = unlist(cohort$conditions))
  sce <- prepSCE(sce, kid = "cluster_id", sid = "sample_id", gid = "group_id",
                 drop = TRUE)

  for (org in organs) {
    org_cells <- sce$organ == org
    if (sum(org_cells) == 0) {
      message(sprintf("[%s] no cells — skip", org))
      next
    }
    org_sce <- sce[, org_cells]
    org_sce$sample_id <- droplevels(org_sce$sample_id)
    org_sce$group_id  <- droplevels(org_sce$group_id)
    org_sce$cluster_id <- droplevels(org_sce$cluster_id)

    message(sprintf("[%s] %d cells, %d cell types, %d donors",
                    org, ncol(org_sce),
                    nlevels(org_sce$cluster_id),
                    nlevels(org_sce$sample_id)))

    pb <- aggregateData(org_sce, assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))

    contrasts <- cohort$organs[[org]]$contrasts
    for (cc in contrasts) {
      cname <- sprintf("%s_vs_%s", cc$test, cc$ref)
      message(sprintf("[%s] contrast %s", org, cname))
      for (m in c("edgeR", "DESeq2", "limma-voom")) {
        msafe <- gsub("-", "_", m)
        df <- run_method(pb, cname, cc$ref, cc$test, m)
        if (!is.null(df)) df$method <- msafe
        write_results(df, args$outdir, org)
      }
    }
  }
  message("[muscat] done")
}

main()
