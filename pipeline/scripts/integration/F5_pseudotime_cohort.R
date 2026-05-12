#!/usr/bin/env Rscript
# M4 — Cohort-level pseudotime: Slingshot + tradeSeq::conditionTest
# ===================================================================
# Status: SCAFFOLD (dry-run executable; full computation requires
# xenium_R_analysis env + integrated h5ad on disk).
#
# Cohort pseudotime is fundamentally different from the per-sample pilot:
# pooling all donors with diverse fibrosis severity blurs trajectories.
# We use the integrated scVI latent space (computed in M1) to align cells
# across donors before fitting Slingshot. Then per-cell-type subsetting
# isolates the lineage of interest (e.g. AT2→AT1, HSC quiescent→activated).
#
# Algorithm (per organ × target lineage):
#   1. Read cohort_integrated.h5ad (zellkonverter).
#   2. Subset to organ + anchor_celltypes (from cohort.pseudotime_cohort).
#   3. Use obsm$X_scVI as the integrated low-dim representation.
#   4. Cluster with HDBSCAN or Leiden at low resolution (anchors define
#      start/end clusters).
#   5. slingshot() on X_scVI with explicit start/end clusters.
#   6. tradeSeq::fitGAM() with `condition` as a factor.
#   7. tradeSeq::conditionTest() to identify genes whose pseudotime
#      profile differs between TBD vs control (or TBD vs Alc_Cirh).
#   8. tradeSeq::patternTest() for monotonic up/down patterns.
#
# Outputs:
#   TBDs/cohort/results/pseudotime_cohort/<organ>/<lineage>/
#       slingshot_curves.rds
#       tradeseq_fit.rds
#       conditionTest__<test>_vs_<ref>.tsv
#       patternTest.tsv
#       figures/curves_<lineage>.png
#       figures/condition_test_top20.png
#
# CLI:
#   Rscript F5_pseudotime_cohort.R
#       --input   TBDs/cohort/results/cohort_integrated.h5ad
#       --cohort  TBDs/cohort/results/cohort_meta.json
#       --organ   lung|liver|all
#       --lineage <name from cohort.pseudotime_cohort.target_lineages>
#       --execute     (else dry-run prints what would run)
#       --outdir  TBDs/cohort/results/pseudotime_cohort

suppressPackageStartupMessages({
  library(jsonlite)
})

parse_args <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    input = NA_character_,
    cohort = NA_character_,
    organ = "all",
    lineage = "all",
    outdir = NA_character_,
    execute = FALSE
  )
  i <- 1
  while (i <= length(argv)) {
    k <- sub("^--", "", argv[i])
    if (k == "execute") {
      defaults[[k]] <- TRUE
      i <- i + 1
      next
    }
    v <- argv[i + 1]
    if (!is.null(defaults[[k]])) defaults[[k]] <- v
    i <- i + 2
  }
  defaults
}

plan_lineages <- function(cohort, organ_filter) {
  lineages <- cohort$pseudotime_cohort$target_lineages
  out <- list()
  for (org in names(lineages)) {
    if (organ_filter != "all" && org != organ_filter) next
    for (l in lineages[[org]]) {
      out[[length(out) + 1]] <- list(
        organ = org,
        name = l$name,
        anchor_celltypes = unlist(l$anchor_celltypes)
      )
    }
  }
  out
}

dry_run <- function(plan, outdir) {
  message(sprintf("[pseudotime_cohort] planned %d lineages", length(plan)))
  for (p in plan) {
    message(sprintf("  - %s/%s (anchors: %s)",
                    p$organ, p$name,
                    paste(p$anchor_celltypes, collapse = ", ")))
  }
  if (!is.na(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    writeLines(toJSON(plan, pretty = TRUE),
               file.path(outdir, "plan.json"))
  }
}

execute_lineage <- function(adata_h5ad, plan_item, cohort, outdir) {
  # ─── REAL IMPLEMENTATION (sketched; uncomment when env is ready) ─────
  # suppressPackageStartupMessages({
  #   library(zellkonverter)
  #   library(SingleCellExperiment)
  #   library(slingshot)
  #   library(tradeSeq)
  #   library(BiocParallel)
  # })
  #
  # sce <- readH5AD(adata_h5ad, X_name = "X", layers = "counts")
  # sce <- sce[, sce$organ == plan_item$organ &
  #              sce$cell_type_L2 %in% plan_item$anchor_celltypes]
  # reducedDim(sce, "scVI") <- as.matrix(reducedDim(sce, "X_scVI"))
  #
  # # Discrete clusters needed for slingshot start/end
  # sce$cluster <- kmeans(reducedDim(sce, "scVI"),
  #                       centers = length(plan_item$anchor_celltypes))$cluster
  #
  # sds <- slingshot(sce, clusterLabels = "cluster", reducedDim = "scVI",
  #                  start.clus = NULL, end.clus = NULL)
  # saveRDS(sds, file.path(outdir, plan_item$organ, plan_item$name,
  #                        "slingshot_curves.rds"))
  #
  # # tradeSeq::fitGAM with condition factor
  # counts <- assay(sce, "counts")
  # pseudotime <- slingPseudotime(sds)
  # weights <- slingCurveWeights(sds)
  # set.seed(0)
  # fit <- fitGAM(counts = counts,
  #               pseudotime = pseudotime,
  #               cellWeights = weights,
  #               conditions = factor(sce$condition,
  #                                   levels = unlist(cohort$conditions)),
  #               nknots = 6,
  #               parallel = TRUE,
  #               BPPARAM = MulticoreParam(workers = 4))
  # saveRDS(fit, file.path(outdir, plan_item$organ, plan_item$name,
  #                        "tradeseq_fit.rds"))
  #
  # for (cc in cohort$organs[[plan_item$organ]]$contrasts) {
  #   ct <- conditionTest(fit, l2fc = log2(1.5),
  #                       global = TRUE, pairwise = TRUE)
  #   write.table(ct, file.path(outdir, plan_item$organ, plan_item$name,
  #                             sprintf("conditionTest__%s_vs_%s.tsv",
  #                                     cc$test, cc$ref)),
  #               sep = "\t", quote = FALSE, row.names = TRUE)
  # }
  message(sprintf(
    "  [pseudotime_cohort] %s/%s — execute path is a SCAFFOLD. The R code is "
    "laid out inside this function as comments — uncomment + load the listed "
    "packages once xenium_R_analysis env is activated and the cohort_integrated.h5ad "
    "is on disk.",
    plan_item$organ, plan_item$name))
}

main <- function() {
  args <- parse_args()

  if (is.na(args$cohort) || !file.exists(args$cohort)) {
    stop("--cohort cohort_meta.json required")
  }
  cohort <- fromJSON(args$cohort, simplifyVector = FALSE)

  plan <- plan_lineages(cohort, args$organ)
  if (length(plan) == 0) {
    stop(sprintf("no lineages match organ='%s'", args$organ))
  }
  if (args$lineage != "all") {
    plan <- Filter(function(p) p$name == args$lineage, plan)
    if (length(plan) == 0) stop(sprintf("lineage %s not in cohort", args$lineage))
  }

  if (!args$execute) {
    dry_run(plan, args$outdir)
    message("[pseudotime_cohort] dry-run complete. Pass --execute to compute.")
    return(invisible(NULL))
  }

  if (is.na(args$input) || !file.exists(args$input)) {
    stop("--input cohort_integrated.h5ad required for --execute")
  }
  for (p in plan) {
    d <- file.path(args$outdir, p$organ, p$name)
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    execute_lineage(args$input, p, cohort, args$outdir)
  }
}

main()
