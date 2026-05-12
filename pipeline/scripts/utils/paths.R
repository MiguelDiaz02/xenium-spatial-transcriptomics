# Centralized path resolution for the Xenium pipeline (R side).
#
# Resolution order for the project root:
#   1. XENIUM_PROJECT_ROOT environment variable (explicit override).
#   2. Walk up from this file's location (Rscript) using --file= arg.
#
# Dataset selection follows the same logic via XENIUM_DATASET
# (default: human_lung_cancer).

.xenium_default_dataset <- "human_lung_cancer"

.xenium_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  }
  # Fall back to sys.frame when sourced interactively.
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  NA_character_
}

xenium_project_root <- function() {
  env <- Sys.getenv("XENIUM_PROJECT_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  # If sourced from pipeline/scripts/utils/paths.R, walk up four levels.
  this <- .xenium_script_path()
  if (!is.na(this) && grepl("pipeline/scripts/utils/paths.R$", this)) {
    return(normalizePath(file.path(dirname(this), "..", "..", ".."),
                         mustWork = FALSE))
  }
  # If sourced from an analysis script, walk up three levels from script.
  if (!is.na(this) && grepl("pipeline/scripts/analysis/", this)) {
    return(normalizePath(file.path(dirname(this), "..", "..", ".."),
                         mustWork = FALSE))
  }
  stop("XENIUM_PROJECT_ROOT must be set when not running via Rscript --file=...")
}

xenium_dataset_root <- function(dataset = NULL) {
  if (is.null(dataset)) {
    dataset <- Sys.getenv("XENIUM_DATASET", unset = .xenium_default_dataset)
  }
  file.path(xenium_project_root(), dataset)
}

xenium_results_root <- function(dataset = NULL) {
  file.path(xenium_dataset_root(dataset), "results")
}
