#!/usr/bin/env Rscript
# F5_regenerate_figures.R
# Re-generates phase5_pseudotime figures from saved CSVs (no Slingshot re-run)
# Fixes: color palette, grid, black background, scatter gradient

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(scater)
  library(SingleCellExperiment)
  library(uwot)
})

.this_file <- normalizePath(sub("^--file=", "",
                                grep("^--file=", commandArgs(trailingOnly = FALSE),
                                     value = TRUE)[1]),
                            mustWork = FALSE)
source(file.path(dirname(.this_file), "..", "utils", "paths.R"))
BASE <- xenium_dataset_root()
PT   <- file.path(BASE, "results/03_phase3_spatial/F5_pseudotime")
FIG  <- file.path(BASE, "results/figures/phase5_pseudotime")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

log_msg <- function(msg) cat(format(Sys.time(), "[%H:%M:%S]"), msg, "\n", flush=TRUE)

# ── Color palette (matches Python ct_color) ──────────────────────────────────
type_colors <- c(
  "Tumor_proliferating"   = "#e41a1c",
  "Tumor_resting"         = "#fb9a99",
  "Macrophage_M1"         = "#ff7f00",
  "Macrophage_M2"         = "#fdbf6f",
  "Monocyte_classical"    = "#1f78b4",
  "CD8_T_cytotoxic"       = "#f4a460",
  "CD8_T_exhausted"       = "#8b4513",
  "CD4_T_helper"          = "#d2691e",
  "Treg"                  = "#6a3d9a",
  "NK_cytotoxic"          = "#cab2d6",
  "cDC1"                  = "#33a02c",
  "cDC2"                  = "#b2df8a",
  "Endothelial_blood"     = "#2ca25f",
  "Endothelial_lymphatic" = "#99d8c9",
  "Epithelial_general"    = "#4169e1",
  "AT1"                   = "#a6cee3",
  "Ciliated"              = "#969696"
)

# Major group assignments for gradient plot
group_map <- c(
  "Tumor_proliferating" = "Tumor", "Tumor_resting" = "Tumor",
  "Macrophage_M1" = "Myeloid", "Macrophage_M2" = "Myeloid",
  "Monocyte_classical" = "Myeloid", "cDC1" = "Myeloid", "cDC2" = "Myeloid",
  "CD8_T_cytotoxic" = "Lymphoid", "CD8_T_exhausted" = "Lymphoid",
  "CD4_T_helper" = "Lymphoid", "Treg" = "Lymphoid", "NK_cytotoxic" = "Lymphoid",
  "Endothelial_blood" = "Stromal", "Endothelial_lymphatic" = "Stromal",
  "Epithelial_general" = "Epithelial", "AT1" = "Epithelial", "Ciliated" = "Epithelial"
)
group_colors <- c(
  "Tumor" = "#e41a1c", "Myeloid" = "#ff7f00",
  "Lymphoid" = "#6a3d9a", "Stromal" = "#2ca25f", "Epithelial" = "#4169e1"
)

# Clean UMAP theme (no grid, white background)
theme_umap <- function() {
  theme_classic(base_size = 13) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.grid        = element_blank(),
      axis.line         = element_line(colour = "#555555", linewidth = 0.4),
      axis.ticks        = element_line(colour = "#555555", linewidth = 0.3),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white"),
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(size = 10, colour = "#555555")
    )
}

# ── Load data ─────────────────────────────────────────────────────────────────
log_msg("Loading pseudotime CSV...")
pt_df <- read.csv(file.path(PT, "pseudotime_with_spatial.csv"), row.names = 1)

log_msg("Loading cell embeddings for UMAP recompute...")
cells <- read.csv(file.path(PT, "cells_for_slingshot.csv"), row.names = 1)
embed_cols <- grep("^novae_dim_", names(cells), value = TRUE)
# keep only rows that are in pt_df (same stratified subset)
cells_sub <- cells[rownames(cells) %in% rownames(pt_df), ]
pt_df <- pt_df[rownames(cells_sub), ]   # align order

log_msg("Computing UMAP (reproducible seed)...")
set.seed(42)
umap_coords <- uwot::umap(as.matrix(cells_sub[, embed_cols]),
                           n_neighbors = 30, min_dist = 0.3,
                           n_threads = 4)
umap_df <- data.frame(
  UMAP1        = umap_coords[, 1],
  UMAP2        = umap_coords[, 2],
  cell_type_L2 = cells_sub$cell_type_L2,
  pseudotime   = pt_df$slingPseudotime_1,
  row.names    = rownames(cells_sub)
)

# ── Fig 1: UMAP colored by cell type (fixed palette + no grid) ───────────────
log_msg("Fig 1: UMAP cell types...")

p1 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cell_type_L2)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_color_manual(values = type_colors, na.value = "grey70") +
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1), ncol = 2,
    title.position = "top")) +
  labs(
    title    = "UMAP — Novae spatial embeddings (TME subset)",
    subtitle = paste0("n = ", nrow(umap_df), " cells  |  17 cell types  |  Slingshot trajectory analysis"),
    x = "UMAP 1", y = "UMAP 2",
    color = "Cell type"
  ) +
  theme_umap()

ggsave(file.path(FIG, "F5_umap_celltypes.png"), p1, width = 12, height = 8, dpi = 200, bg = "white")
log_msg("Saved F5_umap_celltypes.png")

# ── Fig 2: UMAP colored by pseudotime Lineage 1 (fixed grid) ─────────────────
log_msg("Fig 2: UMAP pseudotime L1...")

# Separate assigned and unassigned so assigned renders on top
umap_na  <- umap_df[is.na(umap_df$pseudotime), ]
umap_pt  <- umap_df[!is.na(umap_df$pseudotime), ]
umap_pt  <- umap_pt[order(umap_pt$pseudotime), ]  # low first so high renders on top

p2 <- ggplot() +
  geom_point(data = umap_na, aes(UMAP1, UMAP2), color = "grey85", size = 0.3, alpha = 0.5) +
  geom_point(data = umap_pt, aes(UMAP1, UMAP2, color = pseudotime), size = 0.4, alpha = 0.75) +
  scale_color_viridis(option = "plasma", name = "Pseudotime\n(Lineage 1)") +
  annotate("text", x = min(umap_pt$UMAP1), y = min(umap_pt$UMAP2) - 0.8,
           label = "Start:\nTumor proliferating", color = "#e41a1c",
           size = 3.2, fontface = "bold", hjust = 0) +
  labs(
    title    = "Pseudotime — Lineage 1 (Tumor → Treg)",
    subtitle = paste0(
      sum(!is.na(umap_df$pseudotime)), " cells assigned  |  ",
      sum(is.na(umap_df$pseudotime)), " cells unassigned (grey)  |  Start = Tumor_proliferating  |  End = Treg"),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_umap() +
  theme(legend.position = "right")

ggsave(file.path(FIG, "F5_umap_pseudotime_L1.png"), p2, width = 10, height = 8, dpi = 200, bg = "white")
log_msg("Saved F5_umap_pseudotime_L1.png")

# ── Fig 3: Spatial pseudotime map (white bg, axes, informative) ──────────────
log_msg("Fig 3: Spatial pseudotime map...")

pt_na  <- pt_df[is.na(pt_df$slingPseudotime_1), ]
pt_pt  <- pt_df[!is.na(pt_df$slingPseudotime_1), ]
pt_pt  <- pt_pt[order(pt_pt$slingPseudotime_1), ]

p3 <- ggplot() +
  geom_point(data = pt_na, aes(x, y), color = "grey88", size = 0.15, alpha = 0.4) +
  geom_point(data = pt_pt, aes(x, y, color = slingPseudotime_1), size = 0.2, alpha = 0.7) +
  scale_color_viridis(option = "plasma", name = "Pseudotime\n(Lineage 1)",
                      guide = guide_colorbar(barwidth = 1, barheight = 8,
                                             title.position = "top")) +
  coord_equal() +
  scale_x_continuous(labels = function(x) paste0(round(x/1000, 0), " mm")) +
  scale_y_continuous(labels = function(y) paste0(round(y/1000, 0), " mm")) +
  labs(
    title    = "Spatial distribution of pseudotime — Lineage 1",
    subtitle = "Tumor → Treg gradient  |  Grey = cells not assigned to Lineage 1  |  Purple–yellow = low–high pseudotime",
    x = "X position",
    y = "Y position"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    legend.background = element_rect(fill = "white")
  )

ggsave(file.path(FIG, "F5_spatial_pseudotime.png"), p3, width = 12, height = 9, dpi = 200, bg = "white")
log_msg("Saved F5_spatial_pseudotime.png")

# ── Fig 4: Pseudotime vs spatial distance (hexbin + group LOESS) ─────────────
log_msg("Fig 4: Pseudotime vs spatial gradient...")

pt_grad <- pt_df[!is.na(pt_df$slingPseudotime_1), ]
pt_grad$group <- group_map[pt_grad$cell_type_L2]

p4 <- ggplot(pt_grad, aes(dist_to_tumor, slingPseudotime_1)) +
  # 2D hexbin density background
  geom_hex(bins = 50, alpha = 0.85) +
  scale_fill_distiller(palette = "Greys", direction = 1,
                       name = "Cell density", trans = "log1p") +
  # Per-group smooth trend lines
  geom_smooth(aes(color = group), method = "loess", span = 0.4,
              se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = group_colors, name = "Cell group") +
  # Annotation: tumor centroid region
  annotate("rect", xmin = 0, xmax = 500, ymin = -Inf, ymax = Inf,
           fill = "#e41a1c", alpha = 0.06) +
  annotate("text", x = 250, y = max(pt_grad$slingPseudotime_1, na.rm=TRUE) * 0.95,
           label = "Tumor\ncore", color = "#e41a1c", size = 3, fontface = "bold",
           hjust = 0.5) +
  labs(
    x        = "Distance to tumor centroid (µm)",
    y        = "Pseudotime (Lineage 1)",
    title    = "Pseudotime gradient: tumor core → periphery",
    subtitle = paste0(
      nrow(pt_grad), " cells assigned to Lineage 1  |  ",
      "Hexbin = cell density  |  Lines = LOESS trend per major compartment")
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    legend.background = element_rect(fill = "white"),
    legend.position   = "right"
  )

ggsave(file.path(FIG, "F5_pseudotime_vs_spatial_gradient.png"), p4,
       width = 12, height = 7, dpi = 200, bg = "white")
log_msg("Saved F5_pseudotime_vs_spatial_gradient.png")

log_msg("All figures regenerated successfully.")
