# Public-release version of the image registration workflow
#
# Notes:
# - Reorganized into functions
# - Replaced hard-coded local paths with configurable inputs
# - Removed private Box / local machine paths
# - Added basic argument validation and clearer object names
# - Kept interactive anchor-selection steps, because they are part of the workflow
#
# Suggested repository structure:
# project/
#   R/
#     image_registration_public_release.R
#   data/
#     contour/
#       contour_img.png
#     anchors/
#     outputs/
#   README.md
#   LICENSE
#
# Example usage:
# source("R/image_registration_public_release.R")
# result <- run_image_registration(config)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(imager)
  library(Morpho)
  library(sp)
  library(purrr)
  library(stringr)
  library(fields)
  library(tibble)
})

# =============================================================================
# Configuration
# =============================================================================

config <- list(
  seurat_rds = "data/input/seurat_object.rds",
  contour_image = "inst/extdata/contour_img_example.png",
  contour_anchor_csv = "inst/extdata/contour_anchors_example.csv",
  output_dir = "data/outputs",
  anchor_dir = "data/anchors",
  sample_id = "sample_01",
  slide_id = "fov",
  genes_of_interest = c("Slc17a6", "Slc32a1", "Chat"),
  n_anchor_points = 24,
  lambda_val = 1e-3,
  seed = 42,
  n_plot_centroids = 20000,
  n_plot_molecules = 100000,
  target_cell_types = 0:5,
  cell_id_field_index = 8,
  cell_id_prefix = NULL,
  cell_type_column = "RNA_snn_res.14"
)

# =============================================================================
# Utilities
# =============================================================================

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

date_stamp <- function() format(Sys.Date(), "%Y%m%d")

stop_if_missing <- function(x, label) {
  if (is.null(x) || length(x) == 0 || any(is.na(x))) {
    stop(sprintf("Missing required value: %s", label), call. = FALSE)
  }
}

rename_cells_by_field <- function(seu, field_index = 8, prefix = NULL) {
  old_ids <- Cells(seu)
  split_ids <- strsplit(old_ids, "_", fixed = TRUE)

  new_ids <- vapply(
    split_ids,
    FUN = function(x) {
      if (length(x) < field_index) {
        stop("Some cell IDs do not contain enough underscore-separated fields.", call. = FALSE)
      }
      x[field_index]
    },
    FUN.VALUE = character(1)
  )

  new_ids <- make.unique(new_ids)
  if (!is.null(prefix)) new_ids <- paste0(prefix, new_ids)
  names(new_ids) <- old_ids

  seu <- RenameCells(object = seu, new.names = new_ids)
  stopifnot(all(Cells(seu) == unname(new_ids)))
  stopifnot(all(rownames(seu@meta.data) %in% Cells(seu)))
  seu
}

rotate_pts <- function(x, y, angle_deg, center = NULL) {
  rad <- angle_deg * pi / 180
  R <- matrix(c(cos(rad), sin(rad), -sin(rad), cos(rad)), ncol = 2)

  if (is.null(center)) {
    center <- c(mean(range(x)), mean(range(y)))
  }

  coords <- cbind(x - center[1], y - center[2]) %*% R
  data.frame(
    x = coords[, 1] + center[1],
    y = coords[, 2] + center[2]
  )
}

warp_points <- function(df, mod_x, mod_y) {
  pts <- as.matrix(df[, c("x", "y")])
  new_x <- predict(mod_x, pts)
  new_y <- predict(mod_y, pts)
  out <- data.frame(x = new_x, y = new_y)
  cbind(out, df[, setdiff(names(df), c("x", "y")), drop = FALSE])
}

extract_selected_molecules <- function(seu, slide_id, genes_of_interest) {
  mol_list <- seu@images[[slide_id]]$molecules

  map_dfr(genes_of_interest, function(g) {
    if (!g %in% names(mol_list)) {
      warning(sprintf("Gene '%s' not found in molecule list; skipping.", g))
      return(NULL)
    }
    spdf <- mol_list[[g]]
    df <- as.data.frame(spdf)
    df$gene <- g
    df
  })
}

extract_centroids <- function(seu, slide_id) {
  centroids_df <- as.data.frame(seu@images[[slide_id]]$centroids@coords)
  centroids_df$cell_id <- rownames(centroids_df)
  centroids_df
}

select_polygon_interactively <- function(df, title_text, n_plot = 20000, color = "red") {
  set.seed(42)
  idx <- sample(nrow(df), size = min(n_plot, nrow(df)))

  plot(
    df$x[idx], df$y[idx],
    pch = 16, cex = 0.5, col = color,
    asp = 1,
    xlab = "x", ylab = "y",
    main = title_text
  )

  poly <- locator(type = "l")
  px <- c(poly$x, poly$x[1])
  py <- c(poly$y, poly$y[1])
  lines(px, py, col = "black", lwd = 2)

  list(px = px, py = py)
}

collect_anchor_points <- function(df, n_anchor_points, title_text, n_plot = 100000) {
  set.seed(42)
  idx <- sample(nrow(df), size = min(n_plot, nrow(df)))

  plot(
    df$x[idx], df$y[idx],
    pch = ".", cex = 1, col = "red",
    asp = 1,
    xlab = "x", ylab = "y",
    main = title_text
  )

  message(sprintf("Click ~%d homologous points on spatial plot...", n_anchor_points))
  anchor_spatial <- as.data.frame(locator(n = n_anchor_points))
  colnames(anchor_spatial) <- c("x", "y")

  message("You clicked ", nrow(anchor_spatial), " points.")
  print(anchor_spatial)
  anchor_spatial
}

fit_tps_models <- function(anchor_spatial, anchor_contour, lambda_val = 1e-3) {
  src <- as.matrix(anchor_spatial[, c("x", "y")])
  tgt <- as.matrix(anchor_contour[, c("x", "y")])

  keep <- !duplicated(data.frame(
    src_x = src[, 1], src_y = src[, 2],
    tgt_x = tgt[, 1], tgt_y = tgt[, 2]
  ))

  src <- src[keep, , drop = FALSE]
  tgt <- tgt[keep, , drop = FALSE]

  tps_mod_x <- Tps(x = src, Y = tgt[, 1], lambda = lambda_val)
  tps_mod_y <- Tps(x = src, Y = tgt[, 2], lambda = lambda_val)

  list(src = src, tgt = tgt, tps_mod_x = tps_mod_x, tps_mod_y = tps_mod_y)
}

compute_anchor_residuals <- function(anchors_tgt_df, anchors_warped) {
  if (!all(anchors_tgt_df$anc_id == anchors_warped$anc_id)) {
    anchors <- merge(
      anchors_tgt_df, anchors_warped,
      by = "anc_id",
      suffixes = c("_tgt", "_warped")
    )
    anchors <- anchors[order(anchors$anc_id), ]
  } else {
    anchors <- data.frame(
      anc_id = anchors_tgt_df$anc_id,
      x_tgt = anchors_tgt_df$x,
      y_tgt = anchors_tgt_df$y,
      x_warp = anchors_warped$x,
      y_warp = anchors_warped$y
    )
  }

  anchors$dx <- anchors$x_warp - anchors$x_tgt
  anchors$dy <- anchors$y_warp - anchors$y_tgt
  anchors$dist <- sqrt(anchors$dx^2 + anchors$dy^2)
  anchors
}

summarize_anchor_residuals <- function(anchors) {
  c(
    n_anchors = nrow(anchors),
    mean_dist = mean(anchors$dist),
    rms_dist = sqrt(mean(anchors$dist^2)),
    median_dist = median(anchors$dist),
    max_dist = max(anchors$dist),
    sd_dist = sd(anchors$dist)
  )
}

# =============================================================================
# Plot helpers
# =============================================================================

plot_contour_with_anchors <- function(contour_img, anchors_tgt_df, anchors_warped) {
  img_w <- dim(contour_img)[1]
  img_h <- dim(contour_img)[2]

  ggplot() +
    annotation_raster(
      as.raster(contour_img),
      xmin = 0, xmax = img_w,
      ymin = 0, ymax = img_h
    ) +
    geom_point(
      data = anchors_tgt_df,
      aes(x = x, y = y),
      colour = "blue", size = 2, shape = 16, alpha = 0.9
    ) +
    geom_point(
      data = anchors_warped,
      aes(x = x, y = y),
      colour = "red", size = 2.5, shape = 1, stroke = 1.1
    ) +
    geom_text(
      data = anchors_tgt_df,
      aes(x = x, y = y, label = anc_id),
      colour = "white", size = 3, fontface = "bold"
    ) +
    geom_text(
      data = anchors_warped,
      aes(x = x, y = y, label = anc_id),
      colour = "black", size = 3
    ) +
    coord_fixed(
      xlim = c(0, img_w),
      ylim = c(0, img_h),
      expand = FALSE
    ) +
    theme_void() +
    labs(title = "Contour with target anchors (blue) and warped source anchors (red)")
}

plot_warped_overlay <- function(contour_img, mol_warped, centroids_warped) {
  img_w <- dim(contour_img)[1]
  img_h <- dim(contour_img)[2]

  ggplot() +
    annotation_raster(
      as.raster(contour_img),
      xmin = 0, xmax = img_w,
      ymin = 0, ymax = img_h
    ) +
    geom_point(
      data = mol_warped,
      aes(x = x, y = y, color = gene),
      size = 0.3, alpha = 0.6
    ) +
    geom_point(
      data = centroids_warped,
      aes(x = x, y = y),
      shape = 1, color = "black", size = 0.8
    ) +
    coord_fixed(
      xlim = c(0, img_w),
      ylim = c(0, img_h),
      expand = FALSE
    ) +
    theme_void() +
    labs(
      title = "Warped spatial data overlaid on contour",
      color = "Gene"
    )
}

plot_selected_cell_types <- function(contour_img, cent_plot, target_types) {
  img_w <- dim(contour_img)[1]
  img_h <- dim(contour_img)[2]

  ggplot() +
    annotation_raster(
      as.raster(contour_img),
      xmin = 0, xmax = img_w,
      ymin = 0, ymax = img_h
    ) +
    geom_point(
      data = cent_plot %>% filter(cell_type %in% target_types),
      aes(x = x, y = y, fill = factor(cell_type)),
      shape = 21,
      color = "black",
      size = 2
    ) +
    facet_wrap(~cell_type, ncol = 3, scales = "fixed") +
    coord_fixed(
      xlim = c(0, img_w),
      ylim = c(0, img_h),
      expand = FALSE
    ) +
    theme_void() +
    guides(fill = "none") +
    labs(title = "Spatial distribution of selected cell types")
}

# =============================================================================
# Main pipeline
# =============================================================================

run_image_registration <- function(config) {
  stop_if_missing(config$seurat_rds, "config$seurat_rds")
  stop_if_missing(config$contour_image, "config$contour_image")
  stop_if_missing(config$output_dir, "config$output_dir")
  stop_if_missing(config$anchor_dir, "config$anchor_dir")
  stop_if_missing(config$slide_id, "config$slide_id")
  stop_if_missing(config$n_anchor_points, "config$n_anchor_points")

  ensure_dir(config$output_dir)
  ensure_dir(config$anchor_dir)

  stamp <- date_stamp()

  message("Loading Seurat object...")
  seu <- readRDS(config$seurat_rds)

  message("Renaming cells...")
  seu <- rename_cells_by_field(
    seu,
    field_index = config$cell_id_field_index,
    prefix = config$cell_id_prefix
  )

  message("Loading contour image...")
  contour_img <- load.image(config$contour_image)
  img_w <- dim(contour_img)[1]
  img_h <- dim(contour_img)[2]

  anchor_contour <- read_contour_anchors(config$contour_anchor_csv) %>%
    mutate(y = img_h - y)

  message("Extracting molecules and centroids...")
  mol_df <- extract_selected_molecules(seu, config$slide_id, config$genes_of_interest)
  centroids_df <- extract_centroids(seu, config$slide_id)

  message("Select region polygon interactively...")
  poly <- select_polygon_interactively(
    centroids_df,
    title_text = "Draw polygon to select region; right-click or ESC when done",
    n_plot = config$n_plot_centroids,
    color = "red"
  )

  inside_mol <- sp::point.in.polygon(mol_df$x, mol_df$y, poly$px, poly$py) == 1
  mol_sel <- as.data.frame(mol_df[inside_mol, ])

  cent_coords <- as.data.frame(seu@images[[config$slide_id]]$centroids@coords)
  cent_coords$cell_id <- seu@images[[config$slide_id]]$centroids@cells
  inside_cent <- sp::point.in.polygon(cent_coords$x, cent_coords$y, poly$px, poly$py) == 1
  cent_sel <- cent_coords[inside_cent, ]

  angle_deg <- as.numeric(readline("Rotation angle in degrees (e.g. 90): "))
  if (is.na(angle_deg)) stop("Rotation angle must be numeric.", call. = FALSE)

  rot_center <- c(
    mean(range(c(mol_sel$x, cent_sel$x))),
    mean(range(c(mol_sel$y, cent_sel$y)))
  )

  mol_rot <- rotate_pts(mol_sel$x, mol_sel$y, angle_deg, center = rot_center) %>%
    bind_cols(mol_sel %>% select(-x, -y))
  cent_rot <- rotate_pts(cent_sel$x, cent_sel$y, angle_deg, center = rot_center) %>%
    bind_cols(cent_sel %>% select(-x, -y))

  p_rot <- ggplot() +
    geom_point(data = mol_rot, aes(x, y), color = "pink", alpha = 0.6, size = 0.7) +
    geom_point(data = cent_rot, aes(x, y), shape = 1, color = "blue", size = 0.2) +
    coord_fixed() +
    theme_void() +
    ggtitle(sprintf("%s | rotated by %d°", config$slide_id, angle_deg))

  ggsave(
    filename = file.path(config$output_dir, sprintf("rotated_%s_%s.png", config$sample_id, stamp)),
    plot = p_rot,
    width = 8,
    height = 4,
    dpi = 300
  )

  anchor_spatial_file <- file.path(
    config$anchor_dir,
    sprintf("anchor_spatial_%s_%s.csv", config$sample_id, stamp)
  )

  anchor_spatial <- collect_anchor_points(
    mol_rot,
    n_anchor_points = config$n_anchor_points,
    title_text = sprintf(
      "Click ~%d homologous points on subsampled cloud; right-click or ESC when done",
      config$n_anchor_points
    ),
    n_plot = config$n_plot_molecules
  )

  save_spatial_anchors(anchor_spatial, anchor_spatial_file)

  tps_fit <- fit_tps_models(
    anchor_spatial = anchor_spatial,
    anchor_contour = anchor_contour,
    lambda_val = config$lambda_val
  )

  mol_warped <- warp_points(mol_rot, tps_fit$tps_mod_x, tps_fit$tps_mod_y)
  centroids_warped <- warp_points(cent_rot, tps_fit$tps_mod_x, tps_fit$tps_mod_y)

  mol_warped$y <- img_h - mol_warped$y
  centroids_warped$y <- img_h - centroids_warped$y

  p_warped <- plot_warped_overlay(contour_img, mol_warped, centroids_warped)
  ggsave(
    filename = file.path(config$output_dir, sprintf("warped_overlay_%s_%s.png", config$sample_id, stamp)),
    plot = p_warped,
    width = 8,
    height = 6,
    dpi = 300
  )

  write.csv(
    mol_warped,
    file.path(config$output_dir, sprintf("mol_warped_%s_%s.csv", config$sample_id, stamp)),
    row.names = FALSE
  )
  write.csv(
    centroids_warped,
    file.path(config$output_dir, sprintf("centroids_warped_%s_%s.csv", config$sample_id, stamp)),
    row.names = FALSE
  )

  if (!config$cell_type_column %in% colnames(seu@meta.data)) {
    stop(sprintf("Column '%s' not found in Seurat metadata.", config$cell_type_column), call. = FALSE)
  }

  seu$cell_type <- seu@meta.data[[config$cell_type_column]]
  meta_df <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    select(cell_id, cell_type)

  centroids_warped$cell_id <- cent_sel$cell_id
  cent_plot <- centroids_warped %>% left_join(meta_df, by = "cell_id")

  write.csv(
    cent_plot,
    file.path(config$output_dir, sprintf("cent_plot_%s_%s.csv", config$sample_id, stamp)),
    row.names = FALSE
  )

  anchors_src_df <- as.data.frame(tps_fit$src)
  names(anchors_src_df) <- c("x", "y")
  anchors_warped <- warp_points(anchors_src_df, tps_fit$tps_mod_x, tps_fit$tps_mod_y)
  anchors_warped$y <- img_h - anchors_warped$y

  anchors_tgt_df <- as.data.frame(tps_fit$tgt)
  names(anchors_tgt_df) <- c("x", "y")
  anchors_tgt_df$y <- img_h - anchors_tgt_df$y
  anchors_tgt_df$anc_id <- seq_len(nrow(anchors_tgt_df))
  anchors_warped$anc_id <- seq_len(nrow(anchors_warped))

  p_anchors <- plot_contour_with_anchors(contour_img, anchors_tgt_df, anchors_warped)
  ggsave(
    filename = file.path(config$output_dir, sprintf("anchor_check_%s_%s.png", config$sample_id, stamp)),
    plot = p_anchors,
    width = 8,
    height = 6,
    dpi = 300
  )

  residuals_df <- compute_anchor_residuals(anchors_tgt_df, anchors_warped)
  residual_summary <- summarize_anchor_residuals(residuals_df)

  write.csv(
    residuals_df,
    file.path(config$output_dir, sprintf("anchor_residuals_%s_%s.csv", config$sample_id, stamp)),
    row.names = FALSE
  )

  p_celltypes <- plot_selected_cell_types(contour_img, cent_plot, config$target_cell_types)
  ggsave(
    filename = file.path(config$output_dir, sprintf("cell_types_%s_%s.png", config$sample_id, stamp)),
    plot = p_celltypes,
    width = 10,
    height = 8,
    dpi = 300
  )

  print(residual_summary)

  invisible(list(
    seu = seu,
    mol_warped = mol_warped,
    centroids_warped = centroids_warped,
    cent_plot = cent_plot,
    residuals = residuals_df,
    residual_summary = residual_summary,
    anchor_spatial = anchor_spatial,
    anchor_contour = anchor_contour
  ))
}
