options(stringsAsFactors = FALSE)
set.seed(666)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

parse_cli_options <- function(args) {
  opts <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 1) {
      opts[[kv[1]]] <- TRUE
    } else {
      opts[[kv[1]]] <- paste(kv[-1], collapse = "=")
    }
  }
  opts
}

get_opt <- function(opts, dashed, underscored = gsub("-", "_", dashed), default = NULL) {
  if (!is.null(opts[[dashed]])) return(opts[[dashed]])
  if (!is.null(opts[[underscored]])) return(opts[[underscored]])
  default
}

as_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

opts <- parse_cli_options(commandArgs(trailingOnly = TRUE))

required <- c("input-dir", "output-dir", "intermediate-dir", "lsgi-root", "sample-name")
missing <- required[vapply(required, function(x) is.null(get_opt(opts, x)), FUN.VALUE = logical(1))]
if (length(missing) > 0) {
  stop("Missing required options: ", paste(missing, collapse = ", "), call. = FALSE)
}

input_dir <- normalizePath(get_opt(opts, "input-dir"), winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(get_opt(opts, "output-dir"), winslash = "/", mustWork = FALSE)
intermediate_dir <- normalizePath(get_opt(opts, "intermediate-dir"), winslash = "/", mustWork = FALSE)
lsgi_root <- normalizePath(get_opt(opts, "lsgi-root"), winslash = "/", mustWork = TRUE)
sample_name <- get_opt(opts, "sample-name")
image_path <- get_opt(opts, "image-path", default = "")
if (nzchar(image_path) && file.exists(image_path)) {
  image_path <- normalizePath(image_path, winslash = "/", mustWork = TRUE)
} else {
  image_path <- ""
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("readr", "dplyr", "tibble", "ggplot2", "png", "grid", "viridis", "ComplexHeatmap", "reshape2")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))

if (requireNamespace("anticlust", quietly = TRUE)) {
  suppressPackageStartupMessages(library(anticlust))
  grid_clustering_backend <- paste0("anticlust::balanced_clustering ", as.character(utils::packageVersion("anticlust")))
} else {
  balanced_clustering <- function(x, K) {
    x <- as.data.frame(x)
    K <- max(1L, min(as.integer(K), nrow(x)))
    if (K == 1L) {
      return(rep(1L, nrow(x)))
    }
    scaled_x <- scale(x)
    stats::kmeans(scaled_x, centers = K, iter.max = 100, nstart = 5)$cluster
  }
  grid_clustering_backend <- "stats::kmeans fallback"
}

lsgi_script <- file.path(lsgi_root, "R", "LSGI.R")
if (!file.exists(lsgi_script)) {
  stop("Cannot find LSGI source script: ", lsgi_script, call. = FALSE)
}
source(lsgi_script)

spatial_df <- utils::read.csv(file.path(input_dir, "spatial_coords.csv"), check.names = FALSE)
boundary_df <- utils::read.csv(file.path(input_dir, "boundary_labels.csv"), check.names = FALSE)
embedding_df <- utils::read.csv(file.path(input_dir, "cell_component_embeddings.csv"), check.names = FALSE)

decon_cols <- setdiff(colnames(embedding_df), "cell_ID")
common_ids <- Reduce(intersect, list(spatial_df$cell_ID, boundary_df$cell_ID, embedding_df$cell_ID))
if (length(common_ids) < 10) {
  stop("Too few matched spots across LSGI input files.", call. = FALSE)
}

spatial_df <- spatial_df[match(common_ids, spatial_df$cell_ID), , drop = FALSE]
boundary_df <- boundary_df[match(common_ids, boundary_df$cell_ID), , drop = FALSE]
embedding_df <- embedding_df[match(common_ids, embedding_df$cell_ID), , drop = FALSE]

spatial_coords <- spatial_df[, c("X", "Y"), drop = FALSE]
rownames(spatial_coords) <- spatial_df$cell_ID

embeddings <- as.matrix(embedding_df[, decon_cols, drop = FALSE])
storage.mode(embeddings) <- "numeric"
rownames(embeddings) <- embedding_df$cell_ID
embeddings <- embeddings[, colSums(abs(embeddings), na.rm = TRUE) > 0, drop = FALSE]
if (ncol(embeddings) < 1) {
  stop("No non-zero cell-component columns found.", call. = FALSE)
}

n_grids_scale <- as.numeric(get_opt(opts, "n-grids-scale", default = 10))
n_cells_per_meta <- as.numeric(get_opt(opts, "n-cells-per-meta", default = min(50, nrow(spatial_coords))))
r_squared_thresh <- as.numeric(get_opt(opts, "r-squared-thresh", default = 0.3))
minimum_fctr <- as.numeric(get_opt(opts, "minimum-fctr", default = 3))
arrow_length_scale <- as.numeric(get_opt(opts, "arrow-length-scale", default = 1.4))
arrow_linewidth <- as.numeric(get_opt(opts, "arrow-linewidth", default = 1.0))
arrow_head_cm <- as.numeric(get_opt(opts, "arrow-head-cm", default = 0.20))
arrow_closed <- as_bool(get_opt(opts, "arrow-closed", default = TRUE))
reuse_lsgi <- as_bool(get_opt(opts, "reuse-lsgi", default = FALSE))
arrow_type <- if (arrow_closed) "closed" else "open"

lsgi_result_path <- file.path(intermediate_dir, "11_lsgi_cell_component_result.rds.gz")
did_reuse_lsgi <- FALSE
if (reuse_lsgi && file.exists(lsgi_result_path)) {
  lsgi_res <- readr::read_rds(lsgi_result_path)
  did_reuse_lsgi <- TRUE
} else {
  lsgi_res <- local.traj.preprocessing(
    spatial_coords = spatial_coords,
    embeddings = embeddings,
    n.grids.scale = n_grids_scale,
    n.cells.per.meta = n_cells_per_meta
  )
  readr::write_rds(lsgi_res, lsgi_result_path, compress = "gz")
}
utils::write.csv(lsgi_res$grid.info, file.path(output_dir, "grid_info.csv"), row.names = FALSE)

lin_res <- get.ind.rsqrs(lsgi_res)
lin_res <- stats::na.omit(lin_res)
arrow_df <- lin_res[lin_res$rsquared > r_squared_thresh, , drop = FALSE]
if (nrow(arrow_df) > 0) {
  arrow_df <- arrow_df |>
    dplyr::group_by(fctr) |>
    dplyr::filter(dplyr::n() >= minimum_fctr) |>
    dplyr::ungroup() |>
    as.data.frame()
}
utils::write.csv(arrow_df, file.path(output_dir, "cell_component_gradient_arrows.csv"), row.names = FALSE)

dist_mat <- tryCatch(
  avg.dist.calc(lsgi_res, r_squared_thresh = r_squared_thresh, minimum.fctr = minimum_fctr),
  error = function(e) {
    warning("LSGI distance calculation skipped: ", conditionMessage(e), call. = FALSE)
    NULL
  }
)
if (!is.null(dist_mat) && nrow(dist_mat) > 0) {
  utils::write.csv(dist_mat, file.path(output_dir, "cell_component_gradient_distance.csv"), row.names = FALSE)
  grDevices::pdf(file.path(output_dir, "cell_component_gradient_distance_heatmap.pdf"), width = 7, height = 6)
  ComplexHeatmap::draw(plt.dist.heat(dist_mat))
  grDevices::dev.off()
}

boundary_cols <- c(Mal = "#CB181D", Bdy = "#1f78b4", nMal = "#fdb462")
point_df <- spatial_df |>
  dplyr::left_join(boundary_df[, c("cell_ID", "Location"), drop = FALSE], by = "cell_ID") |>
  dplyr::mutate(Location = factor(Location, levels = names(boundary_cols)))

arrow_plot_df <- arrow_df
if (nrow(arrow_plot_df) > 0) {
  arrow_plot_df <- arrow_plot_df |>
    dplyr::mutate(
      X_end = X + vx.u * arrow_length_scale,
      Y_end = Y + vy.u * arrow_length_scale
    )
}

add_gradient_arrows <- function(p, arrow_data = arrow_plot_df) {
  if (nrow(arrow_data) == 0) {
    return(p + ggplot2::labs(subtitle = paste0("No component gradients passed R2 > ", r_squared_thresh)))
  }
  p +
    ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(
        x = X,
        y = Y,
        xend = X_end,
        yend = Y_end,
        color = fctr
      ),
      inherit.aes = FALSE,
      linewidth = arrow_linewidth,
      lineend = "round",
      arrow = ggplot2::arrow(length = grid::unit(arrow_head_cm, "cm"), type = arrow_type)
    ) +
    ggplot2::labs(color = "LSGI component")
}

boundary_base <- ggplot2::ggplot(point_df, ggplot2::aes(x = X, y = Y, fill = Location)) +
  ggplot2::geom_point(shape = 21, size = 1.8, stroke = 0.1, color = "grey25", alpha = 0.9) +
  ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE) +
  ggplot2::scale_y_reverse() +
  ggplot2::coord_fixed() +
  ggplot2::theme_void() +
  ggplot2::theme(legend.position = "right") +
  ggplot2::ggtitle(paste0(sample_name, " boundary with LSGI cell-component gradients"))

boundary_gradient <- add_gradient_arrows(boundary_base)
ggplot2::ggsave(
  file.path(output_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf")),
  boundary_gradient,
  width = 8,
  height = 7
)

if (nzchar(image_path) && file.exists(image_path)) {
  img <- png::readPNG(image_path)
  img_grob <- grid::rasterGrob(
    img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc")
  )
  point_df_he <- point_df |>
    dplyr::mutate(Y = -Y)
  arrow_plot_df_he <- arrow_plot_df
  if (nrow(arrow_plot_df_he) > 0) {
    arrow_plot_df_he <- arrow_plot_df_he |>
      dplyr::mutate(
        Y = -Y,
        Y_end = -Y_end
      )
  }
  he_base <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob = img_grob,
      xmin = 0,
      xmax = ncol(img),
      ymin = -nrow(img),
      ymax = 0
    ) +
    ggplot2::geom_point(
      data = point_df_he,
      ggplot2::aes(x = X, y = Y, fill = Location),
      shape = 21,
      size = 1.8,
      stroke = 0.1,
      color = "grey20",
      alpha = 0.82
    ) +
    ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, ncol(img)),
      ylim = c(-nrow(img), 0),
      expand = FALSE,
      clip = "on"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::ggtitle(paste0(sample_name, " HE-boundary with LSGI cell-component gradients"))

  he_gradient <- add_gradient_arrows(he_base, arrow_plot_df_he)
  ggplot2::ggsave(
    file.path(output_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf")),
    he_gradient,
    width = 8,
    height = 7
  )
}

grDevices::pdf(file.path(output_dir, "cell_component_gradients_plain_lsgi.pdf"), width = 8, height = 7)
print(plt.factors.gradient.ind(
  info = lsgi_res,
  r_squared_thresh = r_squared_thresh,
  minimum.fctr = minimum_fctr,
  arrow.length.scale = arrow_length_scale
) + ggplot2::ggtitle("LSGI cell-component gradients"))
grDevices::dev.off()

sink(file.path(output_dir, "run_summary.txt"))
cat("Cottrazm Python + R LSGI cell-component gradient analysis\n")
cat("=========================================================\n\n")
cat("Sample:", sample_name, "\n")
cat("Matched spots:", nrow(spatial_coords), "\n")
cat("Cell components:", paste(colnames(embeddings), collapse = ", "), "\n")
cat("Grid clustering backend:", grid_clustering_backend, "\n")
cat("n.grids.scale:", n_grids_scale, "\n")
cat("n.cells.per.meta:", n_cells_per_meta, "\n")
cat("R-squared threshold:", r_squared_thresh, "\n")
cat("Minimum arrows per component:", minimum_fctr, "\n")
cat("Selected gradient arrows:", nrow(arrow_df), "\n\n")
cat("Arrow length scale:", arrow_length_scale, "\n")
cat("Arrow linewidth:", arrow_linewidth, "\n")
cat("Arrow head cm:", arrow_head_cm, "\n")
cat("Arrow closed:", arrow_closed, "\n")
cat("Reused LSGI result:", did_reuse_lsgi, "\n\n")
cat("Inputs:\n")
cat("- ", file.path(input_dir, "spatial_coords.csv"), "\n", sep = "")
cat("- ", file.path(input_dir, "cell_component_embeddings.csv"), "\n", sep = "")
cat("- ", file.path(input_dir, "boundary_labels.csv"), "\n\n", sep = "")
cat("Outputs:\n")
cat("- ", lsgi_result_path, "\n", sep = "")
cat("- ", file.path(output_dir, "grid_info.csv"), "\n", sep = "")
cat("- ", file.path(output_dir, "cell_component_gradient_arrows.csv"), "\n", sep = "")
cat("- ", file.path(output_dir, paste0(sample_name, "_BoundaryDefine_LSGIGradient.pdf")), "\n", sep = "")
cat("- ", file.path(output_dir, paste0(sample_name, "_BoundaryDefine_HE_LSGIGradient.pdf")), "\n", sep = "")
sink()

message("Done. LSGI gradient outputs written to ", normalizePath(output_dir, winslash = "/", mustWork = FALSE), ".")
