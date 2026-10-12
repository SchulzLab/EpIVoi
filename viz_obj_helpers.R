
library(data.table)

get_genes <- function(viz_object) {
  unique(viz_object$feature_importance$file_index$gene_id)
}

get_groups <- function(viz_object) {
  viz_object$feature_importance$groups
}

resolve_epivoi_path <- function(viz_object, path) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Expected one non-empty input file path")
  }

  path <- path.expand(path)

  if (grepl("^/", path)) {
    return(path)
  }

  base_dir <- attr(viz_object, "object_dir", exact = TRUE)

  if (is.null(base_dir) ||
      length(base_dir) != 1L ||
      !nzchar(base_dir)) {
    return(path)
  }

  file.path(base_dir, path)
}

read_feature_importance_file <- function(viz_object, gene) {
  idx <- viz_object$feature_importance$file_index
  path <- idx$file_path[idx$gene_id == gene]

  if (length(path) == 0L) {
    stop("No feature importance file for gene: ", gene)
  }

  path <- resolve_epivoi_path(viz_object, path[1L])

  raw <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )

  groups <- raw[[1L]]
  feature_names <- colnames(raw)[-1L]

  mat <- t(as.matrix(raw[, -1L, drop = FALSE]))
  colnames(mat) <- groups

  out <- as.data.frame(mat, check.names = FALSE)
  out <- cbind(region = feature_names, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL

  out
}

read_epigenetic_signal_file <- function(viz_object, gene) {
  idx <- viz_object$epigenetic_signal$file_index
  path <- idx$file_path[idx$gene_id == gene]

  if (length(path) == 0L) {
    stop("No epigenetic signal file for gene: ", gene)
  }

  path <- resolve_epivoi_path(viz_object, path[1L])

  data.table::fread(path, data.table = FALSE)
}
