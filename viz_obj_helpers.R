library(data.table)

# All gene IDs available in the object
get_genes <- function(viz_object) {
  unique(viz_object$feature_importance$file_index$gene_id)
}

# All groups (e.g. "Excitatory:Alzheimer's") present in the object
get_groups <- function(viz_object) {
  viz_object$feature_importance$groups
}

read_feature_importance_file <- function(viz_object, gene) {
  idx <- viz_object$feature_importance$file_index
  path <- idx$file_path[idx$gene_id == gene]
  if (length(path) == 0) {
    stop("No feature_importance file found for gene: ", gene)
  }
  # the z-score files are written groups-as-rows, features-as-columns (index
  # = group label, e.g. "Excitatory:Alzheimer's") - the app expects the
  # opposite (a "region" column + one column per group), same shape as the
  # epigenetic_signal files. Transpose here once so every downstream app
  # function can stay unchanged.
  raw <- fread(path[1], data.table = FALSE)
  groups <- raw[[1]]
  feature_names <- colnames(raw)[-1]
  mat <- t(as.matrix(raw[, -1, drop = FALSE]))
  colnames(mat) <- groups
  out <- as.data.frame(mat, check.names = FALSE)
  out <- cbind(region = feature_names, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  return(out)
}

read_epigenetic_signal_file <- function(viz_object, gene) {
  idx <- viz_object$epigenetic_signal$file_index
  path <- idx$file_path[idx$gene_id == gene]
  if (length(path) == 0) {
    stop("No epigenetic_signal file found for gene: ", gene)
  }
  fread(path[1], data.table = FALSE)
}


########

# get_genes <- function(viz_object) {
#         unique(viz_object$shap$file_index$gene_id) #SHAP genes and ATAC genes = same gene set
# }


# read_shap_file <- function(viz_object, gene) {

#     idx <- viz_object$shap$file_index
#     path <- idx$file_path[idx$gene_id == gene]
#     if (length(path) == 0) {
#         stop("No SHAP file found for gene: ", gene)
#     }
#     fread(path[1], data.table = FALSE)
# }

# #TODO: more general function for epigenetic data
# read_atac_file <- function(viz_object, gene) {

#     idx <- viz_object$atac$file_index
#     path <- idx$file_path[idx$gene_id == gene]
#     if (length(path) == 0) {
#         stop("No ATAC file found for gene: ", gene)
#     }
#     fread(path[1], data.table = FALSE)
# }