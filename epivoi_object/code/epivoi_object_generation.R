library(readr)
library(data.table)
library(dplyr)
library(tidyr)
library(parallel)

###### Input ######

# helper function to read an env var, fall back to a default if unset/empty
get_param <- function(env_name, default) {
  val <- Sys.getenv(env_name, unset = "")
  if (nzchar(val)) val else default
}

outer_cores <- as.integer(get_param("EPIVOI_BUILD_OUTER_CORES", "1"))

species <- get_param("EPIVOI_BUILD_SPECIES", "human")
genome_version <- get_param("EPIVOI_BUILD_GENOME_VERSION", "hg38")

#no default
gtf_file <- get_param(
  "EPIVOI_BUILD_GTF_FILE",
  "" 
)

#no default
importance_base_dir <- get_param(
  "EPIVOI_BUILD_IMPORTANCE_BASE_DIR",
  ""
)

importance_file_suffix <- get_param("EPIVOI_BUILD_IMPORTANCE_FILE_SUFFIX", "_shap.txt")

# biological groups from the feature-importance matrices #TODO: optionally set as user parameter
metadata_names <- NULL

# file mapping sample id to corresponding group. 
# Optional summarization epigenetic signal: if left empty, the epigenetic signal is unaggregated 
metadata_file <- get_param("EPIVOI_BUILD_METADATA_FILE", "")

metadata_id_col <- get_param("EPIVOI_BUILD_METADATA_ID_COL", "")     # column in metadata_file holding the sample/metacell ID
metadata_group_col <- get_param("EPIVOI_BUILD_METADATA_GROUP_COL", "")  # column in metadata_file holding the group label

aggr_mode <- get_param("EPIVOI_BUILD_AGGR_MODE", "mean") #for epigenetic signal: default = mean, option = sum

# support sparse matrices (.mtx.gz + rows/cols files), one subfolder per gene
# dense format = default: tab-separated
# "pre_formatted": files are already region-as-rows, one column per
#   sample/group, gzipped TSV -- indexed directly (like feature_importance
#   below), with no reading/reshaping/aggregation performed here. Use this
#   if files are already in the required format upstream (aggregation, if
#   any, is then the user's responsibility, same as for feature importance).
epigenetic_dir <- get_param("EPIVOI_BUILD_EPIGENETIC_DIR", "")
epigenetic_format <- get_param("EPIVOI_BUILD_EPIGENETIC_FORMAT", "dense") # "sparse", "dense", or "pre_formatted"

# only used when epigenetic_format = "pre_formatted"
epigenetic_file_suffix <- get_param("EPIVOI_BUILD_EPIGENETIC_FILE_SUFFIX", "_epigenetic_signal_matrix.tsv.gz")


#only needed when epigenetic_format = "sparse"
epigenetic_mtx_suffix <- get_param("EPIVOI_BUILD_EPIGENETIC_MTX_SUFFIX", ".mtx.gz") #only epigenetic signal
epigenetic_rows_suffix <- get_param("EPIVOI_BUILD_EPIGENETIC_ROWS_SUFFIX", "_rows.txt.gz")
epigenetic_cols_suffix <- get_param("EPIVOI_BUILD_EPIGENETIC_COLS_SUFFIX", "_cols.txt.gz")

epigenetic_flat_response_col <- c("Expression")  # drop response from dense activity files, if present
epigenetic_flat_sample_col <- get_param("EPIVOI_BUILD_EPIGENETIC_SAMPLE_COL", "Sample")  # NULL = use column 1 positionally; set a name to rename that column by name instead
epigenetic_signal_out_dir <- get_param(
  "EPIVOI_BUILD_EPIGENETIC_SIGNAL_OUT_DIR",
  ""
)

#performance
#train_err, test_err
error_metric_name <- get_param("EPIVOI_BUILD_ERROR_METRIC_NAME", "MSE") #default
mse_file <- get_param(
  "EPIVOI_BUILD_MSE_FILE",
  ""
)
#any subset of: train_Pearson, test_Pearson, train_Spearman, test_Spearman
corr_file <- get_param(
  "EPIVOI_BUILD_CORR_FILE",
  ""
)

output_rds <- get_param(
  "EPIVOI_BUILD_OUTPUT_RDS",
  ""
)

log_dir <- get_param("EPIVOI_BUILD_LOG_DIR", "")
 
###### Helper functions ######

# get aggregation mode for epigenetic signal
get_aggr_fun <- function(aggr_mode) {
  switch(
    aggr_mode,
    mean = function(x) mean(x, na.rm = TRUE),
    sum  = function(x) sum(x, na.rm = TRUE),
    stop("Unsupported aggr_mode: ", aggr_mode)
  )
}
 
# read dense matrix holding the epigenetic signal
read_epigenetic_matrix_dense <- function(gene_id, epigenetic_dir, response_col, sample_col = NULL) {
  activity_path <- file.path(epigenetic_dir, paste0(gene_id, ".txt.gz"))
  if (!file.exists(activity_path)) {
    stop("Missing activity file for gene: ", gene_id)
  }
  df <- fread(activity_path, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  df <- df %>% select(-any_of(response_col))
  if (is.null(sample_col)) {
    colnames(df)[1] <- "sample_id"
  } else {
    if (!sample_col %in% colnames(df)) {
      stop("Activity file for gene ", gene_id, " is missing the sample id column: ", sample_col)
    }
    colnames(df)[colnames(df) == sample_col] <- "sample_id"
  }
  return(df)
}
 
# sparse-native group aggregation
    # assumption: response stored separately
aggregate_epigenetic_signal_sparse <- function(gene_id, epigenetic_dir, mtx_suffix, rows_suffix, cols_suffix,
                                                metadata_file, metadata_id_col, metadata_group_col, aggr_mode) {
  gene_dir <- file.path(epigenetic_dir, gene_id)
  mtx_file <- file.path(gene_dir, paste0(gene_id, mtx_suffix))
  rows_file <- file.path(gene_dir, paste0(gene_id, rows_suffix))
  cols_file <- file.path(gene_dir, paste0(gene_id, cols_suffix))
  if (!all(file.exists(c(mtx_file, rows_file, cols_file)))) {
    stop("Missing sparse matrix file(s) for gene: ", gene_id)
  }
 
  mat <- Matrix::readMM(gzfile(mtx_file))  # samples x regions, stays sparse
  rows_df <- fread(rows_file, header = TRUE, sep = "\t", data.table = FALSE)
  cols_df <- fread(cols_file, header = TRUE, sep = "\t", data.table = FALSE)
  if (!"cell_id" %in% colnames(rows_df)) {
    stop("rows file for gene ", gene_id, " is missing the expected 'cell_id' column")
  }
  if (!"bin_id" %in% colnames(cols_df)) {
    stop("cols file for gene ", gene_id, " is missing the expected 'bin_id' column")
  }
 
  groups <- assign_groups(rows_df$cell_id, metadata_file, metadata_id_col, metadata_group_col)
  valid <- !is.na(groups)
  mat <- mat[valid, , drop = FALSE]
  groups <- groups[valid]
  if (length(groups) == 0) {
    stop("No annotated samples for gene: ", gene_id)
  }
 
  unique_groups <- unique(groups)
  group_idx <- match(groups, unique_groups)
 
  # sparse indicator matrix: samples x groups
  indicator <- Matrix::sparseMatrix(
    i = seq_along(group_idx),
    j = group_idx,
    x = 1,
    dims = c(length(group_idx), length(unique_groups)),
    dimnames = list(NULL, unique_groups)
  )
 
  # sum per group in one sparse matrix multiply: regions x groups
  group_sums <- Matrix::t(mat) %*% indicator
 
  if (aggr_mode == "mean") { #default
    group_sizes <- as.numeric(table(groups)[unique_groups])
    group_result <- sweep(as.matrix(group_sums), 2, group_sizes, "/")
  } else if (aggr_mode == "sum") {
    group_result <- as.matrix(group_sums)
  } else {
    stop("Unsupported aggr_mode: ", aggr_mode)
  }
 
  signal_mat <- as.data.frame(group_result)
  signal_mat <- cbind(region = cols_df$bin_id, signal_mat, stringsAsFactors = FALSE)
  return(signal_mat)
}
 
# dense-format group aggregation
aggregate_epigenetic_signal_dense <- function(gene_id, epigenetic_dir, response_col, sample_col,
                                               metadata_file, metadata_id_col, metadata_group_col, aggr_mode) {
  activity_df <- read_epigenetic_matrix_dense(gene_id, epigenetic_dir, response_col, sample_col)
  activity_df$group <- assign_groups(activity_df$sample_id, metadata_file, metadata_id_col, metadata_group_col)
  activity_df <- activity_df[!is.na(activity_df$group), ]
 
  if (nrow(activity_df) == 0) {
    stop("No annotated samples for gene: ", gene_id)
  }
 
  regions <- setdiff(colnames(activity_df), c("sample_id", "group"))
  aggr_fun <- get_aggr_fun(aggr_mode)
 
  signal_mat <- activity_df %>%
    group_by(group) %>%
    summarise(across(all_of(regions), aggr_fun), .groups = "drop") %>%
    pivot_longer(cols = -group, names_to = "region", values_to = "signal") %>%
    pivot_wider(names_from = group, values_from = signal)
 
  return(signal_mat)
}
 
# dense-format, no aggregation
# converted to a dense region x sample data frame for output.
read_epigenetic_signal_dense_raw <- function(gene_id, epigenetic_dir, response_col, sample_col) {
  activity_df <- read_epigenetic_matrix_dense(gene_id, epigenetic_dir, response_col, sample_col)

  regions <- setdiff(colnames(activity_df), "sample_id")
  signal_mat <- activity_df %>%
    pivot_longer(cols = all_of(regions), names_to = "region", values_to = "signal") %>%
    pivot_wider(names_from = sample_id, values_from = signal)

  return(signal_mat)
}

# sparse-format, no aggregation: 
# converted to a dense region x sample data frame for output.
read_epigenetic_signal_sparse_raw <- function(gene_id, epigenetic_dir, mtx_suffix, rows_suffix, cols_suffix) {
  gene_dir <- file.path(epigenetic_dir, gene_id)
  mtx_file <- file.path(gene_dir, paste0(gene_id, mtx_suffix))
  rows_file <- file.path(gene_dir, paste0(gene_id, rows_suffix))
  cols_file <- file.path(gene_dir, paste0(gene_id, cols_suffix))
  if (!all(file.exists(c(mtx_file, rows_file, cols_file)))) {
    stop("Missing sparse matrix file(s) for gene: ", gene_id)
  }

  mat <- Matrix::readMM(gzfile(mtx_file))  # samples x regions
  rows_df <- fread(rows_file, header = TRUE, sep = "\t", data.table = FALSE)
  cols_df <- fread(cols_file, header = TRUE, sep = "\t", data.table = FALSE)
  if (!"cell_id" %in% colnames(rows_df)) {
    stop("rows file for gene ", gene_id, " is missing the expected 'cell_id' column")
  }
  if (!"bin_id" %in% colnames(cols_df)) {
    stop("cols file for gene ", gene_id, " is missing the expected 'bin_id' column")
  }

  signal_mat <- as.data.frame(as.matrix(Matrix::t(mat)))  # regions x samples
  colnames(signal_mat) <- rows_df$cell_id
  signal_mat <- cbind(region = cols_df$bin_id, signal_mat, stringsAsFactors = FALSE)

  return(signal_mat)
}

assign_groups <- function(sample_ids, metadata_file, metadata_id_col, metadata_group_col) {
  meta_df <- fread(metadata_file, data.table = FALSE)
 
  if (!metadata_id_col %in% colnames(meta_df)) {
    stop("metadata_file is missing the id column: ", metadata_id_col)
  }
  if (!metadata_group_col %in% colnames(meta_df)) {
    stop("metadata_file is missing the group column: ", metadata_group_col)
  }
  # build the id -> group dictionary. "names" here are the sample/metacell IDs, "values" are their group labels
  id_to_group <- setNames(meta_df[[metadata_group_col]], meta_df[[metadata_id_col]])
  # look up the group for every sample_id, in the same order as sample_ids. Any ID not found in id_to_group becomes NA here.
  matched_groups <- id_to_group[sample_ids]
  return(matched_groups)
}

#logging 
make_progress_logger <- function(progress_log) {
  function(...) {
    cat(paste(Sys.time(), ...), file = progress_log, append = TRUE, sep = "\n")
  }
}
 
 ###### Build object ######
 
dir.create(epigenetic_signal_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
 
# Aggregation only happens if a metadata file was provided; otherwise the
# raw, per-sample/metacell epigenetic signal is stored as-is.
do_aggregate <- nzchar(metadata_file)

if (do_aggregate) {
  aggr_fun <- get_aggr_fun(aggr_mode)  # validates aggr_mode early, only when it's actually used
  print("metadata_file provided -- epigenetic signal will be aggregated by group.")
} else {
  print("No metadata_file provided -- epigenetic signal will be stored per sample/metacell, unaggregated.")
}
 
#logging
log_file <- file.path(log_dir, paste0("epivoi_obj_build_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
progress_log_path <- file.path(log_dir, "epivoi_obj_build_progress.log")
sink(log_file, append = TRUE, split = TRUE)
on.exit(sink(), add = TRUE)
log_progress <- make_progress_logger(progress_log_path)
 
#gtf file: gene ids/coordinates
gtf_df <- fread(gtf_file, header = FALSE, sep = "\t", data.table = FALSE)
# keep only genes
gtf_df <- gtf_df[gtf_df$V3 == "gene", ]
gtf_df$gene_id <- sub('.*gene_id "([^"]+)".*', '\\1', gtf_df$V9)
gtf_df$gene_name <- sub('.*gene_name "([^"]+)".*', '\\1', gtf_df$V9)
gtf_df$gene_id_no_suffix <- sub("\\..*$", "", gtf_df$gene_id)
print(paste("GTF loaded:", nrow(gtf_df), "genes"))
 
# feature_importance index
#gene_dirs <- list.dirs(importance_base_dir, full.names = TRUE, recursive = FALSE)
gene_files <- list.files(importance_base_dir, pattern=importance_file_suffix, full.names = FALSE)
gene_ids <- sub(importance_file_suffix, "", gene_files)
file_paths <- file.path(importance_base_dir, paste0(gene_ids, importance_file_suffix))
exists_mask <- file.exists(file_paths)
 
feature_importance_index <- data.frame(
  gene_id = gene_ids[exists_mask],
  file_path = file_paths[exists_mask],
  stringsAsFactors = FALSE
)
print(paste("feature_importance index:", nrow(feature_importance_index), "genes found"))
log_progress("feature_importance:", nrow(feature_importance_index), "genes indexed")
 
# metadata groups (only relevant when aggregating)
if (do_aggregate) {
  if (is.null(metadata_names)) {
    meta_df_for_groups <- fread(metadata_file, data.table = FALSE)
    if (!metadata_group_col %in% colnames(meta_df_for_groups)) {
      stop("metadata_file is missing the group column: ", metadata_group_col)
    }
    groups <- unique(meta_df_for_groups[[metadata_group_col]])
    print(paste("Detected", length(groups), "groups from metadata_file"))
  } else {
    groups <- metadata_names
    print(paste("Using", length(groups), "user-supplied groups"))
  }
  log_progress("groups:", paste(groups, collapse = "; "))
} else {
  groups <- NULL
  log_progress("groups: none (no metadata_file provided, signal stored per sample/metacell)")
}

 
# epigenetic_signal
if (epigenetic_format == "pre_formatted") {

  # Files already in the required format index directly
  epigenetic_files <- list.files(epigenetic_dir, pattern = epigenetic_file_suffix, full.names = FALSE)
  epigenetic_gene_ids <- sub(epigenetic_file_suffix, "", epigenetic_files)
  epigenetic_file_paths <- file.path(epigenetic_dir, paste0(epigenetic_gene_ids, epigenetic_file_suffix))
  epigenetic_exists_mask <- file.exists(epigenetic_file_paths)

  epigenetic_signal_index <- data.frame(
    gene_id = epigenetic_gene_ids[epigenetic_exists_mask],
    file_path = epigenetic_file_paths[epigenetic_exists_mask],
    stringsAsFactors = FALSE
  )
  print(paste("epigenetic_signal index (pre_formatted):", nrow(epigenetic_signal_index), "genes found"))
  log_progress("epigenetic_signal: indexed", nrow(epigenetic_signal_index), "pre-formatted files directly")

  # Epigenentic signal aggregation was performed by user upstream
  epigenetic_aggregation_label <- "user-provided"

} else {

  process_gene_epigenetic_signal <- function(gene_id) {
    data.table::setDTthreads(1)

    print(gene_id)

    out_path <- file.path(epigenetic_signal_out_dir, paste0(gene_id, "_epigenetic_signal_matrix.tsv.gz"))

    if (file.exists(out_path)) {
      log_progress(gene_id, "skipped (already exists)")
      return(data.frame(gene_id = gene_id, file_path = out_path, stringsAsFactors = FALSE))
    }

    tryCatch({
      if (do_aggregate) {
        if (epigenetic_format == "sparse") {
          signal_mat <- aggregate_epigenetic_signal_sparse(
            gene_id, epigenetic_dir, epigenetic_mtx_suffix, epigenetic_rows_suffix, epigenetic_cols_suffix,
            metadata_file, metadata_id_col, metadata_group_col, aggr_mode
          )
        } else if (epigenetic_format == "dense") {
          signal_mat <- aggregate_epigenetic_signal_dense(
            gene_id, epigenetic_dir, epigenetic_flat_response_col, epigenetic_flat_sample_col,
            metadata_file, metadata_id_col, metadata_group_col, aggr_mode
          )
        } else {
          stop("Unsupported epigenetic_format: ", epigenetic_format)
        }
      } else {
        # no metadata_file provided -- store the raw, per-sample/metacell signal
        if (epigenetic_format == "sparse") {
          signal_mat <- read_epigenetic_signal_sparse_raw(
            gene_id, epigenetic_dir, epigenetic_mtx_suffix, epigenetic_rows_suffix, epigenetic_cols_suffix
          )
        } else if (epigenetic_format == "dense") {
          signal_mat <- read_epigenetic_signal_dense_raw(
            gene_id, epigenetic_dir, epigenetic_flat_response_col, epigenetic_flat_sample_col
          )
        } else {
          stop("Unsupported epigenetic_format: ", epigenetic_format)
        }
      }

      fwrite(signal_mat, file = out_path, sep = "\t", compress = "gzip")

      log_progress(gene_id, "success")
      return(data.frame(gene_id = gene_id, file_path = out_path, stringsAsFactors = FALSE))

    }, error = function(e) {
      warning(paste("Error building epigenetic_signal for", gene_id, ":", conditionMessage(e)))
      log_progress(gene_id, "error:", conditionMessage(e))
      return(NULL)
    })
  }

  res <- mclapply(feature_importance_index$gene_id, process_gene_epigenetic_signal, mc.cores = outer_cores)

  epigenetic_signal_index <- bind_rows(res)
  print(paste("epigenetic_signal index:", nrow(epigenetic_signal_index), "genes built"))

  epigenetic_aggregation_label <- if (do_aggregate) aggr_mode else "none"
}

#model performance
read_performance_file <- function(file){
    perf_df <- read.csv(file, sep = "\t", header = TRUE,check.names = FALSE, stringsAsFactors = FALSE)
    return(perf_df)
}

mse_df <- read_performance_file(mse_file)
corr_df <- read_performance_file(corr_file)

#TODO: allow more than one error metric?
# required columns for error file: error handling
allowed_mse_cols <- c("train_err", "test_err")

available_test_mse_cols <- intersect(c("test_err"), colnames(mse_df))


if (length(available_test_mse_cols) == 0) {
  stop(
    "Error file must contain a test error column: ",
    "'test_err'. Available columns: ",
    paste(colnames(corr_df), collapse = ", ")
  )
}

has_err_pair <- all(
  c("train_err", "test_err") %in% colnames(corr_df)
)

if (!has_err_pair) {
  warning(
    "No matching train/test error pair was found. ",
    "The error-based Train vs Test option will not be available."
  )
}

allowed_corr_cols <- c(
  "train_Pearson",
  "test_Pearson",
  "train_Spearman",
  "test_Spearman"
)

# required columns for correlation file: error handling
available_test_cols <- intersect(
  c("test_Pearson", "test_Spearman"),
  colnames(corr_df)
)

if (length(available_test_cols) == 0) {
  stop(
    "Correlation file must contain at least one test correlation column: ",
    "'test_Pearson' or 'test_Spearman'. Available columns: ",
    paste(colnames(corr_df), collapse = ", ")
  )
}

has_pearson_pair <- all(
  c("train_Pearson", "test_Pearson") %in% colnames(corr_df)
)

has_spearman_pair <- all(
  c("train_Spearman", "test_Spearman") %in% colnames(corr_df)
)

if (!has_pearson_pair && !has_spearman_pair) {
  warning(
    "No matching train/test correlation pair was found. ",
    "The correlation-based Train vs Test option will not be available."
  )
}

## assemble and save ##
epivoi_obj <- list(
 
  species = species,
  genome = genome_version,
  gtf = gtf_df,

  performance = list(
    error_metric_name = error_metric_name,
    mse         = mse_df,
    correlation = corr_df),
 
  feature_importance = list(
    file_index = feature_importance_index,
    groups = groups
  ),
 
  epigenetic_signal = list(
    file_index = epigenetic_signal_index,
    aggregation = epigenetic_aggregation_label
  )
)
 
class(epivoi_obj) <- "epivoi_object"
 
saveRDS(epivoi_obj, file = output_rds)
print(paste("Wrote epivoi object to", output_rds))