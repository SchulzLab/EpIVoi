# app.R
#
# EpIVoi Shiny application
# unified version:
# - one main viz_obj from formatted_shap_matrix_gtf.rds
# - full SHAP cell types in app
# - reduced / merged cell types in IGV (with Tcell_types)

library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(shinyjs)
library(DT)
library(tibble)
library(purrr)
library(stringr)
library(tidyr)
library(jsonlite)
library(later)
# ========= APP CONFIGURATION ================================================

get_app_dir <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }

  normalizePath(getwd())
}

app_dir <- get_app_dir()

helper_path <- file.path(app_dir, "viz_obj_helpers.R")

if (!file.exists(helper_path)) {
  stop(
    "Could not find viz_obj_helpers.R in the app directory: ",
    helper_path
  )
}

source(helper_path)

# ========= 1. LOAD DATA =====================================================

# ========= 1. LOAD DATA =====================================================

default_viz_obj_path <- paste0(
  "/projects/single_cell_stitchit/work/visualization/object/",
  "brain_binned_rf/output/",
  "epivoi_obj_brain_metafr_all_genes_shap_mc100_1MB500bp.rds"
)

user_args <- commandArgs(trailingOnly = TRUE)

if (length(user_args) >= 1 && nzchar(user_args[1])) {

  # User supplied the EpIVoi object when starting the app
  viz_obj_path <- user_args[1]

} else {

  # Fallback: environment variable or bundled default object
  viz_obj_path <- Sys.getenv(
    "EPIVOI_OBJECT_PATH",
    unset = default_viz_obj_path
  )
}

if (!file.exists(viz_obj_path)) {
  stop(
    "EpIVoi object does not exist: ",
    viz_obj_path,
    "\nProvide the object path when starting the app, e.g.:",
    "\nRscript epivoi_app.R /path/to/epivoi_object.rds",
    "\nor set EPIVOI_OBJECT_PATH."
  )
}

viz_obj <- tryCatch(
  readRDS(viz_obj_path),
  error = function(e) {
    stop(
      "Could not load EpIVoi object from ",
      viz_obj_path,
      ": ",
      conditionMessage(e)
    )
  }
) 

# Error metric label used for user-facing text. Older objects without this
# field remain supported and use a generic fallback label.
error_metric_name <- viz_obj$performance$error_metric_name
if (is.null(error_metric_name) ||
    length(error_metric_name) == 0 ||
    is.na(error_metric_name[1]) ||
    !nzchar(trimws(as.character(error_metric_name[1])))) {
  error_metric_name <- "Prediction error"
} else {
  error_metric_name <- trimws(as.character(error_metric_name[1]))
}

train_error_label <- paste("Train", error_metric_name)
test_error_label <- paste("Test", error_metric_name)
train_error_log_label <- paste0("log10 Train ", error_metric_name)
test_error_log_label <- paste0("log10 Test ", error_metric_name)
train_test_error_title <- paste("Train vs Test", error_metric_name)


# =========================
# REFSEQ FOR VECTOR PDF EXPORT
# =========================

refseq_file <- Sys.getenv(
  "EPIVOI_REFSEQ_FILE",
  unset = file.path(
    app_dir,
    "local_files",
    "ncbiRefSeq.txt.gz"
  )
)

if (!file.exists(refseq_file)) {
  stop(
    "RefSeq annotation file was not found or is not accessible: ",
    refseq_file,
    "\nSet EPIVOI_REFSEQ_FILE to a readable ncbiRefSeq.txt.gz file."
  )
}

refseq_df <- read.delim(
  gzfile(refseq_file, open = "rt"),
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(refseq_df)[1:13] <- c(
  "bin", "name", "chrom", "strand", "txStart", "txEnd",
  "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  "score", "gene_name"
)

parse_refseq_exons <- function(exon_starts, exon_ends) {
  starts <- strsplit(exon_starts, ",", fixed = TRUE)[[1]]
  ends   <- strsplit(exon_ends,   ",", fixed = TRUE)[[1]]

  starts <- starts[nzchar(starts)]
  ends   <- ends[nzchar(ends)]

  n <- min(length(starts), length(ends))

  if (n == 0) {
    return(data.frame(start = numeric(0), end = numeric(0)))
  }

  data.frame(
    start = as.numeric(starts[seq_len(n)]),
    end   = as.numeric(ends[seq_len(n)]),
    stringsAsFactors = FALSE
  )
}

get_refseq_models_in_region <- function(refseq_df, chrom_label, region_start, region_end,
                                        max_genes = 40) {
  df <- refseq_df[
    refseq_df$chrom == chrom_label &
      refseq_df$txEnd >= region_start &
      refseq_df$txStart <= region_end,
    , drop = FALSE
  ]

  if (nrow(df) == 0) return(NULL)

  # one representative transcript per gene:
  # prefer NM_ transcripts, then longest transcript
  df <- df[order(
    df$gene_name,
    !grepl("^NM_", df$name),
    -(df$txEnd - df$txStart)
  ), , drop = FALSE]

  df <- df[!duplicated(df$gene_name), , drop = FALSE]

  df <- df[order(df$txStart, df$txEnd), , drop = FALSE]

  if (nrow(df) > max_genes) {
    df <- df[seq_len(max_genes), , drop = FALSE]
  }

  models <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]

    exons <- parse_refseq_exons(row$exonStarts, row$exonEnds)

    if (nrow(exons) > 0) {
      exons <- exons[
        exons$end >= region_start &
          exons$start <= region_end,
        , drop = FALSE
      ]
    }

    models[[i]] <- list(
      transcript_id = row$name,
      gene_name = row$gene_name,
      strand = row$strand,
      start = row$txStart,
      end = row$txEnd,
      exons = exons
    )
  }

  models
}

# ---- performance -----------------------------------------------------------
mse_df   <- viz_obj$performance$mse
corr_df  <- viz_obj$performance$correlation

mse_df  <- mse_df  |> rownames_to_column("gene_symbol")
corr_df <- corr_df |> rownames_to_column("gene_symbol")

perf_df <- mse_df |>
  left_join(corr_df, by = "gene_symbol")



has_finite_metric <- function(df, column) {
  column %in% colnames(df) && any(is.finite(suppressWarnings(as.numeric(df[[column]]))))
}
has_train_error <- has_finite_metric(perf_df, "train_err")
has_test_error  <- has_finite_metric(perf_df, "test_err")

available_test_correlations <- c(
  Pearson = has_finite_metric(corr_df, "test_Pearson"),
  Spearman = has_finite_metric(corr_df, "test_Spearman")
)
available_test_correlations <- names(available_test_correlations)[available_test_correlations]

available_paired_correlations <- available_test_correlations[vapply(
  available_test_correlations,
  function(metric) {
    has_finite_metric(corr_df, paste0("train_", metric)) &&
      has_finite_metric(corr_df, paste0("test_", metric))
  },
  logical(1)
)]

if (length(available_test_correlations) == 0) {
  stop("No finite test Pearson or Spearman correlations are available in the object.")
}

default_test_correlation <- if ("Pearson" %in% available_test_correlations) "Pearson" else available_test_correlations[1]
default_test_column <- paste0("test_", default_test_correlation)

correlation_test_choices <- stats::setNames(
  paste0("test_", available_test_correlations),
  paste("Test", available_test_correlations, "correlation")
)

train_test_choices <- stats::setNames("mse", error_metric_name)
if (length(available_paired_correlations) > 0) {
  train_test_choices <- c(
    train_test_choices,
    stats::setNames(
      available_paired_correlations,
      paste(available_paired_correlations, "correlation (train vs test)")
    )
  )
}

# ===== colors and class ordering ============================================
corr_levels <- c("failed", "low", "medium", "high")

corr_colors <- c(
  "failed" = "grey70",
  "low"    = "#56B4E9",
  "medium" = "#009E73",
  "high"   = "#D55E00"
)

corr_colors_a11y <- c(
  "failed" = "#FFFFFF",
  "low"    = "#00E5FF",
  "medium" = "#00FF6A",
  "high"   = "#FF7A00"
)
mix_hex <- function(hex, mix_with = "#BDBDBD", amount = 0.55) {
  hex <- gsub("^#", "", hex)
  mix_with <- gsub("^#", "", mix_with)

  h2rgb <- function(h) {
    as.integer(strtoi(c(
      substr(h, 1, 2),
      substr(h, 3, 4),
      substr(h, 5, 6)
    ), 16L))
  }

  a <- h2rgb(hex)
  b <- h2rgb(mix_with)

  out <- round((1 - amount) * a + amount * b)

  sprintf("#%02X%02X%02X", out[1], out[2], out[3])
}
assign_corr_class <- function(df, corr_metric = "test_Pearson", low_cut = 0.3, high_cut = 0.5) {
  if (!corr_metric %in% colnames(df)) {
    stop("Correlation metric not found in data frame: ", corr_metric)
  }

  df |>
    mutate(
      corr_value_for_class = suppressWarnings(as.numeric(as.character(.data[[corr_metric]]))),
      corr_class = case_when(
        is.na(corr_value_for_class) | corr_value_for_class <= 0 ~ "failed",
        corr_value_for_class < low_cut ~ "low",
        corr_value_for_class < high_cut ~ "medium",
        TRUE ~ "high"
      ),
      corr_class = factor(corr_class, levels = corr_levels)
    )
}

# ===== GTF / gene mapping ===================================================

get_gtf_df <- function(x) {
  if (!is.null(x$gtf)) {
    g <- x$gtf
  } else if (isS4(x) && methods::hasSlot(x, "gtf")) {
    g <- methods::slot(x, "gtf")
  } else {
    stop("Could not find 'gtf' in the loaded RDS object.")
  }

  if (is.data.frame(g)) {
    return(as.data.frame(g, stringsAsFactors = FALSE))
  }

  out <- tryCatch(
    as.data.frame(g, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (!is.null(out)) return(out)

  stop("Could not convert 'gtf' slot to data.frame.")
}

get_genome_version <- function(x) {
  genome <- NULL

  if (!is.null(x$genome)) genome <- x$genome
  if (is.null(genome) && isS4(x) && methods::hasSlot(x, "genome")) {
    genome <- methods::slot(x, "genome")
  }

  genome <- as.character(genome)[1]

  if (is.na(genome) || !nzchar(genome)) {
    stop("Genome version missing in viz_obj.")
  }

  genome
}

get_species <- function(x) {
  species <- NULL

  if (!is.null(x$species)) species <- x$species

  if (is.null(species) && !is.null(x$metadata) && !is.null(x$metadata$species)) {
    species <- x$metadata$species
  }

  if (is.null(species) && isS4(x) && methods::hasSlot(x, "species")) {
    species <- methods::slot(x, "species")
  }

  species <- as.character(species)[1]

  if (!is.na(species) && nzchar(species)) {
    species_lower <- tolower(trimws(species))

    if (species_lower %in% c(
      "human",
      "homo sapiens",
      "homo_sapiens",
      "hsapiens"
    )) {
      return("Homo sapiens")
    }

    if (species_lower %in% c(
      "mouse",
      "mus musculus",
      "mus_musculus",
      "mmusculus"
    )) {
      return("Mus musculus")
    }

    return(species)
  }

  genome <- get_genome_version(x)

  if (genome %in% c("hg19", "hg38")) {
    return("Homo sapiens")
  }

  if (genome %in% c("mm10", "mm39")) {
    return("Mus musculus")
  }

  NA_character_
}

get_species_resources <- function(species_name, genome_version) {
  list(
    species = as.character(species_name)[1],
    genome = as.character(genome_version)[1]
  )
}

extract_gtf_attr <- function(attr_vec, key) {
  attr_vec <- as.character(attr_vec)

  pat1 <- paste0('(?:^|;\\s*)', key, ' "([^"]+)"')
  m1 <- regexec(pat1, attr_vec, perl = TRUE)
  r1 <- regmatches(attr_vec, m1)

  out <- rep(NA_character_, length(attr_vec))
  ok1 <- lengths(r1) >= 2
  out[ok1] <- vapply(r1[ok1], `[`, character(1), 2)

  miss <- is.na(out)
  if (any(miss)) {
    pat2 <- paste0('(?:^|;\\s*)', key, '=([^;]+)')
    m2 <- regexec(pat2, attr_vec[miss], perl = TRUE)
    r2 <- regmatches(attr_vec[miss], m2)
    ok2 <- lengths(r2) >= 2
    tmp <- rep(NA_character_, sum(miss))
    tmp[ok2] <- vapply(r2[ok2], `[`, character(1), 2)
    out[miss] <- tmp
  }

  trimws(out)
}

strip_ensembl_version <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+$", "", x)
}

make_gene_mapping_from_gtf <- function(gtf_df, gene_ids = NULL) {
  gtf <- as.data.frame(gtf_df, stringsAsFactors = FALSE)
  colnames(gtf) <- trimws(colnames(gtf))

  feature_candidates <- c("type", "feature", "V3", "X3")
  feature_col <- intersect(feature_candidates, colnames(gtf))
  feature_col <- if (length(feature_col) > 0) feature_col[1] else NULL

  if (!is.null(feature_col)) {
    feature_vals <- as.character(gtf[[feature_col]])
    keep <- !is.na(feature_vals) & tolower(feature_vals) == "gene"
    if (any(keep)) {
      gtf <- gtf[keep, , drop = FALSE]
    }
  }

  gene_id_candidates <- c(
    "gene_id", "geneid", "ensembl_gene_id",
    "gene", "geneID", "GENEID"
  )
  gene_name_candidates <- c(
    "gene_name", "gene", "gene_symbol",
    "external_gene_name", "symbol", "geneName", "GENENAME"
  )

  gene_id_col <- intersect(gene_id_candidates, colnames(gtf))
  gene_name_col <- intersect(gene_name_candidates, colnames(gtf))

  gene_id_col <- if (length(gene_id_col) > 0) gene_id_col[1] else NULL
  gene_name_col <- if (length(gene_name_col) > 0) gene_name_col[1] else NULL

  gene_id <- NULL
  gene_name <- NULL

  if (!is.null(gene_id_col)) {
    gene_id <- as.character(gtf[[gene_id_col]])
  }

  if (!is.null(gene_name_col)) {
    gene_name <- as.character(gtf[[gene_name_col]])
  }

  if (is.null(gene_id) || all(is.na(gene_id) | !nzchar(gene_id))) {
    attr_candidates <- c("attribute", "attributes", "group", "V9", "X9")
    attr_col <- intersect(attr_candidates, colnames(gtf))

    if (length(attr_col) == 0 && ncol(gtf) >= 9) {
      attr_col <- colnames(gtf)[ncol(gtf)]
    } else {
      attr_col <- attr_col[1]
    }

    if (!is.null(attr_col) && length(attr_col) == 1 && attr_col %in% colnames(gtf)) {
      attr_vec <- as.character(gtf[[attr_col]])
      gene_id <- extract_gtf_attr(attr_vec, "gene_id")
      gene_name <- extract_gtf_attr(attr_vec, "gene_name")
    }
  }

  if (is.null(gene_id) || all(is.na(gene_id) | !nzchar(gene_id))) {
    stop("Could not find gene_id information in gtf.")
  }

  gene_id <- strip_ensembl_version(gene_id)

  if (is.null(gene_name)) {
    gene_name <- rep(NA_character_, length(gene_id))
  }

  map <- data.frame(
    gene_id = as.character(gene_id),
    gene_symbol_display = as.character(gene_name),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(gene_id), nzchar(gene_id)) |>
    dplyr::mutate(
      gene_symbol_display = dplyr::if_else(
        !is.na(gene_symbol_display) & nzchar(gene_symbol_display),
        gene_symbol_display,
        gene_id
      )
    ) |>
    dplyr::distinct(gene_id, .keep_all = TRUE)

  if (!is.null(gene_ids)) {
    gene_ids2 <- strip_ensembl_version(gene_ids)

    map <- data.frame(
      gene_id = unique(gene_ids2),
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(map, by = "gene_id") |>
      dplyr::mutate(
        gene_symbol_display = dplyr::if_else(
          !is.na(gene_symbol_display) & nzchar(gene_symbol_display),
          gene_symbol_display,
          gene_id
        )
      )
  }

  map
}

# ===== SHAP prep ============================================================

parse_shap_id <- function(id_vec) {
  m <- regexec("^(.*)-([^-]+)-([^-]+)-([^-]+)$", id_vec)
  parts <- regmatches(id_vec, m)

  ok <- lengths(parts) == 5

  gene_id <- rep(NA_character_, length(id_vec))
  interaction <- rep(NA_character_, length(id_vec))

  gene_id[ok] <- vapply(parts[ok], `[`, character(1), 2)
  chrom       <- vapply(parts[ok], `[`, character(1), 3)
  start       <- vapply(parts[ok], `[`, character(1), 4)
  end         <- vapply(parts[ok], `[`, character(1), 5)

  interaction[ok] <- paste(chrom, start, end, sep = "-")

  data.frame(
    gene_id = strip_ensembl_version(gene_id),
    interaction = interaction,
    stringsAsFactors = FALSE
  )
}

parse_interaction_safe <- function(x, default_width = 100) {
  x <- gsub("[:_-]", " ", x)
  parts <- strsplit(x, "\\s+")

  chrom <- sapply(parts, `[`, 1)
  start <- suppressWarnings(as.numeric(sapply(parts, `[`, 2)))
  end   <- suppressWarnings(as.numeric(sapply(parts, `[`, 3)))

  end[is.na(end) & !is.na(start)] <-
    start[is.na(end) & !is.na(start)] + default_width

  data.frame(
    chrom = chrom,
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}

gtf_df <- get_gtf_df(viz_obj)
genome_version <- get_genome_version(viz_obj)
species_name <- get_species(viz_obj)
species_resources <- get_species_resources(species_name, genome_version)
ensembl_to_symbol <- gtf_df |>
  dplyr::transmute(
    gene_id = strip_ensembl_version(as.character(gene_id)),
    gene_symbol = as.character(gene_name)
  ) |>
  dplyr::filter(
    !is.na(gene_id),
    nzchar(gene_id),
    !is.na(gene_symbol),
    nzchar(gene_symbol)
  ) |>
  dplyr::distinct(gene_id, .keep_all = TRUE)

ensembl_to_symbol_map <- stats::setNames(
  ensembl_to_symbol$gene_symbol,
  ensembl_to_symbol$gene_id
)

# ===== NEW SHAP / ATAC on-disk object structure =============================

genes_from_obj <- get_genes(viz_obj)
genes_from_obj <- strip_ensembl_version(genes_from_obj)

gene_map <- make_gene_mapping_from_gtf(gtf_df, genes_from_obj)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

symbol2id <- setNames(
  gene_map$gene_id,
  toupper(gene_map$gene_symbol_display)
)

id2symbol <- setNames(
  gene_map$gene_symbol_display,
  gene_map$gene_id
)
perf_df <- perf_df |>
  dplyr::mutate(
    gene_id = dplyr::case_when(
      gene_symbol %in% gene_map$gene_id ~ gene_symbol,
      toupper(gene_symbol) %in% names(symbol2id) ~
        unname(symbol2id[toupper(gene_symbol)]),
      TRUE ~ NA_character_
    ),
    gene_symbol = dplyr::case_when(
      gene_symbol %in% gene_map$gene_id ~
        unname(id2symbol[gene_symbol]),
      TRUE ~ gene_symbol
    )
  ) |>
  dplyr::relocate(gene_id, gene_symbol)

igv_gene_list <- sort(unique(na.omit(gene_map$gene_symbol_display)))


extract_region_matrix <- function(mat, file_type = "SHAP") {
  mat <- as.data.frame(mat, check.names = FALSE)

  if ("region" %in% colnames(mat)) {
    regions <- as.character(mat$region)
    mat$region <- NULL
  } else if ("id" %in% colnames(mat)) {
    regions <- as.character(mat$id)
    mat$id <- NULL
  } else {
    stop("Could not find region column in ", file_type, " matrix.")
  }

  list(
    regions = regions,
    matrix = mat
  )
}

example_gene <- genes_from_obj[1]
example_shap_raw <- read_feature_importance_file(viz_obj, example_gene)
example_shap <- extract_region_matrix(example_shap_raw, file_type = "SHAP")$matrix

app_cell_types <- colnames(example_shap)
app_cell_types <- app_cell_types[!is.na(app_cell_types) & nzchar(app_cell_types)]
##### Laura: conditions #####
app_conditions <- unique(sub(".*:", "", app_cell_types))
############################

# Laura: Why?
tcell_types <- c(
  "CD4_Naive", "CD4_TCM", "CD4_TEM",
  "CD8_Naive", "CD8_TEM_1", "CD8_TEM_2",
  "MAIT", "Treg", "gdT"
)

tcell_types_present <- intersect(tcell_types, app_cell_types)

igv_cell_types <- app_cell_types

if (length(tcell_types_present) > 0) {
  igv_cell_types <- setdiff(igv_cell_types, tcell_types_present)
  igv_cell_types <- c(igv_cell_types, "Tcell_types")
}

load_shap_df <- function(gene, ct) {
  gene <- strip_ensembl_version(gene)

  if (!gene %in% genes_from_obj) {
    stop("Gene not found in viz object: ", gene)
  }

  raw <- read_feature_importance_file(viz_obj, gene)
  x <- extract_region_matrix(raw, file_type = "SHAP")

  mat <- x$matrix
  interactions <- x$regions

  if (!ct %in% colnames(mat)) {
    stop("Cell type not found in Feature Importance matrix: ", ct)
  }

  gene_symbol <- id2symbol[[gene]] %||% gene

  data.frame(
    gene_id = gene,
    gene_symbol = gene_symbol,
    interaction = interactions,
    score = suppressWarnings(as.numeric(mat[[ct]])),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(gene_id),
      !is.na(interaction),
      !is.na(gene_symbol),
      is.finite(score)
    )
}

load_igv_shap_df <- function(gene, ct) {
  gene <- strip_ensembl_version(gene)

  if (!gene %in% genes_from_obj) {
    stop("Gene not found in viz object: ", gene)
  }

  raw <- read_feature_importance_file(viz_obj, gene)
  x <- extract_region_matrix(raw, file_type = "SHAP")

  mat <- x$matrix
  interactions <- x$regions

  if (ct == "Tcell_types") {
    use_cts <- intersect(tcell_types_present, colnames(mat))

    if (length(use_cts) == 0) {
      stop("No T cell columns found for merged Tcell_types.")
    }

    score <- rowMeans(
      mat[, use_cts, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    if (!ct %in% colnames(mat)) {
      stop("Cell type not found in IGV Feature Importance matrix: ", ct)
    }

    score <- mat[[ct]]
  }

  gene_symbol <- id2symbol[[gene]] %||% gene

  data.frame(
    gene_id = gene,
    gene_symbol = gene_symbol,
    interaction = interactions,
    score = suppressWarnings(as.numeric(score)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(gene_id),
      !is.na(interaction),
      !is.na(gene_symbol),
      is.finite(score)
    )
}

load_atac_df <- function(gene, ct) {
  gene <- strip_ensembl_version(gene)

  if (!gene %in% genes_from_obj) {
    stop("Gene not found in viz object: ", gene)
  }

  raw <- read_epigenetic_signal_file(viz_obj, gene)
  x <- extract_region_matrix(raw, file_type = "ATAC")

  mat <- x$matrix
  interactions <- x$regions

  if (ct == "Tcell_types") {
    use_cts <- intersect(tcell_types_present, colnames(mat))

    if (length(use_cts) == 0) {
      stop("No T cell columns found for merged Tcell_types in Epigenetic Signal.")
    }

    score <- rowMeans(
      mat[, use_cts, drop = FALSE],
      na.rm = TRUE
    )
  } else {
    if (!ct %in% colnames(mat)) {
      stop("Cell type not found in Epigenetic Signal matrix: ", ct)
    }

    score <- mat[[ct]]
  }

  gene_symbol <- id2symbol[[gene]] %||% gene

  data.frame(
    gene_id = gene,
    gene_symbol = gene_symbol,
    interaction = interactions,
    score = suppressWarnings(as.numeric(score)),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(gene_id),
      !is.na(interaction),
      !is.na(gene_symbol),
      is.finite(score)
    )
}

igv_gene_list <- sort(unique(na.omit(gene_map$gene_symbol_display)))



make_safe_id <- function(x) {
  x <- as.character(x %||% "NA")
  gsub("[^A-Za-z0-9_]+", "_", x)
}




download_buttons_ui <- function(prefix, label = "Download table") {
  tagList(
    h5(label),
    fluidRow(
      column(
        3,
        downloadButton(
          paste0(prefix, "_all_csv"),
          "Download CSV"
        )
      ),
      column(
        3,
        downloadButton(
          paste0(prefix, "_all_tsv"),
          "Download TSV"
        )
      )
    ),
    br()
  )
}



resolve_igv_gene_input <- function(gene_input) {
  if (is.null(gene_input) || !nzchar(gene_input)) return(NULL)

  x <- trimws(gene_input)
  x_upper <- toupper(x)

  if (x_upper %in% names(symbol2id)) {
    gene_id <- unname(symbol2id[[x_upper]])
    return(list(gene_id = gene_id, gene_display = id2symbol[[gene_id]] %||% x))
  }

  if (x %in% gene_map$gene_id) {
    return(list(gene_id = x, gene_display = id2symbol[[x]] %||% x))
  }

  NULL
}

# ===== file locations =======================================================

bw_dir <- file.path(app_dir, "shap_bw")
if (!dir.exists(bw_dir)) {
  dir.create(bw_dir, recursive = TRUE)
}
shiny::addResourcePath("shap_bw", normalizePath(bw_dir))

igv_www_dir <- Sys.getenv(
  "EPIVOI_IGV_ASSETS_DIR",
  unset = file.path(
    app_dir,
    "igv"
  )
)

if (!dir.exists(igv_www_dir)) {
  stop(
    "IGV assets directory was not found or is not accessible: ",
    igv_www_dir,
    "\nSet EPIVOI_IGV_ASSETS_DIR to the folder containing igv.min.js and igv.min.css."
  )
}

shiny::addResourcePath(
  "igv",
  normalizePath(igv_www_dir)
)
fimo_www_dir <- file.path(app_dir, "fimo")
if (!dir.exists(fimo_www_dir)) {
  dir.create(fimo_www_dir, recursive = TRUE)
}
shiny::addResourcePath("fimo", normalizePath(fimo_www_dir))
ct_to_file <- function(ct) make_safe_id(ct)

regulatory_build_bed <- Sys.getenv(
  "EPIVOI_REGULATORY_BUILD_BED",
  unset = ""
)
# ===== PASTAA / TRAP locations ==============================================

pastaa_root <- file.path(app_dir, "tf_pastaa_app_files")
pastaa_tmp_dir <- file.path(pastaa_root, "tmp")
pastaa_out_dir <- file.path(pastaa_root, "results")

dir.create(pastaa_tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pastaa_out_dir, recursive = TRUE, showWarnings = FALSE)

pastaa_trap_bin <- Sys.getenv(
  "EPIVOI_TRAP_BIN",
  unset = file.path(
    app_dir,
    "PASTAA_local",
    "TRAP"
  )
)

pastaa_bin <- Sys.getenv(
  "EPIVOI_PASTAA_BIN",
  unset = file.path(
    app_dir,
    "PASTAA_local",
    "PASTAA"
  )
)

pastaa_energy_matrix <- Sys.getenv(
  "EPIVOI_PASTAA_ENERGY_MATRIX",
  unset = file.path(
    app_dir,
    "tf_pastaa_output",
    "matrices",
    "Jaspar_Hocomoco_Kellis_human_energy.txt"
  )
)

pastaa_hg38_fa <- Sys.getenv(
  "EPIVOI_GENOME_FASTA",
  unset = file.path(
    app_dir,
    "local_resources",
    "hg38.fa"
  )
)

FIMO_BIN <- Sys.getenv(
  "EPIVOI_FIMO_BIN",
  unset = Sys.which("fimo")
)
if (!nzchar(FIMO_BIN) || !file.exists(FIMO_BIN)) {
  warning(
    "FIMO executable was not found. ",
    "Set EPIVOI_FIMO_BIN or activate an environment containing FIMO."
  )
}

MEME_MOTIF_FILE <- Sys.getenv(
  "EPIVOI_MEME_MOTIF_FILE",
  unset = file.path(
    app_dir,
    "tf_enrichment",
    "JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.meme"
  )
)
fimo_out_dir <- file.path(pastaa_root, "fimo")
dir.create(fimo_out_dir, recursive = TRUE, showWarnings = FALSE)

fimo_www_dir <- file.path(app_dir, "fimo")
dir.create(fimo_www_dir, recursive = TRUE, showWarnings = FALSE)
shiny::addResourcePath("fimo", normalizePath(fimo_www_dir))

FASTA_GET_MARKOV_BIN <- Sys.getenv(
  "EPIVOI_FASTA_GET_MARKOV_BIN",
  unset = Sys.which("fasta-get-markov")
)
if (
  !nzchar(FASTA_GET_MARKOV_BIN) ||
  !file.exists(FASTA_GET_MARKOV_BIN)
) {
  warning(
    "fasta-get-markov executable was not found. ",
    "Set EPIVOI_FASTA_GET_MARKOV_BIN or activate an environment containing MEME Suite."
  )
}

FIMO_ENCODE_BG <- Sys.getenv(
  "EPIVOI_FIMO_ENCODE_BG",
  unset = file.path(
    app_dir,
    "tf_enrichment",
    "Jaspar_Hocomoco_Kellis_human_meme_bgSCREEN3.txt"
  )
)
FIMO_GENOME_BG <- file.path(fimo_out_dir, "hg38_background.txt")
external_files_to_check <- c(
  "TRAP executable" = pastaa_trap_bin,
  "PASTAA executable" = pastaa_bin,
  "PASTAA energy matrix" = pastaa_energy_matrix,
  "Genome FASTA" = pastaa_hg38_fa,
  "MEME motif file" = MEME_MOTIF_FILE,
  "FIMO ENCODE background" = FIMO_ENCODE_BG
  
)

for (resource_name in names(external_files_to_check)) {
  resource_path <- external_files_to_check[[resource_name]]

  if (
    is.na(resource_path) ||
    !nzchar(resource_path) ||
    !file.exists(resource_path)
  ) {
    warning(
      resource_name,
      " was not found or is not accessible: ",
      resource_path
    )
  }
}

ensure_fimo_bg_file <- function(bg_type) {
  bg_type <- match.arg(bg_type, c("encode", "genome"))

  if (bg_type == "encode") {
    if (!file.exists(FIMO_ENCODE_BG)) {
      stop("ENCODE FIMO background file not found: ", FIMO_ENCODE_BG)
    }
    return(FIMO_ENCODE_BG)
  }

  if (!file.exists(FIMO_GENOME_BG)) {
    if (!file.exists(pastaa_hg38_fa)) {
      stop("hg38 FASTA not found: ", pastaa_hg38_fa)
    }

    if (!file.exists(FASTA_GET_MARKOV_BIN)) {
      stop("fasta-get-markov binary not found: ", FASTA_GET_MARKOV_BIN)
    }

    bg_status <- system2(
      FASTA_GET_MARKOV_BIN,
      args = c(pastaa_hg38_fa, FIMO_GENOME_BG)
    )

    if (!identical(bg_status, 0L) || !file.exists(FIMO_GENOME_BG)) {
      stop("fasta-get-markov failed. Could not create: ", FIMO_GENOME_BG)
    }
  }

  FIMO_GENOME_BG
}

fimo_bg_label <- function(bg_type) {
  dplyr::recode(
    bg_type,
    "encode" = "ENCODE background",
    "genome" = "whole genome background",
    .default = bg_type
  )
}
# ========= 2. USER INTERFACE ===============================================

ui <- tagList(

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "igv/igv.min.css"),
    tags$script(src = "igv/igv.min.js"),
    tags$style(HTML("
      .navbar-default {
        background: linear-gradient(90deg, #243B53, #486581) !important;
        border: none !important;
        border-radius: 0 !important;
      }

      .navbar-default .navbar-brand {
        color: #ffffff !important;
        font-weight: 700 !important;
        letter-spacing: 0.2px;
      }

      .navbar-default .navbar-nav > li > a {
        color: #f0f4f8 !important;
        font-weight: 600 !important;
      }

      .navbar-default .navbar-nav > li > a:hover,
      .navbar-default .navbar-nav > li > a:focus {
        background-color: rgba(255, 255, 255, 0.14) !important;
        color: #ffffff !important;
      }

      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:hover,
      .navbar-default .navbar-nav > .active > a:focus {
        background-color: #102A43 !important;
        color: #ffffff !important;
      }

      .navbar {
        margin-bottom: 18px !important;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
      }

      body.a11y, body.a11y .container-fluid, body.a11y .navbar, body.a11y .tab-content {
        background: #000 !important;
        color: #fff !important;
      }

      body.a11y a { color: #00E5FF !important; text-decoration: underline; }

      body.a11y .well, body.a11y .panel, body.a11y .panel-body,
      body.a11y .sidebar, body.a11y .sidebarPanel, body.a11y .mainPanel {
        background: #000 !important;
        color: #fff !important;
        border-color: #fff !important;
      }

      body.a11y input, body.a11y select, body.a11y textarea {
        background: #000 !important;
        color: #fff !important;
        border: 2px solid #fff !important;
      }

      body.a11y .btn, body.a11y .btn-default, body.a11y .btn-primary {
        background: #000 !important;
        color: #fff !important;
        border: 2px solid #fff !important;
        font-weight: 700;
      }

      body.a11y *:focus {
        outline: 4px solid #FFD400 !important;
        outline-offset: 2px;
      }

      body.a11y {
        font-size: 18px;
        line-height: 1.35;
      }

      body.a11y h3 { font-size: 28px; }
      body.a11y h4 { font-size: 22px; }

      body.a11y table.dataTable {
        background: #000 !important;
        color: #fff !important;
        border-collapse: collapse !important;
      }

      body.a11y table.dataTable thead th {
        background: #111 !important;
        color: #fff !important;
        border-bottom: 2px solid #fff !important;
      }

      body.a11y table.dataTable tbody td {
        background: #000 !important;
        color: #fff !important;
        border-top: 1px solid #333 !important;
      }

      body.a11y table.dataTable tbody tr:nth-child(odd) td {
        background: #070707 !important;
      }

      body.a11y table.dataTable tbody tr:nth-child(even) td {
        background: #000 !important;
      }

      body.a11y table.dataTable tbody tr:hover td {
        background: #1a1a1a !important;
      }

      body.a11y .dataTables_wrapper,
      body.a11y .dataTables_info,
      body.a11y .dataTables_length,
      body.a11y .dataTables_filter,
      body.a11y .dataTables_paginate {
        color: #fff !important;
      }

      body.a11y .dataTables_filter input,
      body.a11y .dataTables_length select {
        background: #000 !important;
        color: #fff !important;
        border: 2px solid #fff !important;
      }

      body.a11y .dataTables_filter input {
        width: 260px;
      }

      body.a11y .dataTables_wrapper .dataTables_paginate .paginate_button {
        color: #fff !important;
        background: #000 !important;
        border: 2px solid #fff !important;
        margin: 0 2px;
      }

      body.a11y .dataTables_wrapper .dataTables_paginate .paginate_button.current,
      body.a11y .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover {
        background: #FFD400 !important;
        color: #000 !important;
        border-color: #FFD400 !important;
      }

      body.a11y .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
        background: #111 !important;
        color: #fff !important;
      }

      body.a11y .shiny-plot-output {
        border: 2px solid #fff;
        padding: 6px;
      }

      .epivoi-track-label,
      .epivoi-track-label-box {
        font-weight: 600 !important;
        font-size: 4px !important;
        line-height: 0.9 !important;

        padding: 0px 1px !important;
        border-radius: 1px !important;

        max-width: 90px !important;
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;

        min-height: 0 !important;
        height: auto !important;

        box-shadow: none !important;
        z-index: 9999 !important;

        transform: translateY(-30px) scale(0.70) !important;
        transform-origin: top left !important;
      }

      .epivoi-track-label-box {
        position: relative !important;
        display: inline-block !important;
      }

      body.a11y .epivoi-track-label,
      body.a11y .epivoi-track-label-box {
        background: #000 !important;
        color: #fff !important;
        border: 1px solid #fff !important;
      }

      body.a11y .epivoi-track-label {
        background: #000 !important;
        color: #fff !important;
        border: 1px solid #fff !important;
      }
    "))
  ),

  useShinyjs(),

  navbarPage(
    id = "main_tabs",
    title = "EpIVoi",

    tabPanel(
      "Overview",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("a11y_mode", "High-contrast mode", value = FALSE),

          h4("Metric distribution plots"),

          selectInput(
            "metric_for_dist",
            "Metric for distribution plot:",
            choices = c(
              correlation_test_choices,
              if (has_train_error) {
                stats::setNames("train_err", train_error_label)
              },
              if (has_test_error) {
                stats::setNames("test_err", test_error_label)
              }
            ),
            selected = default_test_column
          ),

          radioButtons(
            "dist_geom",
            "Distribution type:",
            choices = c(
              "Density curve" = "density",
              "Histogram" = "hist"
            ),
            selected = "density"
          ),

          

          hr(),

          h4("Correlation classes"),

          radioButtons(
            "corr_metric_for_class",
            "Correlation metric for class cut-offs:",
            choices = correlation_test_choices,
            selected = default_test_column
          ),

          numericInput(
            "low_cut",
            "Low/medium cutoff:",
            value = 0.3,
            min = -1,
            max = 1,
            step = 0.05
          ),

          numericInput(
            "high_cut",
            "Medium/high cutoff:",
            value = 0.5,
            min = -1,
            max = 1,
            step = 0.05
          ),

          tags$small(
            "Classes are based on the selected test correlation metric."
          ),

          checkboxGroupInput(
            "overview_classes",
            "Show correlation classes:",
            choices = c("failed", "low", "medium", "high"),
            selected = c("failed", "low", "medium", "high")
          )
        ),
        mainPanel(
          h3("Gene-level performance overview"),

          plotOutput("dist_plot", height = "280px"),
          downloadButton("download_dist_plot_pdf", "Download distribution plot as PDF"),

          br(),
          plotOutput("violin_plot", height = "280px"),
          downloadButton(
            "download_violin_plot_pdf",
            "Download violin plot as PDF"
          ),

          br(),
          plotOutput("class_barplot", height = "280px"),
          downloadButton("download_class_barplot_pdf", "Download class barplot as PDF"),

          br(),
          br(),
          h4("Performance table"),

          textInput(
            "overview_gene_search",
            "Search gene:",
            value = "",
            placeholder = "Gene symbol or Ensembl ID, e.g. ABHD5 or ENSG..."
          ),

          download_buttons_ui("overview", "Download Performance"),
          DTOutput("genes_table")
        )
      )
    ),
    if (has_train_error && has_test_error) {
      tabPanel(
        "Train vs Test",
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "train_test_metric",
              "Metric for train vs test scatter:",
              choices = train_test_choices,
              selected = "mse"
            ),
            checkboxInput("tt_log_mse", paste("Use log10 scale for", error_metric_name), TRUE),
            checkboxGroupInput(
              "tt_classes",
              "Show test-correlation classes:",
              choices = c("low", "medium", "high", "failed"),
              selected = c("low", "medium", "high", "failed")
            ),
            tags$small(
              paste0("For a correlation scatter, point classes use the matching test correlation. For ", error_metric_name, ", they use the selected highlighting metric.")
            ),
            selectInput(
              "highlight_metric",
              "Correlation metric for top-gene highlighting:",
              choices = correlation_test_choices,
              selected = default_test_column
            ),
            

            numericInput(
              "highlight_top_n",
              "Number of top genes to highlight (0 = none):",
              value = 50,
              min = 0,
              max = 500,
              step = 5
            )
          ),
          mainPanel(
            h3("Train vs Test performance across genes"),
            plotOutput("train_test_scatter", height = "350px"),
            br(),
            h4("Top genes among selected correlation classes"),
            download_buttons_ui("train_test", "Download Train vs Test genes"),
            DTOutput("best_genes_table")
          )
        )
      )
    },

    tabPanel(
      "Feature Importance",
      sidebarLayout(
        sidebarPanel(
          selectInput("shap_cell_type", "Cell type:", choices = app_cell_types),

          selectizeInput(
            "shap_gene",
            "Gene symbol:",
            choices = NULL,
            selected = NULL,
            multiple = FALSE,
            options = list(
              placeholder = "Type a gene symbol"
            )
          ),

          numericInput(
            "shap_top_n",
            "Show top N regions by absolute Feature Importance:",
            value = 30,
            min = 5,
            max = 200,
            step = 5
          ),

          numericInput(
            "shap_abs_min",
            "Filter regions by minimum absolute Feature Importance:",
            value = 0,
            min = 0,
            step = 0.01
          ),
          
          tags$small(
            "Filter, rank, and show regions by absolute Feature Importance. The signed Feature Importance value is still preserved in the plot and table."
          ),
        ),
        mainPanel(
          h3("Feature Importance regions for selected gene and cell type"),
          plotOutput("shap_barplot", height = "350px"),

          h4("Feature Importance table"),
          download_buttons_ui("shap", "Download Feature Importance regions"),
          DTOutput("shap_table")
        )
      )
    ),
    tabPanel(
      "IGV",
      sidebarLayout(
        sidebarPanel(
          actionButton("igv_reload", "Load"),
          br(),
          br(),
          selectizeInput(
            "igv_gene",
            "Gene:",
            choices = NULL,
            selected = NULL,
            multiple = FALSE,
            options = list(
              create = TRUE,
              placeholder = "Type gene (e.g. A2M)"
            )
          ),
          uiOutput("igv_hint"),

          selectizeInput(
            "igv_celltypes",
            "Cell types:",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(plugins = list("remove_button"), maxItems = 1000)
          ),

          textInput(
            "feature_importance_type",
            "Feature Importance type:",
            value = "SHAP"
          ),

          textInput(
            "epigenetic_signal_type",
            "Epigenetic Signal type:",
            value = "ATAC"
          ),

          
          
          actionButton("igv_back_gene", "Back to previous gene"),

          br(),
          br(),

          checkboxInput(
            "show_regulatory_track",
            "Show Ensembl regulatory build track",
            value = FALSE
          ),

          br(),

          


          
          radioButtons(
            "fimo_tf_source",
            "TF source for FIMO:",
            choices = c(
              "Use latest saved PASTAA TFs" = "saved",
              "Type TF / motif IDs manually" = "manual"
            ),
            selected = "saved"
          ),

          textAreaInput(
            "fimo_manual_tfs",
            "Manual TF names:",
            value = "",
            placeholder = "Example: TFAP2C, CTCF, JUN, FOS",
            rows = 3
          ),
          tags$small(
            "Enter up to 20 TF names separated by commas, spaces, or new lines. The app will match them to available JASPAR motifs automatically."
          ),
          
        

          tags$small(
            "FIMO uses the ENCODE candidate enhancer background."
          ),

          conditionalPanel(
            condition = "input.fimo_tf_source == 'saved'",
            numericInput(
              "fimo_top_n",
              "Top saved PASTAA TFs to show in FIMO:",
              value = 10,
              min = 1,
              max = 20,
              step = 1
            ),
            tags$small(
              "This only affects FIMO when using the latest saved PASTAA TFs. Manual TF input ignores this number."
            )
          ),

          checkboxInput(
            "fimo_replace_tracks",
            "Replace previous FIMO tracks",
            value = TRUE
          ),
          actionButton("run_fimo", "Run FIMO / show binding sites in IGV"),

          hr(),
          tags$small("Status:"),
          verbatimTextOutput("status", placeholder = TRUE),
          br(),
          br(),

          downloadButton(
            "igv_vector_pdf",
            "Download selected tracks as PDF"
          )
        ),
        mainPanel(
          tags$div(
            tags$div(
            id = "igv-container",
            style = "width:100%; height:1250px; min-height:1250px; border:1px solid #ddd; margin-bottom:40px; overflow:auto;"
            )
          )

          
            

        
        )
      )
    ),
    tabPanel(
      "TF enrichment",
      sidebarLayout(
        sidebarPanel(
          h4("TF enrichment with PASTAA"),

          radioButtons(
            "pastaa_run_mode",
            "PASTAA run mode:",
            choices = c(
              "Run per gene" = "gene",
              "Run per biological group" = "celltype"
            ),
            selected = "gene"
          ),

          conditionalPanel(
            condition = "input.pastaa_run_mode == 'gene'",
            tags$div(
              class = "well",
              tags$strong("Gene-specific analysis"),
              tags$p(
                "Choose one gene. Feature Importance regions assigned to the selected gene are ranked, deduplicated, and limited to the chosen number of top regions before running PASTAA for each selected cell type."
              ),
              selectizeInput(
                "pastaa_gene",
                "Gene:",
                choices = NULL,
                selected = NULL,
                multiple = FALSE,
                options = list(
                  create = TRUE,
                  placeholder = "Type gene, e.g. NFATC2"
                )
              )
            )
          ),

          conditionalPanel(
            condition = "input.pastaa_run_mode == 'celltype'",
            tags$div(
              class = "well",
              tags$strong("Biological-group-level analysis across genes"),
              tags$p(
                "No gene is selected. For biological-group-level analysis, raw Feature Importance values are converted to within-gene z-scores before regions are pooled across genes. The pooled regions are then ranked, deduplicated, and limited to the selected number of top regions."
              ),
              checkboxInput(
                "pastaa_filter_performance",
                "Exclude low-performing gene models",
                value = TRUE
              ),
              conditionalPanel(
                condition = "input.pastaa_filter_performance == true",
                selectInput(
                  "pastaa_performance_metric",
                  "Performance metric:",
                  choices = correlation_test_choices,
                  selected = default_test_column
                ),
                numericInput(
                  "pastaa_min_test_correlation",
                  "Minimum selected test correlation:",
                  value = 0.2,
                  min = -1,
                  max = 1,
                  step = 0.05
                )
              )
            )
          ),

          selectizeInput(
            "pastaa_celltypes",
            "Select biological groups:",
            choices = NULL,
            selected = NULL,
            multiple = TRUE,
            options = list(
              plugins = list("remove_button"),
              maxItems = 1000
            )
          ),

          numericInput(
            "pastaa_save_top_n",
            "Save top N TFs:",
            value = 10,
            min = 1,
            max = 20,
            step = 1
          ),

          numericInput(
            "pastaa_q_cutoff",
            "Max q-value / BH FDR:",
            value = 0.05,
            min = 0,
            max = 1,
            step = 0.01
          ),

          selectInput(
            "pastaa_direction",
            "Feature Importance direction:",
            choices = c(
              "Positive" = "pos",
              "Negative" = "neg",
              "Positive and negative" = "both",
              "Absolute" = "abs"
            ),
            selected = "both"
          ),

          numericInput(
            "pastaa_n_regions",
            "Number of top regions:",
            value = 500,
            min = 10,
            max = 10000,
            step = 50
          ),

      

          actionButton("run_pastaa", "Start PASTAA"),

          hr(),
          tags$small("Status:"),
          verbatimTextOutput("pastaa_status", placeholder = TRUE)
        ),

        mainPanel(
          h3("PASTAA TF enrichment results"),

          download_buttons_ui("pastaa", "Download PASTAA results"),

          uiOutput("pastaa_top_tf_barplot_title"),

          plotOutput(
            "pastaa_top_tf_barplot_current",
            height = "350px"
          ),

          fluidRow(
            column(
              6,
              downloadButton(
                "download_pastaa_top_tf_current_pdf",
                "Download this TF plot as PDF"
              )
            )
          ),

          br(),

          h4("Positive vs. Negative Feature Importance (combined, by group)"),

          plotOutput(
            "pastaa_top_tf_barplot_combined",
            height = "450px"
          ),

          fluidRow(
            column(
              6,
              downloadButton(
                "download_pastaa_top_tf_combined_pdf",
                "Download combined TF plot as PDF"
              )
            )
          ),

          br(),

          h4("Absolute Feature Importance (by group)"),

          plotOutput(
            "pastaa_abs_dotplot",
            height = "450px"
          ),

          fluidRow(
            column(
              6,
              downloadButton(
                "download_pastaa_abs_dotplot_pdf",
                "Download absolute TF dotplot as PDF"
              )
            )
          ),

          hr(),

          uiOutput("pastaa_result_tabs")
        )
      )
    ),

    tabPanel(
      "Reproducibility",
      fluidPage(
        h3("Reproducibility code"),

        p(
          "This code reproduces the current EpIVoi analysis settings and the main analysis steps outside the GUI, including raw Feature Importance loading, biological-group-level z-score calculation, region ranking, BED/FASTA preparation, TRAP, PASTAA, and BH-FDR correction."
        ),
  

        downloadButton(
          "download_repro_code",
          "Download reproducibility code"
        ),

        br(),
        br(),

        verbatimTextOutput("repro_code")
      )
    ),


    tabPanel(
      "Help",
      fluidPage(
        h3("How to interpret EpIVoi"),

        h4("General idea"),
        p(
          "This app visualizes EpIVoi results for selected genes and cell types. ",
          "It combines model performance, Feature Importance-based regulatory regions, epigenetic signals ",
          "from Epigenetic Signal data, transcription factor enrichment with PASTAA, and predicted TF binding sites from FIMO."
        ),
        p(
          "The goal is to make it easier to inspect which genomic regions contribute to the model prediction ",
          "for a gene, whether these regions are accessible in the selected cell type, and which transcription ",
          "factors may be enriched or have predicted binding sites in these regions."
        ),

        hr(),

        h4("Overview tab"),
        p(
          "The Overview tab summarizes model performance across all genes. ",
          "Each gene is assigned to a correlation class based on the selected test correlation metric."
        ),
        tags$ul(
          tags$li(strong("failed:"), " the selected test correlation is missing or could not be calculated."),
          tags$li(strong("low:"), " the selected test correlation is below the low/medium cutoff."),
          tags$li(strong("medium:"), " the selected test correlation is between the low/medium and medium/high cutoffs."),
          tags$li(strong("high:"), " the selected test correlation is above the medium/high cutoff.")
        ),
        p(
          "The default cutoffs are 0.3 and 0.5, but they can be adjusted in the sidebar. ",
          "Changing the cutoffs affects the class labels in the overview plots and tables."
        ),
        p(
          "The distribution plots can be used to inspect the overall quality of the models. ",
          "For example, genes with high test correlation are usually better candidates for detailed inspection ",
          "in the Feature Importance and IGV tabs."
        ),

        hr(),

        h4("Train vs Test tab"),
        p(
          "The Train vs Test tab compares model performance on training and test data. ",
          "This helps to identify genes where the model generalizes well and genes where the model may overfit."
        ),
        tags$ul(
          tags$li(
            strong(paste0(error_metric_name, " mode:")),
            paste(
              "compares train and test", error_metric_name,
              "values. Points far away from the diagonal may indicate different performance on train and test data."
            )
          ),
          tags$li(
            strong("Correlation mode:"),
            " compares matching train and test correlations for the selected metric."
          ),
          tags$li(
            strong("Highlighted genes:"),
            " the best genes by test correlation can be highlighted independently of the currently selected correlation classes."
          )
        ),

        hr(),

        h4("Feature Importance tab"),
        p(
          "The Feature Importance tab shows the genomic regions that contribute most strongly to the model prediction ",
          "for a selected gene and cell type."
        ),
        tags$ul(
          tags$li(
            strong("Positive Feature Importance values"),
            " indicate regions that contribute in the positive direction to the model prediction."
          ),
          tags$li(
            strong("Negative Feature Importance values"),
            " indicate regions that contribute in the negative direction to the model prediction."
          ),
          tags$li(
            strong("Absolute Feature Importance ranking"),
            " is used to identify the strongest regions regardless of direction."
          )
        ),
        p(
          "The table contains the selected regions together with their Feature Importance scores. ",
          "The barplot shows the strongest regions ranked by absolute Feature Importance score, while preserving the sign of the score."
        ),

        hr(),

        h4("IGV tab"),
        p(
          "The IGV tab displays the selected gene locus together with Feature Importance and Epigenetic Signal tracks for the selected cell types. ",
          "The RefSeq track shows gene annotations in the genomic window."
        ),
        tags$ul(
          tags$li(
            strong("Feature Importance tracks:"),
            " show model contribution scores across genomic regions. Positive and negative values indicate opposite contribution directions."
          ),
          tags$li(
            strong("Epigenetic Signal tracks:"),
            " show the epigenetic signal for the same selected cell types. These tracks help compare model-relevant regions with the measured epigenetic signal."
          ),
          tags$li(
            strong("FIMO tracks:"),
            " show predicted TF binding site positions. These tracks are positional/presence tracks and should not be interpreted as quantitative signal tracks."
          )
        ),
        p(
          "For IGV visualization, several related T-cell subtypes may be merged into the combined ",
          strong("Tcell_types"),
          " track if the individual T-cell subtypes are available in the object. ",
          "This merged track represents the average signal across the available T-cell columns."
        ),
        p(
          "The button ",
          strong("Download selected tracks as vector PDF"),
          " exports the currently selected gene and selected cell types as a vector PDF. ",
          "The PDF is useful for figures, posters, and further editing."
        ),

        hr(),

        h4("TF enrichment tab"),
        p(
          "The TF enrichment tab runs PASTAA on regions selected by Feature Importance. ",
          "PASTAA tests whether transcription factor affinity patterns are enriched among the selected regions."
        ),
        tags$ul(
          tags$li(
            strong("Run per gene:"),
            " PASTAA uses regions associated with the selected gene."
          ),
          tags$li(
            strong("Run per biological group:"),
            " no gene is selected. PASTAA pools Feature Importance regions across genes for each selected biological group. Optionally, genes can first be filtered using their selected test-correlation metric."
          ),
          tags$li(
            strong("Positive:"),
            " ranks all regions from the highest to the lowest Feature Importance score and selects the top-ranked regions."
          ),
          tags$li(
            strong("Negative:"),
            " ranks all regions from the lowest to the highest Feature Importance score and selects the top-ranked regions."
          ),
          tags$li(
            strong("Positive and negative:"),
            " performs two separate rankings over all regions: one from highest to lowest Feature Importance and one from lowest to highest."
          ),
          tags$li(
            strong("Absolute:"),
            " ranks regions by absolute Feature Importance score regardless of sign."
          )
        ),
        p(
          "The number of top regions defines how many Feature Importance regions are used as input for each PASTAA run. ",
          "If too few valid regions are available for a selected gene, cell type, or direction, the run may stop with an error message."
        ),
        p(
          "PASTAA results are ranked by q-value. The q-value is calculated using Benjamini-Hochberg FDR correction. ",
          "The top TF barplot shows the strongest enriched TFs, while the dotplot summarizes enrichment patterns ",
          "across selected cell types and Feature Importance directions."
        ),

        hr(),

        h4("FIMO binding-site tracks"),
        p(
          "FIMO can be used to scan selected genomic regions for predicted transcription factor binding sites. ",
          "TFs can either be taken from the latest saved PASTAA results or entered manually."
        ),
        tags$ul(
          tags$li(
            strong("PASTAA-derived TFs:"),
            " useful when the user wants to inspect binding sites for enriched TFs."
          ),
          tags$li(
            strong("Manual TF input:"),
            " useful when the user wants to inspect a specific TF, for example NFATC2, CTCF, JUN, or FOS."
          ),
          tags$li(
            strong("FIMO tracks:"),
            " show predicted positions of motif matches. They are not quantitative signal tracks."
          )
        ),
        p(
          "If a TF cannot be matched to a motif, the app reports this in the status message. ",
          "If FIMO finds motif matches outside the currently displayed locus, tracks may be empty in the visible window."
        ),

        hr(),

        h4("Downloads"),
        p(
          "Most tables can be downloaded as CSV or TSV files. ",
          "The IGV tab can export the selected tracks as a vector PDF. ",
          "PASTAA plots can be exported separately for use in presentations, reports, or posters."
        ),

        hr(),

        h4("Notes and limitations"),
        tags$ul(
          tags$li(
            "The availability of cell types depends on the loaded visualization object."
          ),
          tags$li(
            "Merged tracks such as Tcell_types are averages across available related cell-type columns."
          ),
          tags$li(
            "FIMO binding sites are predictions based on motif matching and should be interpreted as candidate sites."
          ),
          tags$li(
            "Epigenetic Signal and Feature Importance tracks are complementary: Epigenetic Signal describes the measured regulatory signal, while Feature Importance describes model contribution."
          ),
          tags$li(
            "The app is intended for interactive exploration and figure generation, not as a standalone statistical validation pipeline."
          )
        )
      )
    )
  ),

  tags$script(HTML("
    (function () {
      'use strict';

      window._igvBrowser = null;
      let creating = false;
      let pendingCfg = null;
      console.log('[EpIVoi IGV] custom JS loaded');
      console.log('[EpIVoi IGV] typeof igv =', typeof igv);

      if (typeof Shiny !== 'undefined' && Shiny.setInputValue) {
        Shiny.setInputValue('igv_js_loaded', {
          time: new Date().toISOString(),
          igv_type: typeof igv
        }, {priority: 'event'});
      }

      function renderError(msg) {
        console.error('[EpIVoi IGV]', msg);

        if (typeof Shiny !== 'undefined' && Shiny.setInputValue) {
          Shiny.setInputValue('igv_js_error', String(msg), {priority: 'event'});
        }

        const el = document.getElementById('igv-container');
        if (!el) return;

        el.innerHTML = '';

        const div = document.createElement('div');
        div.style.padding = '12px';
        div.style.color = '#b00020';
        div.style.whiteSpace = 'pre-wrap';
        div.textContent = msg;

        el.appendChild(div);
      }

      function hardResetContainer() {
        const oldEl = document.getElementById('igv-container');
        if (!oldEl) return null;

        const parent = oldEl.parentNode;
        if (!parent) return oldEl;

        const newEl = oldEl.cloneNode(false);
        parent.replaceChild(newEl, oldEl);

        return newEl;
      }

      async function destroyBrowser() {
        const b = window._igvBrowser;
        if (!b) return;

        try {
          if (typeof b.dispose === 'function') {
            b.dispose();
          }
        } catch (e) {}

        try {
          if (typeof igv !== 'undefined' && igv.removeBrowser) {
            const r = igv.removeBrowser(b);
            if (r && typeof r.then === 'function') {
              await r;
            }
          }
        } catch (e) {}

        window._igvBrowser = null;
      }

      function asArray(x) {
        if (!x) return [];
        if (Array.isArray(x)) return x;
        if (typeof x === 'object') return Object.values(x);
        return [];
      }

      async function removeBuiltInRefseqAll(browser) {
        if (!browser) return;

        const tracks =
          browser.tracks ||
          (browser.trackViews || [])
            .map(function (v) {
              return v.track;
            })
            .filter(Boolean);

        const toRemove = tracks.filter(function (t) {
          if (!t) return false;

          const id = String(t.id || '').toLowerCase();
          const name = String(t.name || '').toLowerCase();

          // Remove only the broad RefSeq All track. Keep our compact
          // custom RefSeq gene track and keep RefSeq Genes.
          return (
            id === 'refseqall' ||
            id === 'refseq_all' ||
            name === 'refseq all' ||
            name === 'refseqall'
          );
        });

        for (const track of toRemove) {
          try {
            const result = browser.removeTrack(track);

            if (result && typeof result.then === 'function') {
              await result;
            }

            console.log(
              '[EpIVoi IGV] removed RefSeq All track:',
              track.name
            );
          } catch (error) {
            console.warn(
              '[EpIVoi IGV] failed to remove RefSeq All:',
              error
            );
          }
        }
      }

      function hexToRgba(hex, alpha) {
        if (!hex) return 'rgba(255,255,255,' + alpha + ')';

        hex = String(hex).replace('#', '');

        if (hex.length !== 6) {
          return 'rgba(255,255,255,' + alpha + ')';
        }

        const r = parseInt(hex.substring(0, 2), 16);
        const g = parseInt(hex.substring(2, 4), 16);
        const b = parseInt(hex.substring(4, 6), 16);

        return 'rgba(' + r + ',' + g + ',' + b + ',' + alpha + ')';
      }

      

      function styleEpIVoiTrackLabels() {
        const container = document.getElementById('igv-container');
        if (!container) return;

        const candidates = container.querySelectorAll('div, span, label');

        candidates.forEach(function (el) {
          const txt = (el.textContent || '').trim();

          const isTrackLabel =
            txt.startsWith('Feature Importance') ||
            txt.startsWith('Epigenetic Signal') ||
            txt.startsWith('FIMO - ');

          if (!isTrackLabel) return;

          let trackColor = '#666666';

          const browser = window._igvBrowser;
          if (browser) {
            const trackList =
              browser.tracks ||
              (browser.trackViews || []).map(function (v) {
                return v.track;
              }).filter(Boolean);

            const hit = trackList.find(function (t) {
              return t && t.name && String(t.name) === txt;
            });

            if (hit) {
              trackColor = hit.color || hit.labelColor || '#666666';
            }
          }

          let labelBox = el.closest(
            \".igv-track-label, .igv-track-label-container, [class*='track-label'], [class*='TrackLabel']\"
          );

          if (!labelBox) {
            labelBox = el;
          }

          el.classList.add('epivoi-track-label');
          labelBox.classList.add('epivoi-track-label-box');

          [el, labelBox].forEach(function (node) {
            node.style.setProperty('font-size', '4px', 'important');
            node.style.setProperty('font-weight', '600', 'important');
            node.style.setProperty('line-height', '0.9', 'important');

            node.style.setProperty('padding', '0px 1px', 'important');
            node.style.setProperty('border-radius', '1px', 'important');

            node.style.setProperty('max-width', '90px', 'important');
            node.style.setProperty('white-space', 'nowrap', 'important');
            node.style.setProperty('overflow', 'hidden', 'important');
            node.style.setProperty('text-overflow', 'ellipsis', 'important');

            node.style.setProperty('min-height', '0', 'important');
            node.style.setProperty('height', 'auto', 'important');

            node.style.setProperty('box-shadow', 'none', 'important');
            node.style.setProperty('z-index', '9999', 'important');

            node.style.setProperty('position', 'relative', 'important');
            node.style.setProperty('display', 'inline-block', 'important');

            const txt2 = (node.textContent || '').trim();

            let trackDiv = node.closest('.igv-track-div');

            if (!trackDiv) {
              let p = node.parentElement;

              while (p && p.id !== 'igv-container') {
                const cls = String(p.className || '');

                if (
                  cls.includes('track') ||
                  cls.includes('Track')
                ) {
                  trackDiv = p;
                  break;
                }

                p = p.parentElement;
              }
            }

            node.style.setProperty('transform', 'translateY(-30px) scale(0.70)', 'important');
            node.style.setProperty('transform-origin', 'top left', 'important');
            node.style.setProperty('z-index', '9999', 'important');

            if (
              trackDiv &&
              /^(Feature Importance|Epigenetic Signal|FIMO)/.test(txt2)
            ) {
              trackDiv.style.setProperty('margin-top', '32px', 'important');
              trackDiv.style.setProperty('margin-bottom', '6px', 'important');
            }
          });

          labelBox.style.setProperty('color', trackColor, 'important');
          labelBox.style.setProperty('background-color', hexToRgba(trackColor, 0.10), 'important');
          labelBox.style.setProperty('border', '0.5px solid ' + trackColor, 'important');

          el.style.setProperty('color', trackColor, 'important');
          el.style.setProperty('background-color', 'transparent', 'important');
          el.style.setProperty('border', 'none', 'important');
        });
      }

      let epivoiLabelObserver = null;

      function startepivoiLabelObserver() {
        const container = document.getElementById('igv-container');
        if (!container) return;

        if (epivoiLabelObserver) {
          epivoiLabelObserver.disconnect();
        }

        epivoiLabelObserver = new MutationObserver(function () {
          window.requestAnimationFrame(styleEpIVoiTrackLabels);
        });

        epivoiLabelObserver.observe(container, {
          childList: true,
          subtree: true,
          characterData: true
        });
      }

      function waitForVisibleIgvContainer(maxWaitMs = 8000) {
        return new Promise(function (resolve, reject) {
          const start = Date.now();

          function check() {
            const el = document.getElementById('igv-container');

            if (
              el &&
              el.offsetParent !== null &&
              el.clientWidth > 50 &&
              el.clientHeight > 50
            ) {
              resolve(el);
              return;
            }

            if (Date.now() - start > maxWaitMs) {
              reject(new Error('IGV container is not visible or has zero size.'));
              return;
            }

            setTimeout(check, 100);
          }

          check();
        });
      }

      async function createBrowserFromCfg(cfg) {
        try {
          console.log('[EpIVoi IGV] received configuration:', cfg);

          if (typeof igv === 'undefined') {
            throw new Error('IGV library is undefined.');
          }

          if (typeof igv.createBrowser !== 'function') {
            throw new Error('igv.createBrowser is not a function.');
          }

          await waitForVisibleIgvContainer();
          await destroyBrowser();

          const el = hardResetContainer();
          if (!el) {
            throw new Error('IGV container could not be created.');
          }

          const tracks = asArray(cfg.tracks);

          console.log('[EpIVoi IGV] creating browser', {
            genome: cfg.genome,
            locus: cfg.locus,
            n_tracks: tracks.length,
            tracks: tracks
          });

          if (typeof Shiny !== 'undefined' && Shiny.setInputValue) {
            Shiny.setInputValue('igv_config_received', {
              genome: cfg.genome,
              locus: cfg.locus,
              n_tracks: tracks.length
            }, {priority: 'event'});
          }

          const browser = await igv.createBrowser(el, {
            genome: cfg.genome,
            locus: cfg.locus,
            tracks: tracks
          });

          window._igvBrowser = browser;
          startepivoiLabelObserver();

          console.log('[EpIVoi IGV] browser successfully created');

          if (typeof Shiny !== 'undefined' && Shiny.setInputValue) {
            Shiny.setInputValue('igv_browser_created', {
              locus: cfg.locus,
              n_tracks: tracks.length
            }, {priority: 'event'});
          }

          // Remove only the broad built-in RefSeq All track.
          // Keep the compact Genes / RefSeq annotation track.
          await removeBuiltInRefseqAll(browser);

          styleEpIVoiTrackLabels();
          setTimeout(styleEpIVoiTrackLabels, 300);
          setTimeout(styleEpIVoiTrackLabels, 800);
          setTimeout(styleEpIVoiTrackLabels, 1500);

        } catch (err) {
          console.error('[EpIVoi IGV] createBrowser failed:', err);

          renderError(
            'IGV createBrowser failed: ' +
            String(err && err.stack ? err.stack : err)
          );
        }
      }

      function registerEpIVoiIgvHandlers() {
        if (window._epivoiIgvHandlersRegistered) {
          return true;
        }

        if (
          typeof Shiny === 'undefined' ||
          typeof Shiny.addCustomMessageHandler !== 'function'
        ) {
          return false;
        }

        Shiny.addCustomMessageHandler('igv-create', function (cfg) {
          console.log('[EpIVoi IGV] igv-create message received', cfg);
          pendingCfg = cfg;

          if (creating) return;

          creating = true;

          (async function loop() {
            while (pendingCfg) {
              const c = pendingCfg;
              pendingCfg = null;

              await createBrowserFromCfg(c);
            }

            creating = false;
          })();
        });

        Shiny.addCustomMessageHandler('igv-search', async function (msg) {
          const b = window._igvBrowser;
          if (!b) return;

          if (msg && msg.locus && typeof b.search === 'function') {
            try {
              const r = b.search(msg.locus);
              if (r && typeof r.then === 'function') {
                await r;
              }
            } catch (e) {}
          }
        });

        Shiny.addCustomMessageHandler('igv-add-tracks', async function (msg) {

          const tracks = msg && msg.tracks ? msg.tracks : [];

          if (!tracks.length) {
            return;
          }

          // Wait until the IGV browser is actually available.
          let b = window._igvBrowser;
          let attempts = 0;
          const maxAttempts = 50;

          while (!b && attempts < maxAttempts) {
            await new Promise(function(resolve) {
              setTimeout(resolve, 100);
            });

            b = window._igvBrowser;
            attempts += 1;
          }

          if (!b) {
            console.error(
              '[EpIVoi IGV] Could not add tracks: IGV browser is not ready.',
              tracks
            );

            if (
              typeof Shiny !== 'undefined' &&
              typeof Shiny.setInputValue === 'function'
            ) {
              Shiny.setInputValue(
                'igv_track_add_error',
                {
                  error: 'IGV browser was not ready after waiting.',
                  n_tracks: tracks.length,
                  time: new Date().toISOString()
                },
                {priority: 'event'}
              );
            }

            return;
          }

          for (const tr of tracks) {
            try {
              const r = b.loadTrack(tr);

              if (r && typeof r.then === 'function') {
                await r;
              }

              console.log(
                '[EpIVoi IGV] Added track:',
                tr && tr.name ? tr.name : tr
              );

            } catch (e) {
              console.warn(
                '[EpIVoi IGV] Failed to add track',
                tr,
                e
              );

              if (
                typeof Shiny !== 'undefined' &&
                typeof Shiny.setInputValue === 'function'
              ) {
                Shiny.setInputValue(
                  'igv_track_add_error',
                  {
                    track: tr && tr.name ? tr.name : 'unknown',
                    error: String(e),
                    time: new Date().toISOString()
                  },
                  {priority: 'event'}
                );
              }
            }
          }

          setTimeout(styleEpIVoiTrackLabels, 300);
          setTimeout(styleEpIVoiTrackLabels, 800);
          setTimeout(styleEpIVoiTrackLabels, 1500);
        });

        Shiny.addCustomMessageHandler('igv-remove-fimo-tracks', async function (msg) {
          const b = window._igvBrowser;
          if (!b) return;

          const tracks =
            b.tracks ||
            (b.trackViews || []).map(function (v) {
              return v.track;
            }).filter(Boolean);

          const toRemove = tracks.filter(function (t) {
            return t && t.id && String(t.id).startsWith('fimo__');
          });

          for (const t of toRemove) {
            try {
              const r = b.removeTrack(t);
              if (r && typeof r.then === 'function') {
                await r;
              }
            } catch (e) {
              console.warn('[IGV] failed to remove FIMO track', e);
            }
          }
        });

        window._epivoiIgvHandlersRegistered = true;
        console.log('[EpIVoi IGV] Shiny message handlers registered');

        if (typeof Shiny.setInputValue === 'function') {
          Shiny.setInputValue(
            'igv_handlers_ready',
            {time: new Date().toISOString()},
            {priority: 'event'}
          );
        }

        return true;
      }

      function ensureEpIVoiIgvHandlers() {
        if (registerEpIVoiIgvHandlers()) return;
        window.setTimeout(ensureEpIVoiIgvHandlers, 100);
      }

      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', ensureEpIVoiIgvHandlers);
      } else {
        ensureEpIVoiIgvHandlers();
      }

      if (window.jQuery) {
        window.jQuery(document).on('shiny:connected', ensureEpIVoiIgvHandlers);
      }


    })();
  "))
  
)



make_download_name <- function(prefix, ext, gene = NULL, ct = NULL) {
  parts <- c(prefix)

  if (!is.null(gene) && nzchar(as.character(gene))) {
    parts <- c(parts, make_safe_id(gene))
  }

  if (!is.null(ct) && nzchar(as.character(ct))) {
    parts <- c(parts, make_safe_id(ct))
  }

  paste0(
    paste(parts, collapse = "_"),
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".",
    ext
  )
}

format_table_number <- function(x, digits = 3) {
  out <- rep(NA_character_, length(x))

  x_num <- suppressWarnings(as.numeric(x))

  is_ok <- !is.na(x_num) & is.finite(x_num)
  is_true_zero <- is_ok & x_num == 0
  is_small <- is_ok & x_num != 0 & abs(x_num) < 10^(-digits)
  is_normal <- is_ok & !is_true_zero & !is_small

  out[is_true_zero] <- formatC(0, format = "f", digits = digits)
  out[is_small] <- formatC(x_num[is_small], format = "e", digits = digits)
  out[is_normal] <- formatC(x_num[is_normal], format = "f", digits = digits)

  out
}

format_numeric_columns_for_display <- function(df, digits = 3) {
  df |>
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ format_table_number(.x, digits = digits)
      )
    )
}

write_download_table <- function(df, file, sep = ",") {
  write.table(
    df,
    file = file,
    sep = sep,
    dec = ".",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}



# ========= 3. SERVER LOGIC ================================================

server <- function(input, output, session) {
  regulatory_build_cache <- reactiveVal(list())

  observeEvent(input$a11y_mode, {
    if (isTRUE(input$a11y_mode)) {
      shinyjs::runjs("document.body.classList.add('a11y');")
    } else {
      shinyjs::runjs("document.body.classList.remove('a11y');")
    }
  }, ignoreInit = TRUE)

  overview_table_all <- reactive({
    sort_metric <- input$corr_metric_for_class %||% default_test_column

    shiny::validate(
      shiny::need(
        sort_metric %in% unname(correlation_test_choices),
        "Please select test Pearson or test Spearman."
      )
    )

    perf_overview_filtered() |>
      dplyr::mutate(
        sort_value = suppressWarnings(as.numeric(as.character(.data[[sort_metric]])))
      ) |>
      dplyr::arrange(
        dplyr::desc(sort_value),
        test_err
      ) |>
      dplyr::select(
        -sort_value
      ) |>
      dplyr::rename(
        correlation_class = corr_class
      ) |>
      dplyr::select(
        -corr_value_for_class
      ) |>
      format_numeric_columns_for_display(3)
  })

  overview_table_top <- reactive({
    sort_metric <- input$corr_metric_for_class %||% default_test_column

    perf_overview_filtered() |>
      dplyr::mutate(
        sort_value = suppressWarnings(as.numeric(as.character(.data[[sort_metric]])))
      ) |>
      dplyr::arrange(
        dplyr::desc(sort_value),
        test_err
      ) |>
      dplyr::slice_head(n = 100) |>
      dplyr::select(
        -sort_value
      ) |>
      dplyr::rename(
        correlation_class = corr_class
      ) |>
      dplyr::select(
        -corr_value_for_class
      ) |>
      format_numeric_columns_for_display(3)
  })

  # ======================================================================
  # 1. PERFORMANCE TABLE + OVERVIEW
  # ======================================================================

  perf_with_class <- reactive({
    shiny::validate(
      shiny::need(
        input$corr_metric_for_class %in% unname(correlation_test_choices),
        "Please select test Pearson or test Spearman for correlation classes."
      ),
      shiny::need(
        is.finite(input$low_cut) && is.finite(input$high_cut),
        "Correlation cutoffs must be numeric."
      ),
      shiny::need(
        input$low_cut >= -1 && input$low_cut <= 1 &&
          input$high_cut >= -1 && input$high_cut <= 1,
        "Correlation cutoffs must be between -1 and 1."
      ),
      shiny::need(
        input$high_cut > input$low_cut,
        "Medium/high cutoff must be greater than low/medium cutoff."
      )
    )

    assign_corr_class(
      perf_df,
      corr_metric = input$corr_metric_for_class,
      low_cut = input$low_cut,
      high_cut = input$high_cut
    )
  })

  tt_class_metric <- reactive({
    if (!is.null(input$train_test_metric) && input$train_test_metric %in% available_paired_correlations) {
      paste0("test_", input$train_test_metric)
    } else {
      input$highlight_metric %||% default_test_column
    }
  })

  perf_with_tt_class <- reactive({
    metric <- tt_class_metric()

    shiny::validate(
      shiny::need(
        metric %in% unname(correlation_test_choices),
        "The selected test correlation is not available."
      ),
      shiny::need(
        is.finite(input$low_cut) && is.finite(input$high_cut),
        "Correlation cutoffs must be numeric."
      ),
      shiny::need(
        input$low_cut >= -1 && input$low_cut <= 1 &&
          input$high_cut >= -1 && input$high_cut <= 1,
        "Correlation cutoffs must be between -1 and 1."
      ),
      shiny::need(
        input$high_cut > input$low_cut,
        "Medium/high cutoff must be greater than low/medium cutoff."
      )
    )

    assign_corr_class(
      perf_df,
      corr_metric = metric,
      low_cut = input$low_cut,
      high_cut = input$high_cut
    )
  })

  perf_overview_filtered <- reactive({
    df <- perf_with_class()

    selected_classes <- input$overview_classes

    if (is.null(selected_classes) || length(selected_classes) == 0) {
      return(df[0, , drop = FALSE])
    }

    df |>
      dplyr::filter(corr_class %in% selected_classes)
  })

  corr_colors_active <- reactive({
    if (isTRUE(input$a11y_mode)) corr_colors_a11y else corr_colors
  })

  legend_title_for_corr_classes <- reactive({
    paste0(
      "Correlation class\n",
      "Metric: ", input$corr_metric_for_class, "\n",
      "low < ", input$low_cut, "\n",
      "medium: ", input$low_cut, "–", input$high_cut, "\n",
      "high ≥ ", input$high_cut
    )
  })

  legend_title_for_tt_classes <- reactive({
    paste0(
      paste0(sub("^test_", "Test ", tt_class_metric()), " class\n"),
      "low < ", input$low_cut, "\n",
      "medium: ", input$low_cut, "–", input$high_cut, "\n",
      "high ≥ ", input$high_cut
    )
  })

  make_dist_plot <- function() {
    df <- perf_overview_filtered()
    metric <- input$metric_for_dist
    req(metric %in% names(df))

    g <- ggplot(df, aes(x = .data[[metric]]))

    if (input$dist_geom == "density") {
      g <- g + geom_density(na.rm = TRUE, fill = "steelblue", alpha = 0.4)
    } else {
      g <- g + geom_histogram(na.rm = TRUE, bins = 60, fill = "steelblue", alpha = 0.8)
    }

    g +
      labs(
        x = metric,
        y = ifelse(input$dist_geom == "density", "Density", "Count"),
        title = paste("Distribution of", metric, "across genes")
      ) +
      theme_bw(base_size = 14)
  }


  make_violin_plot <- function() {
    df_all <- perf_overview_filtered()
    metric <- input$metric_for_dist
    req(metric %in% names(df_all))

    df_plot <- df_all |>
      dplyr::filter(
        is.finite(suppressWarnings(as.numeric(.data[[metric]]))),
        corr_class != "failed"
      )

    failed_n <- sum(df_all$corr_class == "failed", na.rm = TRUE)

    p <- ggplot(
      df_plot,
      aes(
        x = corr_class,
        y = .data[[metric]],
        fill = corr_class
      )
    ) +
      geom_violin(
        na.rm = TRUE,
        alpha = 0.7
      ) +
      scale_x_discrete(
        limits = corr_levels,
        drop = FALSE
      ) +
      scale_fill_manual(
        values = corr_colors_active(),
        drop = FALSE
      ) +
      labs(
        x = paste("Class based on", input$corr_metric_for_class),
        y = metric,
        title = paste("Distribution of", metric, "by correlation class"),
        fill = "Correlation class"
      ) +
      theme_bw(base_size = 14)

    if (failed_n > 0) {
      y_min <- min(df_plot[[metric]], na.rm = TRUE)

      p <- p +
        geom_point(
          data = data.frame(
            corr_class = factor("failed", levels = corr_levels),
            plot_value = y_min
          ),
          aes(
            x = corr_class,
            y = plot_value
          ),
          inherit.aes = FALSE,
          size = 3
        ) +
        annotate(
          "text",
          x = "failed",
          y = y_min,
          label = paste0("n = ", failed_n),
          vjust = -1
        )
    }

    p
  }


  make_class_barplot <- function() {
    df <- perf_overview_filtered() |>
      dplyr::count(corr_class)

    ggplot(
      df,
      aes(
        x = corr_class,
        y = n,
        fill = corr_class
      )
    ) +
      geom_col(alpha = 0.8) +
      scale_fill_manual(values = corr_colors_active(), drop = FALSE) +
      labs(
        x = paste("Class based on", input$corr_metric_for_class),
        y = "Number of genes",
        title = paste("Number of genes per", input$corr_metric_for_class, "correlation class"),
        fill = "Correlation class"
      ) +
      theme_bw(base_size = 14)
  }

  output$dist_plot <- renderPlot({
    make_dist_plot()
  })

  output$violin_plot <- renderPlot({
    make_violin_plot()
  })

  output$class_barplot <- renderPlot({
    make_class_barplot()
  })
  save_plot_pdf <- function(file, plot_fun, width = 8, height = 5) {
    p <- plot_fun()

    grDevices::cairo_pdf(
      filename = file,
      width = width,
      height = height,
      onefile = TRUE
    )

    print(p)

    grDevices::dev.off()
  }

  output$download_dist_plot_pdf <- downloadHandler(
    filename = function() {
      make_download_name("overview_distribution_plot", "pdf")
    },
    content = function(file) {
      save_plot_pdf(file, make_dist_plot, width = 8, height = 5)
    }
  )

  output$download_violin_plot_pdf <- downloadHandler(
    filename = function() {
      make_download_name("overview_violin_plot", "pdf")
    },
    content = function(file) {
      save_plot_pdf(file, make_violin_plot, width = 8, height = 5)
    }
  )

  output$download_class_barplot_pdf <- downloadHandler(
    filename = function() {
      make_download_name("overview_class_barplot", "pdf")
    },
    content = function(file) {
      save_plot_pdf(file, make_class_barplot, width = 8, height = 5)
    }
  )

  output$genes_table <- renderDT({

    df <- overview_table_all()

    search_value <- trimws(input$overview_gene_search %||% "")

    if (nzchar(search_value)) {

      search_clean <- strip_ensembl_version(search_value)
      search_upper <- toupper(search_value)

      # If the user entered an Ensembl gene ID, resolve it to the gene symbol.
      resolved_symbol <- if (search_clean %in% names(ensembl_to_symbol_map)) {
        unname(ensembl_to_symbol_map[[search_clean]])
      } else {
        search_value
      }

      df <- df |>
        dplyr::filter(
          grepl(
            toupper(resolved_symbol),
            toupper(as.character(gene_symbol)),
            fixed = TRUE
          )
        )
    }

    datatable(
      df,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        ordering = TRUE,
        dom = "ltip"
      ),
      rownames = FALSE
    )
  })


  # ======================================================================
  # 2. TRAIN vs TEST
  # ======================================================================

  tt_df <- reactive({
    df <- perf_with_tt_class()

    if (!is.null(input$tt_classes)) {
      df <- df |> dplyr::filter(corr_class %in% input$tt_classes)
    }

    df
  })

  best_genes_symbols <- reactive({
    if (is.null(input$highlight_top_n) || input$highlight_top_n <= 0) {
      return(character(0))
    }

    perf_with_tt_class() |>
      dplyr::mutate(
        top_gene_metric_value = suppressWarnings(
          as.numeric(
            as.character(
              .data[[input$highlight_metric %||% default_test_column]]
            )
          )
        )
      ) |>
      dplyr::filter(
        !is.na(top_gene_metric_value),
        is.finite(top_gene_metric_value)
      ) |>
      dplyr::arrange(
        dplyr::desc(top_gene_metric_value),
        test_err
      ) |>
      dplyr::slice_head(n = input$highlight_top_n) |>
      dplyr::pull(gene_symbol)
  })

  highlight_genes_df <- reactive({
    top_symbols <- best_genes_symbols()

    perf_with_tt_class() |>
      dplyr::filter(gene_symbol %in% top_symbols) |>
      dplyr::mutate(
        top_gene_metric_value = suppressWarnings(as.numeric(as.character(.data[[input$highlight_metric %||% default_test_column]])))
      )
  })

  output$train_test_scatter <- renderPlot({
    df <- tt_df()

    pt_alpha <- if (isTRUE(input$a11y_mode)) 0.9 else 0.5
    pt_size  <- if (isTRUE(input$a11y_mode)) 2.6 else 1.6

    shape_map <- c(
      failed = 4,
      low = 16,
      medium = 17,
      high = 15
    )

    # ============================================================
    # Error metric: train_err vs test_err
    # ============================================================
    if (input$train_test_metric == "mse") {
      req(all(c("train_err", "test_err") %in% names(df)))

      if (input$tt_log_mse) {
        df <- df |>
          dplyr::mutate(
            train_err_plot = ifelse(train_err > 0, log10(train_err), NA),
            test_err_plot  = ifelse(test_err > 0, log10(test_err), NA)
          )

        x_col <- "train_err_plot"
        y_col <- "test_err_plot"
        x_lab <- train_error_log_label
        y_lab <- test_error_log_label
      } else {
        x_col <- "train_err"
        y_col <- "test_err"
        x_lab <- train_error_label
        y_lab <- test_error_label
      }

      g <- ggplot(
        df,
        aes(
          x = .data[[x_col]],
          y = .data[[y_col]],
          colour = corr_class,
          shape = corr_class
        )
      ) +
        geom_point(alpha = pt_alpha, size = pt_size, na.rm = TRUE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        scale_color_manual(values = corr_colors_active(), drop = FALSE) +
        scale_shape_manual(values = shape_map, drop = FALSE) +
        labs(
          title = train_test_error_title,
          x = x_lab,
          y = y_lab,
          colour = legend_title_for_tt_classes(),
          shape = legend_title_for_tt_classes()
        ) +
        theme_bw(base_size = 14)

      # Highlight top N genes ALWAYS from full data,
      # even if their correlation class is not selected in tt_classes.
      if (!is.null(input$highlight_top_n) && input$highlight_top_n > 0) {
        hdf <- highlight_genes_df()

        if (input$tt_log_mse) {
          hdf <- hdf |>
            dplyr::mutate(
              train_err_plot = ifelse(train_err > 0, log10(train_err), NA),
              test_err_plot  = ifelse(test_err > 0, log10(test_err), NA)
            )
        }

        g <- g +
          geom_point(
            data = hdf,
            mapping = aes(
              x = .data[[x_col]],
              y = .data[[y_col]]
            ),
            inherit.aes = FALSE,
            size = if (isTRUE(input$a11y_mode)) 3.6 else 2.6,
            colour = "red",
            shape = 21,
            fill = NA,
            stroke = 1.4,
            na.rm = TRUE
          )
      }

      g

    # ============================================================
    # Correlation: train correlation vs selected test correlation
    # ============================================================
    } else {
      corr_type <- input$train_test_metric
      req(corr_type %in% available_paired_correlations)

      train_corr_metric <- paste0("train_", corr_type)
      test_corr_metric <- paste0("test_", corr_type)

      req(all(c(train_corr_metric, test_corr_metric) %in% names(df)))

      g <- ggplot(
        df,
        aes(
          x = .data[[train_corr_metric]],
          y = .data[[test_corr_metric]],
          colour = corr_class,
          shape = corr_class
        )
      ) +
        geom_point(alpha = pt_alpha, size = pt_size, na.rm = TRUE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        scale_color_manual(values = corr_colors_active(), drop = FALSE) +
        scale_shape_manual(values = shape_map, drop = FALSE) +
        labs(
          title = paste("Train vs Test", corr_type, "correlation"),
          x = paste("Train", corr_type, "correlation"),
          y = paste("Test", corr_type, "correlation"),
          colour = legend_title_for_tt_classes(),
          shape = legend_title_for_tt_classes()
        ) +
        theme_bw(base_size = 14)

      # Highlight top N genes ALWAYS from full data,
      # even if their correlation class is not selected in tt_classes.
      if (!is.null(input$highlight_top_n) && input$highlight_top_n > 0) {
        hdf <- highlight_genes_df()

        g <- g +
          geom_point(
            data = hdf,
            mapping = aes(
              x = .data[[train_corr_metric]],
              y = .data[[test_corr_metric]]
            ),
            inherit.aes = FALSE,
            size = if (isTRUE(input$a11y_mode)) 3.6 else 2.6,
            colour = "red",
            shape = 21,
            fill = NA,
            stroke = 1.4,
            na.rm = TRUE
          )
      }

      g
    }
  })

  train_test_table_top <- reactive({
    if (is.null(input$highlight_top_n) || input$highlight_top_n <= 0) {
      return(data.frame())
    }
    metric <- input$highlight_metric %||% default_test_column

    tt_df() |>
      dplyr::mutate(
        table_metric_value = suppressWarnings(as.numeric(as.character(.data[[metric]])))
      ) |>
      dplyr::filter(is.finite(table_metric_value)) |>
      dplyr::arrange(dplyr::desc(table_metric_value), test_err) |>
      dplyr::slice_head(n = input$highlight_top_n) |>
      dplyr::select(-table_metric_value) |>
      dplyr::rename(correlation_class = corr_class) |>
      dplyr::select(-corr_value_for_class) |>
      format_numeric_columns_for_display(3)
  })

  train_test_table_all <- reactive({
    tt_df() |>
      dplyr::rename(correlation_class = corr_class) |>
      dplyr::select(-corr_value_for_class) |>
      format_numeric_columns_for_display(3)
  })

  output$best_genes_table <- renderDT({
    datatable(
      train_test_table_top(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        ordering = TRUE
      ),
      rownames = FALSE
    )
  })

  observeEvent(input$igv_js_error, {
    status_msg(paste("IGV JavaScript error:", input$igv_js_error))
  })

  # ======================================================================
  # 3. SHAP
  # ======================================================================

  observeEvent(TRUE, {
    genes <- igv_gene_list
    genes <- genes[!is.na(genes) & nzchar(genes)]

    if (length(genes) > 0) {
      updateSelectizeInput(
        session,
        "shap_gene",
        choices = genes,
        selected = genes[1],
        server = TRUE
      )
    }
  }, once = TRUE)

  
  shap_top_current <- reactive({
    shiny::validate(
      shiny::need(
        is.finite(input$shap_top_n) && input$shap_top_n >= 1,
        "Top N Feature Importance regions must be at least 1."
      ),
      shiny::need(
        is.finite(input$shap_abs_min) && input$shap_abs_min >= 0,
        "Minimum absolute Feature Importance score must be >= 0."
      )
    )

    ct <- input$shap_cell_type
    gene_input <- input$shap_gene

    gene_resolved <- resolve_igv_gene_input(gene_input)

    shiny::validate(
      shiny::need(!is.null(gene_resolved), "Selected gene was not found.")
    )

    shap_df <- load_shap_df(gene_resolved$gene_id, ct)

    shap_df <- shap_df |>
      dplyr::mutate(
        abs_score = abs(score)
      )

    if (!is.null(input$shap_abs_min) && input$shap_abs_min > 0) {
      shap_df <- shap_df |>
        dplyr::filter(abs_score >= input$shap_abs_min)
    }

    shap_df |>
      dplyr::arrange(dplyr::desc(abs_score)) |>
      dplyr::slice_head(n = input$shap_top_n)
  })

  shap_all_current <- reactive({
    ct <- input$shap_cell_type
    gene_input <- input$shap_gene

    gene_resolved <- resolve_igv_gene_input(gene_input)

    shiny::validate(
      shiny::need(!is.null(gene_resolved), "Selected gene was not found.")
    )

    shap_df <- load_shap_df(gene_resolved$gene_id, ct)

    shap_df <- shap_df |>
      dplyr::mutate(
        abs_score = abs(score)
      )

    if (!is.null(input$shap_abs_min) && input$shap_abs_min > 0) {
      shap_df <- shap_df |>
        dplyr::filter(abs_score >= input$shap_abs_min)
    }

    shap_df |>
      dplyr::arrange(dplyr::desc(abs_score))
  })

  

  output$shap_table <- renderDT({
    df <- shap_top_current() |>
      dplyr::select(
        gene_id,
        gene_symbol,
        interaction,
        score,
        abs_score
      ) |>
      dplyr::mutate(
        score = round(score, 3),
        abs_score = round(abs_score, 3)
      )

    datatable(
      df,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        ordering = TRUE
      ),
      rownames = FALSE
    )
  })
  
  output$shap_barplot <- renderPlot({
    df <- shap_top_current()
    req(nrow(df) > 0)

    df <- df |>
      dplyr::mutate(
        interaction_short = stringr::str_trunc(as.character(interaction), width = 30)
      )

    ggplot(df, aes(
      x = reorder(interaction_short, score),
      y = score,
      fill = score
    )) +
      geom_col() +
      scale_fill_gradient2(
        low = "red",
        mid = "grey90",
        high = "blue",
        midpoint = 0
      ) +
      labs(
        x = "Region, ranked by absolute Feature Importance score",
        y = "Signed Feature Importance score",
        fill = "Feature Importance",
        title = paste(
          "Feature Importance scores for",
          input$shap_gene,
          "in",
          input$shap_cell_type
        ),
        subtitle = paste0(
          "Showing top ",
          input$shap_top_n,
          " regions after absolute Feature Importance filter ≥ ",
          input$shap_abs_min
        )
      ) +
      theme_bw(base_size = 14) +
      coord_flip()
  })


  # ======================================================================
  # PASTAA helper functions
  # ======================================================================

  load_meme_motif_map <- function(meme_file) {
    if (!file.exists(meme_file)) {
      stop("MEME motif file not found: ", meme_file)
    }

    lines <- readLines(meme_file, warn = FALSE)
    motif_lines <- grep("^MOTIF\\s+", lines, value = TRUE)

    if (length(motif_lines) == 0) {
      return(data.frame(
        motif_id = character(0),
        tf_name = character(0),
        stringsAsFactors = FALSE
      ))
    }

    x <- strsplit(motif_lines, "\\s+")

    motif_id <- vapply(x, function(z) {
      if (length(z) >= 2) z[2] else NA_character_
    }, character(1))

    tf_name <- vapply(x, function(z) {
      if (length(z) >= 3) {
        paste(z[3:length(z)], collapse = "_")
      } else {
        z[2]
      }
    }, character(1))

    data.frame(
      motif_id = motif_id,
      tf_name = tf_name,
      stringsAsFactors = FALSE
    ) |>
      dplyr::filter(!is.na(motif_id), nzchar(motif_id)) |>
      dplyr::distinct(motif_id, .keep_all = TRUE)
  }

  clean_tf_name <- function(x) {
    x <- as.character(x)
    x <- stringr::str_replace(x, "^#\\s*", "")
    x <- stringr::str_replace(x, "\\s*-\\s*$", "")
    x <- stringr::str_trim(x)
    x
  }

  tf_display_name <- function(x) {
    x <- clean_tf_name(x)

    # Example: RORA(MA0071.1) -> RORA
    x <- ifelse(
      grepl("\\(MA[0-9]{4}\\.[0-9]+\\)$", x),
      sub("\\(MA[0-9]{4}\\.[0-9]+\\)$", "", x),
      x
    )

    x <- stringr::str_trim(x)
    x
  }

  clean_tf_key <- function(x) {
    x <- clean_tf_name(x)
    x <- toupper(x)
    gsub("[^A-Z0-9]+", "", x)
  }

  map_tf_or_motif_to_motif_ids <- function(tf_values, motif_map) {
    tf_values <- unique(as.character(tf_values))
    tf_values <- trimws(tf_values)
    tf_values <- tf_values[nzchar(tf_values)]

    if (length(tf_values) == 0) {
      return(list(
        motif_ids = character(0),
        missing_tfs = character(0),
        mapping_table = data.frame()
      ))
    }

    motif_map2 <- motif_map |>
      dplyr::mutate(
        tf_name = as.character(tf_name),
        motif_id = as.character(motif_id),
        tf_key = clean_tf_key(tf_name)
      )

    motif_ids <- character(0)
    missing_tfs <- character(0)
    mapping_rows <- list()

    for (tf in tf_values) {
      tf_clean <- clean_tf_name(tf)

      direct_ma <- stringr::str_extract(tf_clean, "MA[0-9]{4}\\.[0-9]+")

      if (!is.na(direct_ma) && nzchar(direct_ma)) {
        motif_ids <- c(motif_ids, direct_ma)

        mapping_rows[[length(mapping_rows) + 1]] <- data.frame(
          TF_input = tf,
          TF_clean = tf_clean,
          motif_id = direct_ma,
          stringsAsFactors = FALSE
        )

        next
      }

      key <- clean_tf_key(tf_clean)

      hit <- motif_map2 |>
        dplyr::filter(tf_key == key)

      if (nrow(hit) == 0) {
        missing_tfs <- c(missing_tfs, tf_clean)
      } else {
        motif_ids <- c(motif_ids, hit$motif_id)

        mapping_rows[[length(mapping_rows) + 1]] <- data.frame(
          TF_input = tf,
          TF_clean = tf_clean,
          motif_id = hit$motif_id,
          stringsAsFactors = FALSE
        )
      }
    }

    mapping_table <- dplyr::bind_rows(mapping_rows)

    list(
      motif_ids = unique(motif_ids),
      missing_tfs = unique(missing_tfs),
      mapping_table = mapping_table
    )
  }

  pastaa_results <- reactiveVal(list())

  

  save_top_pastaa_tfs <- function(res_list, top_n = 10, q_cutoff = 0.05) {
    if (length(res_list) == 0) return(NULL)

    top_n <- min(as.integer(top_n), 20)

    if (is.na(top_n) || top_n < 1) {
      stop("Save top N TFs must be between 1 and 20.")
    }

    df <- dplyr::bind_rows(
      lapply(names(res_list), function(nm) {
        x <- res_list[[nm]]$result
        x$run_id <- nm
        x
      })
    )

    if (nrow(df) == 0) return(NULL)

    top_df <- df |>
      dplyr::filter(TF != "DUMMY_MOTIF") |>
      dplyr::mutate(
        TF = as.character(TF),
        TF_clean = clean_tf_name(TF),
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR))
      ) |>
      dplyr::filter(
        is.finite(BH_FDR),
        BH_FDR <= q_cutoff,
        !is.na(TF_clean),
        nzchar(TF_clean)
      ) |>
      dplyr::group_by(TF_clean) |>
      dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(BH_FDR, p_value) |>
      dplyr::slice_head(n = top_n)

    if (nrow(top_df) == 0) return(NULL)

    save_file <- file.path(
      pastaa_out_dir,
      paste0("saved_top_pastaa_TFs_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    )

    write.csv(top_df, save_file, row.names = FALSE)

    save_file
  }



  get_fimo_motif_ids_from_input <- function(top_n = 10,
                                          manual_text = "",
                                          source = "saved",
                                          current_gene = NULL) {
    motif_map <- load_meme_motif_map(MEME_MOTIF_FILE)

    top_n <- min(as.integer(top_n), 20)

    if (is.na(top_n) || top_n < 1) {
      stop("Top N for FIMO must be between 1 and 20.")
    }

    if (source == "manual") {
      ids <- unlist(strsplit(manual_text, "[,;\\n\\t ]+"))
      ids <- trimws(ids)
      ids <- ids[nzchar(ids)]
      ids <- unique(ids)

      if (length(ids) == 0) {
        stop("No TF names provided. Please enter TF names such as TFAP2C, CTCF, JUN.")
      }

      if (length(ids) > 20) {
        stop("Too many TF names provided for FIMO. Please enter at most 20 TFs.")
      }

      mapped <- map_tf_or_motif_to_motif_ids(ids, motif_map)

      if (length(mapped$missing_tfs) > 0) {
        stop(
          "No JASPAR motifs found for TF name(s): ",
          paste(mapped$missing_tfs, collapse = ", "),
          ". Please check the spelling or try a different TF name."
        )
      }

      if (length(mapped$motif_ids) == 0) {
        stop("No valid motifs found for the provided TF names.")
      }

      motif_label_map <- mapped$mapping_table |>
      dplyr::mutate(
        motif_id = as.character(motif_id),
        TF_label = tf_display_name(TF_clean)
      ) |>
      dplyr::filter(
        !is.na(motif_id),
        nzchar(motif_id),
        !is.na(TF_label),
        nzchar(TF_label)
      ) |>
      dplyr::distinct(motif_id, .keep_all = TRUE)

      motif_label_map <- stats::setNames(
        motif_label_map$TF_label,
        motif_label_map$motif_id
      )

      return(list(
        motif_ids = mapped$motif_ids,
        mapping_table = mapped$mapping_table,
        source_table = data.frame(
          TF = ids,
          TF_clean = tf_display_name(ids),
          stringsAsFactors = FALSE
        ),
        missing_tfs = mapped$missing_tfs,
        motif_label_map = motif_label_map
      ))
    }

    saved_file <- get_latest_saved_pastaa_tfs(
      current_gene = current_gene
    )

    if (is.null(saved_file)) {
      stop("No saved PASTAA TF file found. Run PASTAA first.")
    }

    df <- read.csv(saved_file, stringsAsFactors = FALSE)

    required_cols <- c("TF", "gene_id", "cell_type", "direction", "BH_FDR")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
      stop(
        "Saved PASTAA file is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    if (!is.null(current_gene) && nzchar(current_gene)) {
      genes_in_file <- unique(as.character(df$gene_id))

      if (!current_gene %in% genes_in_file && !"ALL_GENES" %in% genes_in_file) {
        stop(
          "Latest saved PASTAA TF file is not for the current IGV gene. ",
          "Current IGV gene: ", current_gene,
          ". Saved file genes: ", paste(head(genes_in_file, 10), collapse = ", "),
          if (length(genes_in_file) > 10) " ..." else "",
          ". Run PASTAA for this gene first or use manual TF input."
        )
      }

      df <- df |>
        dplyr::filter(gene_id == current_gene | gene_id == "ALL_GENES")
    }

    df <- df |>
      dplyr::mutate(
        TF = as.character(TF),
        TF_clean = clean_tf_name(TF),
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR))
      ) |>
      dplyr::filter(
        !is.na(TF_clean),
        nzchar(TF_clean),
        is.finite(BH_FDR)
      ) |>
      dplyr::arrange(BH_FDR, p_value) |>
      dplyr::distinct(TF_clean, .keep_all = TRUE) |>
      dplyr::slice_head(n = top_n)

    if (nrow(df) == 0) {
      stop("No saved PASTAA TFs passed the current filtering.")
    }

    mapped <- map_tf_or_motif_to_motif_ids(df$TF_clean, motif_map)

    message(
      "FIMO selected TFs from saved PASTAA: ",
      paste(df$TF_clean, collapse = ", ")
    )

    message(
      "FIMO mapped motif IDs: ",
      if (length(mapped$motif_ids) > 0) paste(mapped$motif_ids, collapse = ", ") else "NONE"
    )

    message(
      "FIMO missing TFs: ",
      if (length(mapped$missing_tfs) > 0) paste(mapped$missing_tfs, collapse = ", ") else "NONE"
    )

    if (length(mapped$motif_ids) == 0) {
      stop(
        "None of the selected PASTAA TFs could be mapped to JASPAR motifs. ",
        "Selected TFs: ",
        paste(df$TF_clean, collapse = ", ")
      )
    }

    motif_label_map <- mapped$mapping_table |>
      dplyr::mutate(
        motif_id = as.character(motif_id),
        TF_label = tf_display_name(TF_clean)
      ) |>
      dplyr::filter(
        !is.na(motif_id),
        nzchar(motif_id),
        !is.na(TF_label),
        nzchar(TF_label)
      ) |>
      dplyr::distinct(motif_id, .keep_all = TRUE)

    motif_label_map <- stats::setNames(
      motif_label_map$TF_label,
      motif_label_map$motif_id
    )

    list(
      motif_ids = mapped$motif_ids,
      mapping_table = mapped$mapping_table,
      source_table = df |>
        dplyr::mutate(TF_clean = tf_display_name(TF_clean)),
      missing_tfs = mapped$missing_tfs,
      motif_label_map = motif_label_map
    )
  }

  get_latest_saved_pastaa_combinations <- function(current_gene = NULL) {
    saved_file <- get_latest_saved_pastaa_tfs(
      current_gene = current_gene
    )

    if (is.null(saved_file)) {
      stop("No saved PASTAA TF file found. Run PASTAA first.")
    }

    df <- read.csv(saved_file, stringsAsFactors = FALSE)

    required_cols <- c("gene_id", "cell_type", "direction", "TF", "BH_FDR")
    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
      stop(
        "Saved PASTAA file is missing columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    if (!is.null(current_gene) && nzchar(current_gene)) {
      genes_in_file <- unique(as.character(df$gene_id))

      if (!current_gene %in% genes_in_file && !"ALL_GENES" %in% genes_in_file) {
        stop(
          "Latest saved PASTAA TF file is not for the current IGV gene. ",
          "Current IGV gene: ", current_gene,
          ". Saved file genes: ", paste(head(genes_in_file, 10), collapse = ", "),
          if (length(genes_in_file) > 10) " ..." else "",
          ". Run PASTAA for this gene first or use manual TF input."
        )
      }

      df <- df |>
        dplyr::filter(gene_id == current_gene | gene_id == "ALL_GENES")
    }

    df |>
      dplyr::distinct(gene_id, cell_type)
  }

  get_latest_saved_pastaa_tfs <- function(current_gene = NULL) {
    files <- list.files(
      pastaa_out_dir,
      pattern = "^saved_top_pastaa_TFs_.*\\.csv$",
      full.names = TRUE
    )

    if (length(files) == 0) {
      return(NULL)
    }

    # Newest files first
    files <- files[
      order(file.info(files)$mtime, decreasing = TRUE)
    ]

    # If no gene was requested, keep the old behaviour
    if (
      is.null(current_gene) ||
      length(current_gene) == 0 ||
      is.na(current_gene[1]) ||
      !nzchar(as.character(current_gene[1]))
    ) {
      return(files[1])
    }

    current_gene <- as.character(current_gene[1])

    # Find the newest saved PASTAA file that actually belongs
    # to the current IGV gene.
    for (f in files) {
      df <- tryCatch(
        read.csv(f, stringsAsFactors = FALSE),
        error = function(e) NULL
      )

      if (
        is.null(df) ||
        nrow(df) == 0 ||
        !"gene_id" %in% colnames(df)
      ) {
        next
      }

      genes_in_file <- unique(as.character(df$gene_id))

      if (
        current_gene %in% genes_in_file ||
        "ALL_GENES" %in% genes_in_file
      ) {
        return(f)
      }
    }

    NULL
  }

  get_pastaa_genes_by_performance <- function(
    filter_performance = TRUE,
    performance_metric = default_test_column,
    min_test_correlation = 0.2
  ) {
    if (!isTRUE(filter_performance)) {
      return(genes_from_obj)
    }

    threshold <- suppressWarnings(as.numeric(min_test_correlation)[1])

    if (!is.finite(threshold) || threshold < -1 || threshold > 1) {
      stop("Minimum test correlation must be between -1 and 1.")
    }

    if (!performance_metric %in% unname(correlation_test_choices)) {
      stop("The selected test-correlation metric is not available: ", performance_metric)
    }

    perf <- corr_df |>
      dplyr::transmute(
        performance_gene = strip_ensembl_version(as.character(gene_symbol)),
        test_correlation = suppressWarnings(as.numeric(as.character(.data[[performance_metric]])))
      ) |>
      dplyr::filter(is.finite(test_correlation), test_correlation >= threshold)

    direct_ids <- perf$performance_gene[
      perf$performance_gene %in% genes_from_obj
    ]

    symbol_ids <- unname(
      symbol2id[toupper(perf$performance_gene)]
    )

    eligible <- unique(c(direct_ids, symbol_ids))
    eligible <- eligible[
      !is.na(eligible) & eligible %in% genes_from_obj
    ]

    if (length(eligible) == 0) {
      stop(
        "No genes passed the selected ", performance_metric, " threshold (>= ",
        threshold,
        ")."
      )
    }

    eligible
  }

  # Cache completed per-cell-type region rankings within the current Shiny session.
  # This avoids reopening thousands of per-gene SHAP files when the same
  # cell type, direction, threshold and top-N setting are requested again.
  pastaa_region_cache <- new.env(parent = emptyenv())

  make_pastaa_cache_key <- function(
    ct,
    direction,
    n_regions,
    use_abs_shap,
    filter_performance,
    performance_metric,
    min_test_correlation
  ) {
    paste(
      make_safe_id(ct),
      direction,
      as.integer(n_regions),
      isTRUE(use_abs_shap),
      isTRUE(filter_performance),
      make_safe_id(performance_metric),
      format(as.numeric(min_test_correlation), scientific = FALSE, trim = TRUE),
      sep = "__"
    )
  }

  prepare_pastaa_regions <- function(gene = NULL, ct,
                                   direction = c("pos", "neg", "abs"),
                                   n_regions = 500,
                                   run_mode = c("gene", "celltype"),
                                   use_abs_shap = FALSE,
                                   filter_performance = TRUE,
                                   performance_metric = default_test_column,
                                   min_test_correlation = 0.2) {
    direction <- match.arg(direction)
    run_mode <- match.arg(run_mode)

    cache_key <- NULL

    if (run_mode == "gene") {
      df <- load_igv_shap_df(gene, ct) |>
        dplyr::filter(is.finite(score))
    } else {
      cache_key <- make_pastaa_cache_key(
        ct = ct,
        direction = direction,
        n_regions = n_regions,
        use_abs_shap = use_abs_shap,
        filter_performance = filter_performance,
        performance_metric = performance_metric,
        min_test_correlation = min_test_correlation
      )

      if (exists(cache_key, envir = pastaa_region_cache, inherits = FALSE)) {
        message("PASTAA cell-type mode: using cached region ranking for ", ct)
        return(get(cache_key, envir = pastaa_region_cache, inherits = FALSE))
      }

      # IMPORTANT: performance filtering happens before any SHAP files are read.
      genes_for_enrichment <- get_pastaa_genes_by_performance(
        filter_performance = filter_performance,
        performance_metric = performance_metric,
        min_test_correlation = min_test_correlation
      )

      message(
        "PASTAA cell-type mode: reading Feature Importance files for ",
        length(genes_for_enrichment),
        " of ",
        length(genes_from_obj),
        " genes",
        if (isTRUE(filter_performance)) {
          paste0(" with ", performance_metric, " >= ", min_test_correlation)
        } else {
          " without performance filtering"
        }
      )

      pastaa_n_cores <- {
        env_cores <- suppressWarnings(
          as.integer(
            Sys.getenv(
              "EPIVOI_PASTAA_CORES",
              ""
            )
          )
        )

        if (
          is.finite(env_cores) &&
          env_cores >= 1L
        ) {
          env_cores
        } else {
          max(
            1L,
            floor(
              parallel::detectCores(
                logical = TRUE
              ) / 2L
            )
          )
        }
      }

      use_parallel <-
        pastaa_n_cores > 1L &&
        .Platform$OS.type == "unix"


      read_one_gene_fi <- function(g) {
        tryCatch(
          {
            gene_df <- load_igv_shap_df(
              g,
              ct
            ) |>
              dplyr::filter(
                is.finite(score)
              )

            if (nrow(gene_df) == 0) {
              NULL
            } else {
              score_sd <- stats::sd(
                gene_df$score,
                na.rm = TRUE
              )

              if (
                !is.finite(score_sd) ||
                score_sd == 0
              ) {
                NULL
              } else {
                gene_df |>
                  dplyr::mutate(
                    score_z = as.numeric(
                      scale(score)
                    )
                  )
              }
            }
          },
          error = function(e) NULL
        )
      }


      df_list <- shiny::withProgress(
        message = paste0(
          "Preparing TF enrichment: ",
          ct
        ),

        detail = paste0(
          "Reading ",
          length(genes_for_enrichment),
          " gene-level Feature Importance files",
          if (use_parallel) {
            paste0(
              " using ",
              pastaa_n_cores,
              " cores"
            )
          } else {
            ""
          }
        ),

        value = 0,

        {
          n_genes <- length(
            genes_for_enrichment
          )

          if (use_parallel) {

            n_chunks <- min(
              20L,
              n_genes
            )

            chunk_ids <- ceiling(
              seq_along(
                genes_for_enrichment
              ) /
                (
                  n_genes /
                    n_chunks
                )
            )

            gene_chunks <- split(
              genes_for_enrichment,
              chunk_ids
            )

            out_list <- vector(
              "list",
              0
            )

            n_done <- 0L

            for (chunk in gene_chunks) {

              chunk_results <-
                parallel::mclapply(
                  chunk,
                  read_one_gene_fi,
                  mc.cores = pastaa_n_cores
                )

              out_list <- c(
                out_list,
                chunk_results
              )

              n_done <-
                n_done +
                length(chunk)

              shiny::setProgress(
                value = n_done / n_genes,

                detail = paste0(
                  "Read ",
                  format(
                    n_done,
                    big.mark = ","
                  ),
                  " / ",
                  format(
                    n_genes,
                    big.mark = ","
                  ),
                  " genes (",
                  pastaa_n_cores,
                  " cores)"
                )
              )
            }

          } else {

            progress_every <- max(
              1L,
              floor(
                n_genes / 100L
              )
            )

            out_list <- vector(
              "list",
              n_genes
            )

            for (
              i in seq_along(
                genes_for_enrichment
              )
            ) {

              out_list[[i]] <-
                read_one_gene_fi(
                  genes_for_enrichment[[i]]
                )

              if (
                i == 1L ||
                i == n_genes ||
                i %% progress_every == 0L
              ) {

                shiny::setProgress(
                  value = i / n_genes,

                  detail = paste0(
                    "Read ",
                    format(
                      i,
                      big.mark = ","
                    ),
                    " / ",
                    format(
                      n_genes,
                      big.mark = ","
                    ),
                    " genes"
                  )
                )
              }
            }
          }

          out_list
        }
      )

      df <- dplyr::bind_rows(df_list) |>
        dplyr::filter(
          is.finite(score),
          is.finite(score_z)
        )
      message(
        "PASTAA cell-type mode: calculated within-gene SHAP z-scores for ",
        dplyr::n_distinct(df$gene_id),
        " genes | raw SHAP range: ",
        paste(signif(range(df$score, na.rm = TRUE), 4), collapse = " to "),
        " | z-score range: ",
        paste(signif(range(df$score_z, na.rm = TRUE), 4), collapse = " to ")
      )
    }
    
    if (run_mode == "celltype") {
      df <- df |>
        dplyr::mutate(score_for_ranking = score_z)
    } else {
      df <- df |>
        dplyr::mutate(score_for_ranking = score)
    }

    if (isTRUE(use_abs_shap)) {

      # Absolute: rank all regions by absolute Feature Importance
      df <- df |>
        dplyr::mutate(
          rank_score = abs(score_for_ranking)
        ) |>
        dplyr::arrange(
          dplyr::desc(rank_score)
        )

    } else if (direction == "pos") {

      # Positive: rank ALL regions from highest to lowest Feature Importance
      df <- df |>
        dplyr::mutate(
          rank_score = score_for_ranking
        ) |>
        dplyr::arrange(
          dplyr::desc(rank_score)
        )

    } else if (direction == "neg") {

      # Negative: rank ALL regions from lowest to highest Feature Importance
      df <- df |>
        dplyr::mutate(
          rank_score = score_for_ranking
        ) |>
        dplyr::arrange(
          rank_score
        )
    }
    message(
      "PASTAA ranking check | mode=", run_mode,
      " | ct=", ct,
      " | direction=", direction,
      " | first scores=",
      paste(
        signif(head(df$score_for_ranking, 5), 4),
        collapse = ", "
      )
    )

    if (nrow(df) == 0) return(NULL)

    coords <- parse_interaction_safe(df$interaction)

    df2 <- dplyr::bind_cols(df, coords) |>
      dplyr::mutate(
        chrom = as.character(chrom),
        chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
        start = suppressWarnings(as.integer(start)),
        end = suppressWarnings(as.integer(end)),
        start2 = pmin(start, end),
        end2 = pmax(start, end)
      ) |>
      dplyr::filter(
        !is.na(chrom),
        !is.na(start2),
        !is.na(end2),
        end2 > start2
      ) |>
      dplyr::mutate(
        region_id = paste0(chrom, ":", start2, "-", end2)
      ) |>
      dplyr::distinct(region_id, .keep_all = TRUE)

    requested_n <- as.integer(n_regions)

    if (is.na(requested_n) || requested_n <= 0) {
      stop("Number of top regions must be a positive integer.")
    }

    if (nrow(df2) < requested_n) {
      stop(
        "Not enough valid regions for PASTAA. Requested ",
        requested_n,
        " regions, but only ",
        nrow(df2),
        " are available for cell type ",
        ct,
        " and direction ",
        direction,
        ". Please lower the number of regions or use absolute Feature Importance values."
      )
    }

    result <- df2 |> dplyr::slice(seq_len(requested_n))

    if (run_mode == "celltype" && !is.null(cache_key)) {
      assign(cache_key, result, envir = pastaa_region_cache)
      message(
        "PASTAA cell-type mode: cached ",
        nrow(result),
        " ranked unique regions for ",
        ct
      )
    }

    result
  }
  
  run_pastaa_for_gene_ct_direction <- function(gene, ct, direction,
                                             n_regions = 500,
                                             run_mode = "gene",
                                             use_abs_shap = FALSE,
                                             filter_performance = TRUE,
                                             performance_metric = default_test_column,
                                             min_test_correlation = 0.2) {

    message("Running PASTAA: ", gene, " | ", ct, " | ", direction, " | mode=", run_mode)

    region_df <- prepare_pastaa_regions(
      gene = gene,
      ct = ct,
      direction = direction,
      n_regions = n_regions,
      run_mode = run_mode,
      use_abs_shap = use_abs_shap,
      filter_performance = filter_performance,
      performance_metric = performance_metric,
      min_test_correlation = min_test_correlation
    )

    if (is.null(region_df) || nrow(region_df) == 0) {
      stop("No regions found for ", ct, " / ", direction)
    }

    gene_safe <- if (run_mode == "celltype") {
      "ALL_GENES"
    } else {
      make_safe_id(gene)
    }

    ct_safe <- make_safe_id(ct)
    
    metric_safe <- if (isTRUE(filter_performance)) {
      make_safe_id(performance_metric)
    } else {
      "noPerfFilter"
    }

    corr_safe <- if (isTRUE(filter_performance)) {
      paste0(
        "corr",
        gsub(
          "\\.",
          "p",
          format(min_test_correlation, trim = TRUE, scientific = FALSE)
        )
      )
    } else {
      "corrNA"
    }

    if (identical(run_mode, "gene")) {

      run_id <- paste(
        gene_safe,
        direction,
        paste0("n", n_regions),
        sep = "_"
      )

    } else {

      run_id_parts <- c(
        ct_safe,
        direction,
        paste0("n", n_regions)
      )

      if (isTRUE(filter_performance)) {
        run_id_parts <- c(
          run_id_parts,
          corr_safe
        )
      }

      run_id <- paste(
        run_id_parts,
        collapse = "_"
      )

      run_id <- paste(
        run_id_parts,
        collapse = "_"
      )
    }

    bed_file <- file.path(pastaa_tmp_dir, paste0(run_id, ".bed"))
    fasta_file <- file.path(pastaa_tmp_dir, paste0(run_id, ".fa"))
    ranked_file <- file.path(pastaa_tmp_dir, paste0(run_id, "_ranked_regions.txt"))
    affinity_file <- file.path(pastaa_tmp_dir, paste0(run_id, "_Affinity.txt"))
    affinity_patched_file <- file.path(pastaa_tmp_dir, paste0(run_id, "_Affinity_patched.txt"))
    raw_out <- file.path(pastaa_out_dir, paste0(run_id, "_pastaa_raw.txt"))
    sorted_out <- file.path(pastaa_out_dir, paste0(run_id, "_pastaa_sorted.txt"))

    bed <- region_df |>
      dplyr::select(chrom, start = start2, end = end2)

    write.table(
      bed, bed_file,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )

    ranked <- region_df |>
      dplyr::mutate(rank = dplyr::row_number()) |>
      dplyr::transmute(region_id, rank)

    write.table(
      ranked, ranked_file,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )

    if (!file.exists(pastaa_hg38_fa)) {
      stop("hg38 FASTA not found: ", pastaa_hg38_fa)
    }

    getfasta_status <- system2(
      "bedtools",
      args = c(
        "getfasta",
        "-fi", pastaa_hg38_fa,
        "-bed", bed_file,
        "-fo", fasta_file
      )
    )

    if (!identical(getfasta_status, 0L) || !file.exists(fasta_file)) {
      stop("bedtools getfasta failed for ", run_id)
    }

    trap_status <- system2(
      pastaa_trap_bin,
      args = c(pastaa_energy_matrix, fasta_file),
      stdout = affinity_file
    )

    if (!identical(trap_status, 0L) || !file.exists(affinity_file)) {
      stop("TRAP failed for ", run_id)
    }

    # ===== CHANGED: restore complete motif names in TRAP affinity header =====
    #
    # The TRAP binary used here currently removes the first two characters
    # from motif labels in the affinity output (for example RUNX1 -> NX1).
    # The affinity values and motif-column order remain unchanged. Therefore,
    # restore the labels from the original energy matrix in the same order
    # before passing the affinity matrix to PASTAA.
    lines <- readLines(affinity_file, warn = FALSE)

    if (length(lines) == 0L || !nzchar(lines[1])) {
      stop("TRAP affinity output is empty or has no header: ", affinity_file)
    }

    energy_lines <- readLines(pastaa_energy_matrix, warn = FALSE)
    motif_header_lines <- energy_lines[grepl("^>", energy_lines)]

    correct_motif_names <- sub("^>", "", motif_header_lines)
    correct_motif_names <- sub("\\t.*$", "", correct_motif_names)
    correct_motif_names <- trimws(correct_motif_names)
    correct_motif_names <- correct_motif_names[
      !is.na(correct_motif_names) & nzchar(correct_motif_names)
    ]

    affinity_header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]

    if (length(affinity_header) < 2L) {
      stop("Malformed TRAP affinity header: ", lines[1])
    }

    n_affinity_motifs <- length(affinity_header) - 1L

    if (length(correct_motif_names) != n_affinity_motifs) {
      stop(
        "Motif-count mismatch between energy matrix and TRAP affinity output. ",
        "Energy matrix motifs: ", length(correct_motif_names),
        "; TRAP affinity motifs: ", n_affinity_motifs,
        ". The affinity header was not modified."
      )
    }

    # Replace only the damaged labels. The score columns remain untouched.
    restored_header <- c(affinity_header[1], correct_motif_names)

    # Keep the existing PASTAA workaround: add one dummy motif header.
    lines[1] <- paste(
      c(restored_header, "DUMMY_MOTIF"),
      collapse = "\t"
    )

    writeLines(lines, affinity_patched_file)

    message(
      "Restored ", length(correct_motif_names),
      " complete motif names in TRAP affinity header."
    )
    # ===== END CHANGED ========================================================

    pastaa_time <- system.time({
      pastaa_status <- system2(
        pastaa_bin,
        args = c(affinity_patched_file, ranked_file),
        stdout = raw_out
      )
    })

    message("PASTAA runtime for ", run_id, ": ", round(pastaa_time[["elapsed"]], 2), " sec")

    if (!identical(pastaa_status, 0L) || !file.exists(raw_out)) {
      stop("PASTAA failed for ", run_id)
    }

    res <- read.table(
      raw_out,
      sep = "\t",
      header = FALSE,
      quote = "",
      stringsAsFactors = FALSE,
      fill = TRUE
    )

    res <- res |>
      dplyr::arrange(as.numeric(V2))

    write.table(
      res, sorted_out,
      sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE
    )

    colnames(res)[seq_len(min(ncol(res), 7))] <- c(
      "TF",
      "p_value",
      "optimal_targets_in_tissue",
      "optimal_genes_in_tissue",
      "optimal_all_targets",
      "num_genes",
      "num_user_genes"
    )[seq_len(min(ncol(res), 7))]

    gene_label_for_result <- if (run_mode == "celltype") "ALL_GENES" else gene

    res <- res |>
      dplyr::filter(TF != "DUMMY_MOTIF") |>
      dplyr::mutate(
        direction = direction,
        cell_type = ct,
        gene_id = gene_label_for_result,
        run_mode = run_mode,
        BH_FDR = stats::p.adjust(as.numeric(p_value), method = "BH")
      ) |>
      dplyr::select(
        gene_id,
        cell_type,
        direction,
        run_mode,
        TF,
        p_value,
        BH_FDR,
        optimal_targets_in_tissue,
        optimal_genes_in_tissue,
        optimal_all_targets,
        num_genes,
        num_user_genes
      )

    list(
      result = res,
      sorted_file = sorted_out,
      ranked_file = ranked_file,
      fasta_file = fasta_file,
      affinity_file = affinity_patched_file,
      n_regions = nrow(region_df),
      run_id = run_id,
      run_mode = run_mode
    )
  }

  ensure_fimo_fasta <- function(
    gene,
    ct,
    direction,
    run_id,
    n_regions = 500,
    run_mode = "gene",
    use_abs_shap = TRUE
  ) {
    fasta_file <- file.path(pastaa_tmp_dir, paste0(run_id, ".fa"))

    if (file.exists(fasta_file)) {
      return(fasta_file)
    }

    region_direction <- direction

    if (isTRUE(use_abs_shap)) {
      region_direction <- "abs"
    }

    region_df <- prepare_pastaa_regions(
      gene = gene,
      ct = ct,
      direction = region_direction,
      n_regions = n_regions,
      run_mode = run_mode,
      use_abs_shap = use_abs_shap
    )

    if (is.null(region_df) || nrow(region_df) == 0) {
      stop("No regions found for FIMO FASTA: ", run_id)
    }

    bed_file <- file.path(pastaa_tmp_dir, paste0(run_id, ".bed"))

    bed <- region_df |>
      dplyr::select(chrom, start = start2, end = end2)

    write.table(
      bed,
      bed_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    if (!file.exists(pastaa_hg38_fa)) {
      stop("hg38 FASTA not found: ", pastaa_hg38_fa)
    }

    getfasta_status <- system2(
      "bedtools",
      args = c(
        "getfasta",
        "-fi", pastaa_hg38_fa,
        "-bed", bed_file,
        "-fo", fasta_file
      )
    )

    if (!identical(getfasta_status, 0L) || !file.exists(fasta_file)) {
      stop("bedtools getfasta failed for FIMO FASTA: ", run_id)
    }

    fasta_file
  }

  map_motif_ids_for_background <- function(
    motif_ids,
    motif_file,
    bg_type,
    motif_label_map = NULL
  ) {
    motif_ids <- unique(as.character(motif_ids))
    motif_ids <- trimws(motif_ids)
    motif_ids <- motif_ids[nzchar(motif_ids)]

    if (bg_type == "genome") {
      motif_ids <- motif_ids[grepl("^MA[0-9]{4}\\.[0-9]+$", motif_ids)]

      if (length(motif_ids) == 0) {
        stop("No valid JASPAR MA motif IDs found for whole-genome FIMO.")
      }

      return(motif_ids)
    }

    if (bg_type == "encode") {
      if (!file.exists(motif_file)) {
        stop("ENCODE motif file not found: ", motif_file)
      }

      lines <- readLines(motif_file, warn = FALSE)
      motif_lines <- grep("^MOTIF[[:space:]]+", lines, value = TRUE)

      if (length(motif_lines) == 0) {
        stop("No MOTIF lines found in ENCODE motif file: ", motif_file)
      }

      # Take exactly the first token after MOTIF.
      # Example:
      # "MOTIF RORA(MA0071.1) " -> "RORA(MA0071.1)"
      encode_motif_ids <- vapply(
        strsplit(trimws(motif_lines), "[[:space:]]+"),
        function(x) {
          if (length(x) >= 2) x[2] else NA_character_
        },
        character(1)
      )

      encode_motif_ids <- unique(encode_motif_ids)
      encode_motif_ids <- encode_motif_ids[!is.na(encode_motif_ids)]
      encode_motif_ids <- trimws(encode_motif_ids)
      encode_motif_ids <- encode_motif_ids[nzchar(encode_motif_ids)]

      mapped <- character(0)
      missing <- character(0)

      for (id in motif_ids) {
        id <- trimws(id)
        if (!nzchar(id)) next

        # 1) Exact match
        exact_hit <- encode_motif_ids[encode_motif_ids == id]

        if (length(exact_hit) > 0) {
          mapped <- c(mapped, exact_hit)
          next
        }

        # 2) Match exact MA ID inside ENCODE motif name.
        # Example:
        # MA0071.1 -> RORA(MA0071.1)
        ma_hit <- encode_motif_ids[
          grepl(id, encode_motif_ids, fixed = TRUE)
        ]

        if (length(ma_hit) > 0) {
          mapped <- c(mapped, ma_hit)
          next
        }

        # 3) Fallback: ignore JASPAR motif version.
        # Example:
        # MA0139.2 can match CTCF(MA0139.1)
        ma_base <- sub("\\.[0-9]+$", "", id)

        ma_base_hit <- encode_motif_ids[
          grepl(ma_base, encode_motif_ids, fixed = TRUE)
        ]

        if (length(ma_base_hit) > 0) {
          mapped <- c(mapped, ma_base_hit)
          next
        }
        # 4) Fallback: use the TF name associated with this JASPAR motif.
        # This is especially important for saved PASTAA results because
        # the ENCODE motif file may contain TF names without JASPAR IDs.
        tf_label <- NULL

        if (
          !is.null(motif_label_map) &&
          length(motif_label_map) > 0 &&
          id %in% names(motif_label_map)
        ) {
          tf_label <- as.character(motif_label_map[[id]])
          tf_label <- trimws(tf_label)
        }

        if (!is.null(tf_label) && nzchar(tf_label)) {

          # First try exact TF-name match.
          tf_hit <- encode_motif_ids[
            toupper(encode_motif_ids) == toupper(tf_label)
          ]

          # Then allow ENCODE names such as TFNAME(MAxxxx.x).
          if (length(tf_hit) == 0) {
            tf_hit <- encode_motif_ids[
              startsWith(
                toupper(encode_motif_ids),
                paste0(toupper(tf_label), "(")
              )
            ]
          }

          if (length(tf_hit) > 0) {
            mapped <- c(mapped, tf_hit)
            next
          }
        }

        missing <- c(missing, id)
      }

      mapped <- unique(mapped)

      message(
        "ENCODE motif mapping check | input: ",
        paste(motif_ids, collapse = ", "),
        " | mapped: ",
        if (length(mapped) > 0) paste(mapped, collapse = ", ") else "NONE",
        " | missing: ",
        if (length(missing) > 0) paste(missing, collapse = ", ") else "NONE",
        " | example ENCODE motifs: ",
        paste(head(encode_motif_ids, 10), collapse = ", ")
      )

      if (length(mapped) == 0) {
        stop(
          "None of the selected JASPAR motif IDs could be mapped to ENCODE motif IDs. Missing: ",
          paste(missing, collapse = ", "),
          ". Example ENCODE motifs parsed from file: ",
          paste(head(encode_motif_ids, 10), collapse = ", ")
        )
      }

      return(mapped)
    }

    stop("Unknown FIMO background type: ", bg_type)
  }

  run_fimo_for_pastaa <- function(run_id,
                                top_n = 10,
                                motif_ids_override = NULL,
                                motif_label_map = NULL,
                                bg_type = "encode",
                                bg_file = NULL) {
    bg_type <- match.arg(bg_type, c("encode", "genome"))

    fasta_file <- file.path(pastaa_tmp_dir, paste0(run_id, ".fa"))

    if (!file.exists(fasta_file)) {
      stop("FASTA file not found: ", fasta_file)
    }

    if (!file.exists(FIMO_BIN)) {
      stop("FIMO binary not found: ", FIMO_BIN)
    }

    # ----------------------------------------------------------------------
    # Choose motif file and background logic
    # ----------------------------------------------------------------------
    # ENCODE:
    #   Uses the ENCODE-specific MEME file with embedded background.
    #   FIMO background option: --bgfile motif-file
    #
    # Whole genome:
    #   Uses normal JASPAR MEME file plus external hg38 Markov background.
    #   FIMO background option: --bgfile hg38_background.txt
    # ----------------------------------------------------------------------

    if (bg_type == "encode") {
      motif_file_for_fimo <- FIMO_ENCODE_BG
      bg_arg <- "motif-file"

      if (!file.exists(motif_file_for_fimo)) {
        stop("ENCODE MEME motif/background file not found: ", motif_file_for_fimo)
      }

    } else {
      motif_file_for_fimo <- MEME_MOTIF_FILE

      if (!file.exists(motif_file_for_fimo)) {
        stop("MEME motif file not found: ", motif_file_for_fimo)
      }

      if (is.null(bg_file)) {
        bg_file <- ensure_fimo_bg_file(bg_type)
      }

      if (!file.exists(bg_file)) {
        stop("FIMO background file not found: ", bg_file)
      }

      bg_arg <- bg_file
    }

    out_dir <- file.path(fimo_out_dir, paste0(run_id, "_", bg_type))

    if (dir.exists(out_dir)) {
      unlink(out_dir, recursive = TRUE, force = TRUE)
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    # ----------------------------------------------------------------------
    # Get motif IDs
    # ----------------------------------------------------------------------

    if (!is.null(motif_ids_override)) {
      motif_ids <- unique(as.character(motif_ids_override))
      motif_ids <- motif_ids[nzchar(motif_ids)]

    } else {
      pastaa_res <- file.path(pastaa_out_dir, paste0(run_id, "_pastaa_sorted.txt"))

      if (!file.exists(pastaa_res)) {
        stop("PASTAA result not found: ", pastaa_res)
      }

      pastaa_df <- read.table(
        pastaa_res,
        header = FALSE,
        sep = "",
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        fill = TRUE
      )

      if (nrow(pastaa_df) == 0 || ncol(pastaa_df) < 2) {
        stop("PASTAA result file is empty or malformed: ", pastaa_res)
      }

      pastaa_df$V2_num <- suppressWarnings(as.numeric(pastaa_df$V2))
      pastaa_df <- pastaa_df[is.finite(pastaa_df$V2_num), , drop = FALSE]
      pastaa_df <- pastaa_df[order(pastaa_df$V2_num), , drop = FALSE]

      motif_ids <- stringr::str_extract(pastaa_df$V1, "MA[0-9]{4}\\.[0-9]+")
      motif_ids <- motif_ids[!is.na(motif_ids)]
      motif_ids <- unique(head(motif_ids, top_n))
    }

    motif_ids <- motif_ids[grepl("^MA[0-9]{4}\\.[0-9]+$", motif_ids)]

    if (length(motif_ids) == 0) {
      stop("No valid JASPAR motif IDs found for FIMO.")
    }

    # ----------------------------------------------------------------------
    # Convert motif IDs if needed
    # ----------------------------------------------------------------------
    # genome:
    #   MA0071.1
    #
    # encode:
    #   MA0071.1 -> RORA(MA0071.1), if present in ENCODE motif file
    # ----------------------------------------------------------------------

    motif_ids_for_this_bg <- map_motif_ids_for_background(
      motif_ids = motif_ids,
      motif_file = motif_file_for_fimo,
      bg_type = bg_type,
      motif_label_map = motif_label_map
    )

    writeLines(
      motif_ids_for_this_bg,
      file.path(out_dir, "top_motif_ids_exact.txt")
    )

    # ----------------------------------------------------------------------
    # Run FIMO
    # ----------------------------------------------------------------------

    motif_args <- unlist(lapply(motif_ids_for_this_bg, function(x) {
      c("--motif", shQuote(x))
    }))

    args <- c(
      "--oc", out_dir,
      "--bgfile", bg_arg,
      "--motif-pseudo", "0.1",
      "--thresh", "1e-4",
      motif_args,
      motif_file_for_fimo,
      fasta_file
    )

    message("FIMO command: ", FIMO_BIN, " ", paste(args, collapse = " "))

    fimo_run <- system2(
      FIMO_BIN,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    )

    fimo_status <- attr(fimo_run, "status")

    if (is.null(fimo_status)) {
      fimo_status <- 0L
    }

    if (!identical(fimo_status, 0L)) {
      stop(
        "FIMO failed for ",
        run_id,
        " with background ",
        bg_type,
        ". Command: ",
        FIMO_BIN,
        " ",
        paste(args, collapse = " "),
        ". FIMO output: ",
        paste(fimo_run, collapse = " | ")
      )
    }

    # ----------------------------------------------------------------------
    # Find FIMO result file
    # ----------------------------------------------------------------------

    fimo_txt <- file.path(out_dir, "fimo.tsv")

    if (!file.exists(fimo_txt)) {
      fimo_txt <- file.path(out_dir, "fimo.txt")
    }

    if (!file.exists(fimo_txt)) {
      stop("FIMO did not create fimo.tsv or fimo.txt")
    }

    # ----------------------------------------------------------------------
    # Parse FIMO result
    # ----------------------------------------------------------------------
    # Modern fimo.tsv files contain a tab-separated header and may end with
    # comment lines beginning with '#'. Reading with comment.char = '#'
    # removes the footer command and avoids coercion warnings from header or
    # comment rows.
    # ----------------------------------------------------------------------

    fimo_df <- read.delim(
      fimo_txt,
      header = TRUE,
      sep = "\t",
      comment.char = "#",
      quote = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (nrow(fimo_df) == 0) {
      stop("FIMO finished, but no motif sites were found.")
    }

    # Standardize names such as p-value -> p_value and q-value -> q_value.
    names(fimo_df) <- gsub("-", "_", names(fimo_df), fixed = TRUE)

    required_fimo_cols <- c(
      "motif_id",
      "sequence_name",
      "start",
      "stop",
      "strand",
      "score"
    )

    missing_fimo_cols <- setdiff(required_fimo_cols, names(fimo_df))

    if (length(missing_fimo_cols) > 0) {
      stop(
        "FIMO output is missing required columns: ",
        paste(missing_fimo_cols, collapse = ", "),
        ". Available columns: ",
        paste(names(fimo_df), collapse = ", "),
        ". File: ",
        fimo_txt
      )
    }

    fimo_df <- fimo_df |>
      dplyr::mutate(
        motif_id = as.character(motif_id),
        sequence_name = as.character(sequence_name),
        start = suppressWarnings(as.integer(start)),
        stop = suppressWarnings(as.integer(stop)),
        strand = as.character(strand),
        score = suppressWarnings(as.numeric(score))
      ) |>
      dplyr::filter(
        !is.na(motif_id),
        nzchar(motif_id),
        !is.na(sequence_name),
        nzchar(sequence_name),
        !is.na(start),
        !is.na(stop),
        stop >= start
      )

    if (nrow(fimo_df) == 0) {
      stop("FIMO output contained no rows with valid coordinates.")
    }

    # ----------------------------------------------------------------------
    # Convert FIMO output to BED
    # ----------------------------------------------------------------------
    # bedtools getfasta produces headers such as chr17:44900000-44900500.
    # FIMO parses these genomic headers by default and reports:
    #   sequence_name = chromosome
    #   start/stop    = genomic 1-based inclusive coordinates
    # BED requires a 0-based inclusive start and an exclusive end, hence:
    #   BED start = FIMO start - 1
    #   BED end   = FIMO stop
    # ----------------------------------------------------------------------

    motif_map <- load_meme_motif_map(motif_file_for_fimo)

    if (is.null(motif_label_map)) {
      motif_label_map <- character(0)
    }

    bed_df <- fimo_df |>
      dplyr::left_join(motif_map, by = c("motif_id" = "motif_id")) |>
      dplyr::mutate(
        chrom = sequence_name,
        hit_start = start - 1L,
        hit_end = stop,

        motif_id_clean = as.character(motif_id),

        ma_id_from_fimo = stringr::str_extract(
          motif_id_clean,
          "MA[0-9]{4}\\.[0-9]+"
        ),

        tf_from_selected_pastaa = unname(
          motif_label_map[ma_id_from_fimo]
        ),

        tf_from_encode_id = dplyr::if_else(
          grepl("\\(MA[0-9]{4}\\.[0-9]+\\)$", motif_id_clean),
          sub("\\(MA[0-9]{4}\\.[0-9]+\\)$", "", motif_id_clean),
          NA_character_
        ),

        name_raw = dplyr::case_when(
          !is.na(tf_from_selected_pastaa) &
            nzchar(tf_from_selected_pastaa) ~ tf_from_selected_pastaa,

          !is.na(tf_from_encode_id) &
            nzchar(tf_from_encode_id) ~ tf_from_encode_id,

          !is.na(tf_name) &
            nzchar(tf_name) &
            tf_name != motif_id_clean ~ as.character(tf_name),

          TRUE ~ motif_id_clean
        ),

        name = paste0(
          gsub("[^A-Za-z0-9_.:+-]+", "_", name_raw),
          "|",
          strand
        ),

        bed_score = pmin(
          1000L,
          pmax(
            0L,
            as.integer(round(score * 10))
          )
        ),

        bed_strand = strand
      ) |>
      dplyr::filter(
        !is.na(chrom),
        nzchar(chrom),
        !is.na(hit_start),
        !is.na(hit_end),
        hit_start >= 0L,
        hit_end > hit_start
      ) |>
      dplyr::select(
        chrom,
        hit_start,
        hit_end,
        name,
        bed_score,
        bed_strand
      ) |>
      dplyr::distinct() |>
      dplyr::arrange(chrom, hit_start, hit_end)

    if (nrow(bed_df) == 0) {
      stop("FIMO produced hits, but no valid BED rows remained.")
    }

    message(
      "Parsed ",
      nrow(fimo_df),
      " valid FIMO rows and created ",
      nrow(bed_df),
      " BED entries. Example: ",
      paste(
        bed_df$chrom[1],
        bed_df$hit_start[1],
        bed_df$hit_end[1],
        sep = ":"
      )
    )

    # ----------------------------------------------------------------------
    # Write BED files
    # ----------------------------------------------------------------------

    bed_file <- file.path(out_dir, paste0("fimo_sites_", bg_type, ".bed"))

    write.table(
      bed_df,
      bed_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    www_bed <- file.path(
      fimo_www_dir,
      paste0(run_id, "_", bg_type, "_fimo_sites.bed")
    )

    file.copy(bed_file, www_bed, overwrite = TRUE)

    list(
      run_id = run_id,
      bg_type = bg_type,
      bg_file = bg_arg,
      motif_file = motif_file_for_fimo,
      bed_url = paste0("fimo/", run_id, "_", bg_type, "_fimo_sites.bed"),
      n_sites = nrow(bed_df),
      motifs = motif_ids_for_this_bg,
      fimo_file = fimo_txt,
      bed_file = bed_file
    )
  }
  observeEvent(input$fimo_manual_tfs, {
    txt <- trimws(input$fimo_manual_tfs %||% "")

    if (nzchar(txt)) {
      updateRadioButtons(
        session,
        "fimo_tf_source",
        selected = "manual"
      )
    }
  }, ignoreInit = TRUE)

  observeEvent(input$show_regulatory_track, {
    gene_input <- isolate(input$igv_gene)
    selected_cts <- isolate(input$igv_celltypes)

    if (is.null(gene_input) || !nzchar(gene_input)) {
      return(NULL)
    }

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      return(NULL)
    }

    load_igv_for_gene(
      gene_input = gene_input,
      selected_cts = selected_cts,
      reason = "Updating regulatory track"
    )
  }, ignoreInit = TRUE)

  # ======================================================================
  # 4. IGV
  # ======================================================================
  parse_locus_string <- function(locus) {
    # Example: chr3:12345-67890
    m <- regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", locus)
    x <- regmatches(locus, m)[[1]]

    if (length(x) != 4) {
      return(NULL)
    }

    chrom <- x[2]
    start <- as.integer(gsub(",", "", x[3]))
    end <- as.integer(gsub(",", "", x[4]))

    if (is.na(start) || is.na(end) || end <= start) {
      return(NULL)
    }

    list(
      chrom = chrom,
      start = start,
      end = end
    )
  }

  get_gtf_col <- function(gtf, candidates) {
    hit <- intersect(candidates, colnames(gtf))
    if (length(hit) > 0) hit[1] else NULL
  }

  make_reference_track_df <- function(gene, gene_label, locus) {
    loc <- parse_locus_string(locus)

    if (is.null(loc)) {
      stop("Could not parse locus for reference track: ", locus)
    }

    transcript_models <- get_refseq_models_in_region(
      refseq_df = refseq_df,
      chrom_label = loc$chrom,
      region_start = loc$start,
      region_end = loc$end,
      max_genes = 40
    )

    if (is.null(transcript_models) || length(transcript_models) == 0) {
      return(data.frame(
        chrom = character(0),
        start = numeric(0),
        end = numeric(0),
        name = character(0),
        score = numeric(0),
        strand = character(0),
        track = character(0),
        track_type = character(0),
        feature_type = character(0),
        stringsAsFactors = FALSE
      ))
    }

    rows <- list()

    for (tx in transcript_models) {
      tx_start <- max(tx$start, loc$start)
      tx_end   <- min(tx$end, loc$end)

      if (tx_start < tx_end) {
        rows[[length(rows) + 1]] <- data.frame(
          chrom = loc$chrom,
          start = tx_start,
          end = tx_end,
          name = tx$gene_name,
          score = 0,
          strand = tx$strand,
          track = "RefSeq",
          track_type = "REFSEQ",
          feature_type = "transcript",
          stringsAsFactors = FALSE
        )
      }

      exons <- tx$exons

      if (!is.null(exons) && nrow(exons) > 0) {
        ex_start <- pmax(exons$start, loc$start)
        ex_end   <- pmin(exons$end, loc$end)
        ok <- ex_start < ex_end

        if (any(ok)) {
          rows[[length(rows) + 1]] <- data.frame(
            chrom = loc$chrom,
            start = ex_start[ok],
            end = ex_end[ok],
            name = tx$gene_name,
            score = 0,
            strand = tx$strand,
            track = "RefSeq",
            track_type = "REFSEQ",
            feature_type = "exon",
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(rows) == 0) {
      return(NULL)
    }

    dplyr::bind_rows(rows)
  }

 

  make_pdf_reference_df <- function(gene, gene_label, locus) {
    loc <- parse_locus_string(locus)

    if (is.null(loc)) {
      stop("Could not parse locus for reference track: ", locus)
    }

    ref_models <- get_refseq_models_in_region(
      refseq_df = refseq_df,
      chrom_label = loc$chrom,
      region_start = loc$start,
      region_end = loc$end,
      max_genes = 40
    )

    if (is.null(ref_models) || length(ref_models) == 0) {
      return(data.frame(
        track = "RefSeq",
        track_type = "REF_EMPTY",
        xmin = loc$start,
        xmax = loc$end,
        ymin = 0.48,
        ymax = 0.52,
        color = "#2C2C8A",
        label = gene_label,
        stringsAsFactors = FALSE
      ))
    }

    # optional: keep NFATC2 + neighboring genes, like browser screenshot
    ref_models <- ref_models[
      order(
        vapply(ref_models, function(x) x$start, numeric(1)),
        vapply(ref_models, function(x) x$end, numeric(1))
      )
    ]

    n <- length(ref_models)

    # y rows from top to bottom
    y_rows <- seq(0.82, 0.18, length.out = n)

    rows <- list()

    for (i in seq_along(ref_models)) {
      tx <- ref_models[[i]]
      y <- y_rows[i]

      tx_start <- max(tx$start, loc$start)
      tx_end   <- min(tx$end, loc$end)

      if (tx_start < tx_end) {
        rows[[length(rows) + 1]] <- data.frame(
          track = "RefSeq",
          track_type = "REF_LINE",
          xmin = tx_start,
          xmax = tx_end,
          ymin = y,
          ymax = y,
          color = "#2C2C8A",
          label = tx$gene_name,
          stringsAsFactors = FALSE
        )
      }

      exons <- tx$exons

      if (!is.null(exons) && nrow(exons) > 0) {
        exon_half_height <- 0.025

        exon_df <- data.frame(
          track = "RefSeq",
          track_type = "REF_EXON",
          xmin = pmax(exons$start, loc$start),
          xmax = pmin(exons$end, loc$end),
          ymin = y - exon_half_height,
          ymax = y + exon_half_height,
          color = "#2C2C8A",
          label = tx$gene_name,
          stringsAsFactors = FALSE
        ) |>
          dplyr::filter(xmax > xmin) |>
          dplyr::distinct()

        if (nrow(exon_df) > 0) {
          rows[[length(rows) + 1]] <- exon_df
        }
      }
    }

    if (length(rows) == 0) {
      return(NULL)
    }

    dplyr::bind_rows(rows)
  }

  make_vector_track_df <- function(gene, gene_label, selected_cts, locus) {
    loc <- parse_locus_string(locus)

    if (is.null(loc)) {
      stop("Could not parse IGV locus: ", locus)
    }

    out <- list()

    # Reference / gene model track on top
    out[[length(out) + 1]] <- make_pdf_reference_df(
      gene = gene,
      gene_label = gene_label,
      locus = locus
    )

    for (ct in selected_cts) {
      col <- get_ct_color(ct, fallback = "#666666")

      # SHAP
      shap_df <- tryCatch(
        load_igv_shap_df(gene, ct),
        error = function(e) NULL
      )

      if (!is.null(shap_df) && nrow(shap_df) > 0) {
        coords <- parse_interaction_safe(shap_df$interaction)

        shap_plot_df <- dplyr::bind_cols(shap_df, coords) |>
          dplyr::mutate(
            chrom = as.character(chrom),
            chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
            start = suppressWarnings(as.integer(start)),
            end = suppressWarnings(as.integer(end)),
            start2 = pmin(start, end),
            end2 = pmax(start, end),
            score = suppressWarnings(as.numeric(score)),
            track = paste0(
              "Feature Importance (",
              input$feature_importance_type,
              ") - ",
              ct
            ),
            track_type = "SHAP",
            color = col,
            label = NA_character_
          ) |>
          dplyr::filter(
            chrom == loc$chrom,
            !is.na(start2),
            !is.na(end2),
            is.finite(score),
            end2 > loc$start,
            start2 < loc$end
          ) |>
          dplyr::mutate(
            xmin = pmax(start2, loc$start),
            xmax = pmin(end2, loc$end),
            ymin = pmin(0, score),
            ymax = pmax(0, score)
          ) |>
          dplyr::select(track, track_type, xmin, xmax, ymin, ymax, color, label)

        out[[length(out) + 1]] <- shap_plot_df
      }

      # ATAC
      atac_df <- tryCatch(
        load_atac_df(gene, ct),
        error = function(e) NULL
      )

      if (!is.null(atac_df) && nrow(atac_df) > 0) {
        coords <- parse_interaction_safe(atac_df$interaction)

        atac_plot_df <- dplyr::bind_cols(atac_df, coords) |>
          dplyr::mutate(
            chrom = as.character(chrom),
            chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
            start = suppressWarnings(as.integer(start)),
            end = suppressWarnings(as.integer(end)),
            start2 = pmin(start, end),
            end2 = pmax(start, end),
            score = suppressWarnings(as.numeric(score)),
            track = paste0(
              "Epigenetic Signal (",
              input$epigenetic_signal_type,
              ") - ",
              ct
            ),
            track_type = "ATAC",
            color = mix_hex(col, mix_with = "#BDBDBD", amount = 0.60),
            label = NA_character_
          ) |>
          dplyr::filter(
            chrom == loc$chrom,
            !is.na(start2),
            !is.na(end2),
            is.finite(score),
            end2 > loc$start,
            start2 < loc$end
          ) |>
          dplyr::mutate(
            xmin = pmax(start2, loc$start),
            xmax = pmin(end2, loc$end),
            ymin = 0,
            ymax = score
          ) |>
          dplyr::select(track, track_type, xmin, xmax, ymin, ymax, color, label)

        out[[length(out) + 1]] <- atac_plot_df
      }
    }

    # FIMO tracks at the bottom, if FIMO was run in the current session
    fimo_df <- current_fimo_export_df()

    if (!is.null(fimo_df) && nrow(fimo_df) > 0) {
      fimo_plot_df <- fimo_df |>
        dplyr::mutate(
          chrom = as.character(chrom),
          start = suppressWarnings(as.integer(start)),
          end = suppressWarnings(as.integer(end)),
          fimo_label = ifelse(
            !is.na(fimo_label) & nzchar(fimo_label),
            fimo_label,
            gsub("\\|", " ", as.character(name))
          )
        ) |>
        dplyr::filter(
          chrom == loc$chrom,
          !is.na(start),
          !is.na(end),
          end > loc$start,
          start < loc$end
        ) |>
        dplyr::mutate(
          xmin = pmax(start, loc$start),
          xmax = pmin(end, loc$end),
          ymin = 0.43,
          ymax = 0.57,
          color = "#7B1FA2",
          label = fimo_label
        ) |>
        dplyr::filter(xmax > xmin) |>
        dplyr::select(track, track_type, xmin, xmax, ymin, ymax, color, label)

      if (nrow(fimo_plot_df) > 0) {
        out[[length(out) + 1]] <- fimo_plot_df
      }
    }

    df <- dplyr::bind_rows(out)

    if (nrow(df) == 0) {
      stop("No data available for the selected locus.")
    }

    signal_levels <- unlist(lapply(selected_cts, function(ct) {
      c(
        paste0(
          "Feature Importance (",
          input$feature_importance_type,
          ") - ",
          ct
        ),
        paste0(
          "Epigenetic Signal (",
          input$epigenetic_signal_type,
          ") - ",
          ct
        )
      )
    }))

    fimo_levels <- df |>
      dplyr::filter(track_type == "FIMO") |>
      dplyr::pull(track) |>
      unique()

    track_levels <- c("RefSeq", signal_levels, fimo_levels)
    track_levels <- track_levels[track_levels %in% unique(as.character(df$track))]

    df$track <- factor(as.character(df$track), levels = track_levels)

    df
  }

  output$igv_vector_pdf <- downloadHandler(
    filename = function() {
      gene_resolved <- resolve_igv_gene_input(input$igv_gene)

      gene_label <- if (!is.null(gene_resolved)) {
        gene_resolved$gene_display
      } else {
        "selected_gene"
      }

      paste0(
        "IGV_tracks_",
        make_safe_id(gene_label),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      gene_resolved <- resolve_igv_gene_input(input$igv_gene)

      shiny::validate(
        shiny::need(!is.null(gene_resolved), "Selected gene was not found.")
      )

      gene <- gene_resolved$gene_id
      gene_label <- gene_resolved$gene_display

      selected_cts <- input$igv_celltypes

      if (is.null(selected_cts) || length(selected_cts) == 0) {
        selected_cts <- if (length(igv_cell_types) > 0) sort(igv_cell_types)[1] else character(0)
      }

      selected_cts <- intersect(selected_cts, igv_cell_types)

      shiny::validate(
        shiny::need(length(selected_cts) > 0, "No valid cell types selected.")
      )

      ensure_gene_tracks(gene, selected_cts)

      locus <- calc_gene_locus_from_shap(gene, pad = 50000)

      locus <- calc_gene_locus_from_shap(gene, pad = 50000)

      if (is.null(locus) || !nzchar(locus)) {
        stop("Could not determine locus for gene: ", gene)
      }

      plot_df <- make_vector_track_df(
        gene = gene,
        gene_label = gene_label,
        selected_cts = selected_cts,
        locus = locus
      )

      loc <- parse_locus_string(locus)

      if (is.null(loc)) {
        stop("Could not parse locus: ", locus)
      }

      ref_line_df <- plot_df |> dplyr::filter(track_type == "REF_LINE")
      ref_exon_df <- plot_df |> dplyr::filter(track_type == "REF_EXON")
      ref_empty_df <- plot_df |> dplyr::filter(track_type == "REF_EMPTY")
      signal_df <- plot_df |> dplyr::filter(track_type %in% c("SHAP", "ATAC"))
      fimo_df <- plot_df |> dplyr::filter(track_type == "FIMO")

      ref_label_df <- ref_line_df |>
        dplyr::mutate(
          x_label = xmin,
          y_label = pmin(ymin + 0.05, 0.86),
          label_text = as.character(label)
        ) |>
        dplyr::filter(
          !is.na(x_label),
          !is.na(y_label),
          !is.na(label_text),
          nzchar(label_text)
        ) |>
        dplyr::group_by(label_text) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup()

      if (nrow(ref_label_df) == 0) {
        ref_label_df <- data.frame(
          x_label = numeric(0),
          y_label = numeric(0),
          label_text = character(0)
        )
      }

      fimo_label_df <- fimo_df |>
        dplyr::mutate(
          x_label = (xmin + xmax) / 2,
          y_label = 0.78,
          label_text = as.character(label)
        ) |>
        dplyr::filter(
          !is.na(x_label),
          !is.na(label_text),
          nzchar(label_text)
        )

      if (nrow(fimo_label_df) > 120) {
        fimo_label_df <- fimo_label_df |>
          dplyr::group_by(track) |>
          dplyr::slice_head(n = 40) |>
          dplyr::ungroup()
      }

      p <- ggplot()

      # SHAP / ATAC
      if (nrow(signal_df) > 0) {
        p <- p +
          geom_rect(
            data = signal_df,
            aes(
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax,
              fill = color
            ),
            colour = NA,
            alpha = 0.75
          )
      }

      # Reference line
      if (nrow(ref_line_df) > 0) {
        p <- p +
          geom_segment(
            data = ref_line_df,
            aes(
              x = xmin,
              xend = xmax,
              y = ymin,
              yend = ymax
            ),
            colour = "#2C2C8A",
            size = 0.28
          )
      }

      # Reference exons
      if (nrow(ref_exon_df) > 0) {
        p <- p +
          geom_rect(
            data = ref_exon_df,
            aes(
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax
            ),
            fill = "#2C2C8A",
            colour = "#2C2C8A",
            size = 0.15
          )
      }

      # Reference fallback
      if (nrow(ref_empty_df) > 0 && nrow(ref_line_df) == 0 && nrow(ref_exon_df) == 0) {
        p <- p +
          geom_segment(
            data = ref_empty_df,
            aes(
              x = xmin,
              xend = xmax,
              y = 0.5,
              yend = 0.5
            ),
            colour = "#2C2C8A",
            size = 0.28,
            linetype = "dashed"
          )
      }

      # Reference label
      if (nrow(ref_label_df) > 0) {
        p <- p +
          geom_text(
            data = ref_label_df,
            aes(
              x = x_label,
              y = y_label,
              label = label_text
            ),
            hjust = 0,
            size = 1.5,
            colour = "#2C2C8A",
            check_overlap = TRUE
          )
      }

      # FIMO short bars
      if (nrow(fimo_df) > 0) {
        p <- p +
          geom_rect(
            data = fimo_df,
            aes(
              xmin = xmin,
              xmax = xmax,
              ymin = ymin,
              ymax = ymax
            ),
            fill = "#7B1FA2",
            colour = "#7B1FA2",
            size = 0.12,
            alpha = 0.9
          )
      }

      # FIMO labels
      if (nrow(fimo_label_df) > 0) {
        p <- p +
          geom_text(
            data = fimo_label_df,
            aes(
              x = x_label,
              y = y_label,
              label = label_text
            ),
            angle = 90,
            hjust = 0,
            vjust = 0.5,
            size = 1.2,
            colour = "#7B1FA2"
          )
      }

      p <- p +
        scale_fill_identity() +
        scale_x_continuous(
          labels = function(x) format(round(x), big.mark = ",", scientific = FALSE),
          sec.axis = dup_axis(
            name = NULL,
            labels = function(x) format(round(x), big.mark = ",", scientific = FALSE)
          )
        ) +
        facet_grid(
          track ~ .,
          scales = "free_y",
          switch = "y"
        ) +
        coord_cartesian(
          xlim = c(loc$start, loc$end),
          expand = FALSE
        ) +
        labs(
          title = NULL,
          subtitle = NULL,
          x = paste0(loc$chrom, " position"),
          y = NULL
        ) +
        theme_bw(base_size = 8) +
        theme(
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, size = 6),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          panel.spacing.y = grid::unit(0.10, "lines"),

          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),

          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),

          axis.text.x.top = element_text(size = 6),
          axis.ticks.x.top = element_line(),
          axis.text.x.bottom = element_text(size = 6),
          axis.title.x = element_text(size = 7),

          plot.title = element_blank(),
          plot.subtitle = element_blank()
        )

      n_tracks <- length(unique(plot_df$track))

      pdf_height <- max(5, n_tracks * 0.50)
      pdf_width <- 11

      grDevices::cairo_pdf(
        filename = file,
        width = pdf_width,
        height = pdf_height,
        onefile = TRUE
      )

      print(p)

      grDevices::dev.off()
    }
  )

  okabe_ito <- c(
    orange        = "#E69F00",
    skyblue       = "#56B4E9",
    bluishgreen   = "#009E73",
    yellow        = "#F0E442",
    blue          = "#0072B2",
    vermillion    = "#D55E00",
    reddishpurple = "#CC79A7",
    black         = "#000000"
  )

  

  lighten <- function(hex, amount = 0.35) mix_hex(hex, "#FFFFFF", amount)
  darken  <- function(hex, amount = 0.20) mix_hex(hex, "#000000", amount)

  # Automatically assign one qualitative, colorblind-friendly color per
  # existing base cell type. Conditions such as "Astrocytes:Alzheimer's"
  # and "Astrocytes:Unaffected" intentionally share the Astrocytes color.
  base_cell_types <- unique(sub(":.*$", "", igv_cell_types))
  base_cell_types <- sort(base_cell_types[!is.na(base_cell_types) & nzchar(base_cell_types)])

  ct_palette <- stats::setNames(
    grDevices::hcl.colors(
      n = length(base_cell_types),
      palette = "Dark 3"
    ),
    base_cell_types
  )

  get_ct_color <- function(ct, fallback = "#666666") {
    base_ct <- sub(":.*$", "", as.character(ct)[1])
    col <- unname(ct_palette[base_ct])

    if (length(col) == 0 || is.na(col) || !nzchar(col)) {
      return(fallback)
    }

    col
  }

  status_msg <- reactiveVal("")
  output$status <- renderText(status_msg())
  observeEvent(input$igv_handlers_ready, {
    status_msg("IGV JavaScript handlers ready.")
  })


  observeEvent(input$igv_js_loaded, {
    status_msg(
      paste0(
        "IGV JavaScript loaded | typeof igv: ",
        input$igv_js_loaded$igv_type
      )
    )
  })

  observeEvent(input$igv_config_received, {
    status_msg(
      paste0(
        "IGV configuration received by JavaScript",
        "\nGenome: ", input$igv_config_received$genome,
        "\nLocus: ", input$igv_config_received$locus,
        "\nTracks: ", input$igv_config_received$n_tracks
      )
    )
  })
  observeEvent(input$igv_track_add_error, {
    err <- input$igv_track_add_error

    status_msg(
      paste0(
        "Could not add track to IGV",
        if (!is.null(err$track)) paste0("\nTrack: ", err$track) else "",
        if (!is.null(err$error)) paste0("\nError: ", err$error) else ""
      )
    )
  })

  observeEvent(input$igv_browser_created, {
    igv_loaded(TRUE)

    status_msg(
      paste0(
        "IGV browser created successfully",
        "\nLocus: ", input$igv_browser_created$locus,
        "\nTracks: ", input$igv_browser_created$n_tracks
      )
    )
  })

  observeEvent(input$igv_track_load_error, {
    igv_loaded(FALSE)

    status_msg(
      paste0(
        "IGV track failed",
        "\nTrack: ", input$igv_track_load_error$track,
        "\nURL: ", input$igv_track_load_error$url,
        "\nError: ", input$igv_track_load_error$error
      )
    )
  })

  observeEvent(input$igv_js_error, {
    igv_loaded(FALSE)

    status_msg(
      paste0(
        "IGV JavaScript failed",
        "\n",
        input$igv_js_error
      )
    )
  })

  repro_code <- reactive({
    pastaa_gene_value <- input$pastaa_gene

    if (is.null(pastaa_gene_value) || !nzchar(pastaa_gene_value)) {
      pastaa_gene_value <- input$igv_gene
    }

    if (is.null(pastaa_gene_value) || !nzchar(pastaa_gene_value)) {
      pastaa_gene_value <- "selected_gene"
    }

    gene_resolved <- tryCatch(
      resolve_igv_gene_input(pastaa_gene_value),
      error = function(e) NULL
    )

    gene_id <- if (
      !is.null(gene_resolved) &&
      !is.null(gene_resolved$gene_id)
    ) {
      gene_resolved$gene_id
    } else {
      pastaa_gene_value
    }

    gene_label <- if (
      !is.null(gene_resolved) &&
      !is.null(gene_resolved$gene_display)
    ) {
      gene_resolved$gene_display
    } else {
      pastaa_gene_value
    }
    selected_cts <- input$igv_celltypes

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      selected_cts <- character(0)
    }

    celltypes_code <- if (length(selected_cts) > 0) {
      paste0('"', selected_cts, '"', collapse = ", ")
    } else {
      ""
    }

    fimo_tf_source_value <- input$fimo_tf_source
    if (is.null(fimo_tf_source_value) || !nzchar(fimo_tf_source_value)) {
      fimo_tf_source_value <- "saved"
    }

    fimo_manual_tfs_value <- input$fimo_manual_tfs
    if (is.null(fimo_manual_tfs_value)) {
      fimo_manual_tfs_value <- ""
    }
    fimo_manual_tfs_value <- gsub('"', '\\"', fimo_manual_tfs_value)

    fimo_top_n_value <- input$fimo_top_n
    if (is.null(fimo_top_n_value) || !is.finite(fimo_top_n_value)) {
      fimo_top_n_value <- 10
    }

    

    pastaa_n_regions_value <- input$pastaa_n_regions
    if (is.null(pastaa_n_regions_value) || !is.finite(pastaa_n_regions_value)) {
      pastaa_n_regions_value <- 500
    }

    pastaa_direction_value <- input$pastaa_direction
    if (is.null(pastaa_direction_value) || !nzchar(pastaa_direction_value)) {
      pastaa_direction_value <- "both"
    }
    pastaa_run_mode_value <- input$pastaa_run_mode
    if (is.null(pastaa_run_mode_value) || !nzchar(pastaa_run_mode_value)) {
      pastaa_run_mode_value <- "gene"
    }

    pastaa_celltypes_value <- input$pastaa_celltypes
    if (is.null(pastaa_celltypes_value)) {
      pastaa_celltypes_value <- character(0)
    }

    pastaa_celltypes_code <- if (length(pastaa_celltypes_value) > 0) {
      paste0('"', pastaa_celltypes_value, '"', collapse = ", ")
    } else {
      ""
    }

    pastaa_filter_performance_value <- isTRUE(input$pastaa_filter_performance)

    pastaa_performance_metric_value <-
      input$pastaa_performance_metric %||% default_test_column

    pastaa_min_correlation_value <- suppressWarnings(
      as.numeric(input$pastaa_min_test_correlation)[1]
    )

    if (!is.finite(pastaa_min_correlation_value)) {
      pastaa_min_correlation_value <- 0.2
    }

    pastaa_q_cutoff_value <- suppressWarnings(
      as.numeric(input$pastaa_q_cutoff)[1]
    )

    if (!is.finite(pastaa_q_cutoff_value)) {
      pastaa_q_cutoff_value <- 0.05
    }

    pastaa_save_top_n_value <- suppressWarnings(
      as.integer(input$pastaa_save_top_n)[1]
    )

    if (is.na(pastaa_save_top_n_value)) {
      pastaa_save_top_n_value <- 10L
    }

    

    paste0(
    '# EpIVoi reproducibility script
    # Generated from the current EpIVoi settings.

    library(dplyr)
    library(tibble)
    library(ggplot2)

    # =====================================================================
    # 1. INPUT OBJECT
    # =====================================================================

    viz_obj_path <- "', viz_obj_path, '"
    viz_obj <- readRDS(viz_obj_path)

    # EpIVoi object helper functions
    source("viz_obj_helpers.R")

    species <- "', species_name, '"
    genome_version <- "', genome_version, '"

    # =====================================================================
    # 2. SELECTED SETTINGS
    # =====================================================================

    gene_id <- "', gene_id, '"
    gene_symbol <- "', gene_label, '"

    igv_cell_types <- c(', celltypes_code, ')

    pastaa_run_mode <- "', pastaa_run_mode_value, '"
    pastaa_cell_types <- c(', pastaa_celltypes_code, ')
    pastaa_direction <- "', pastaa_direction_value, '"
    pastaa_n_regions <- ', pastaa_n_regions_value, '

    pastaa_filter_performance <- ', 
    if (pastaa_filter_performance_value) "TRUE" else "FALSE", '
    pastaa_performance_metric <- "', pastaa_performance_metric_value, '"
    pastaa_min_test_correlation <- ', pastaa_min_correlation_value, '

    pastaa_q_cutoff <- ', pastaa_q_cutoff_value, '
    pastaa_save_top_n <- ', pastaa_save_top_n_value, '

    fimo_tf_source <- "', fimo_tf_source_value, '"
    fimo_manual_tfs <- "', fimo_manual_tfs_value, '"
    fimo_top_n <- ', fimo_top_n_value, '

    # =====================================================================
    # 3. LOAD RAW FEATURE IMPORTANCE
    # =====================================================================

    strip_ensembl_version <- function(x) {
      sub("\\\\.[0-9]+$", "", as.character(x))
    }

    extract_region_matrix <- function(mat) {
      mat <- as.data.frame(mat, check.names = FALSE)

      if ("region" %in% colnames(mat)) {
        regions <- as.character(mat$region)
        mat$region <- NULL
      } else if ("id" %in% colnames(mat)) {
        regions <- as.character(mat$id)
        mat$id <- NULL
      } else {
        stop("Could not find region column in Feature Importance matrix.")
      }

      list(
        regions = regions,
        matrix = mat
      )
    }

    load_raw_shap <- function(gene, cell_type) {
      gene <- strip_ensembl_version(gene)

      raw <- read_feature_importance_file(viz_obj, gene)
      x <- extract_region_matrix(raw)

      if (!cell_type %in% colnames(x$matrix)) {
        stop("Cell type not found: ", cell_type)
      }

      data.frame(
        gene_id = gene,
        interaction = x$regions,
        score = suppressWarnings(
          as.numeric(x$matrix[[cell_type]])
        ),
        stringsAsFactors = FALSE
      ) |>
        dplyr::filter(
          !is.na(interaction),
          is.finite(score)
        )
    }

    # =====================================================================
    # 4. PREPARE FEATURE IMPORTANCE FOR PASTAA
    # =====================================================================
    #
    # The EpIVoi object contains RAW Feature Importance values.
    #
    # Gene-level PASTAA:
    #   raw Feature Importance values are ranked directly.
    #
    # Cell-type-level PASTAA:
    #   Feature Importance values are converted to within-gene z-scores
    #   before regions from different genes are compared.
    #
    #   z = (raw_score - mean(raw_score for gene)) /
    #       sd(raw_score for gene)
    #
    # =====================================================================

    prepare_gene_shap <- function(gene, cell_type) {

      df <- load_raw_shap(
        gene = gene,
        cell_type = cell_type
      )

      df |>
        dplyr::mutate(
          score_for_ranking = score
        )
    }

    prepare_celltype_shap <- function(
      genes,
      cell_type
    ) {

      df_list <- lapply(genes, function(g) {

        gene_df <- tryCatch(
          load_raw_shap(
            gene = g,
            cell_type = cell_type
          ),
          error = function(e) NULL
        )

        if (is.null(gene_df) || nrow(gene_df) == 0) {
          return(NULL)
        }

        score_sd <- stats::sd(
          gene_df$score,
          na.rm = TRUE
        )

        if (
          !is.finite(score_sd) ||
          score_sd == 0
        ) {
          return(NULL)
        }

        gene_df |>
          dplyr::mutate(
            score_z = as.numeric(
              scale(score)
            )
          )
      })

      dplyr::bind_rows(df_list) |>
        dplyr::filter(
          is.finite(score),
          is.finite(score_z)
        ) |>
        dplyr::mutate(
          score_for_ranking = score_z
        )
    }

    # =====================================================================
    # 5. RANK REGIONS
    # =====================================================================

    rank_regions <- function(
      df,
      direction = c(
        "pos",
        "neg",
        "abs"
      ),
      n_regions = 500
    ) {

      direction <- match.arg(direction)

      if (direction == "abs") {

        df <- df |>
          dplyr::mutate(
            rank_score = abs(score_for_ranking)
          ) |>
          dplyr::arrange(
            dplyr::desc(rank_score)
          )

      } else if (direction == "pos") {

        # Rank ALL regions from highest to lowest.
        df <- df |>
          dplyr::mutate(
            rank_score = score_for_ranking
          ) |>
          dplyr::arrange(
            dplyr::desc(rank_score)
          )

      } else {

        # Rank ALL regions from lowest to highest.
        df <- df |>
          dplyr::mutate(
            rank_score = score_for_ranking
          ) |>
          dplyr::arrange(
            rank_score
          )
      }

      df |>
        dplyr::slice_head(
          n = n_regions
        )
    }

    # =====================================================================
    # Current analysis settings
    # =====================================================================

    message(
      "EpIVoi reproducibility script loaded."
    )

    message(
      "PASTAA mode: ",
      pastaa_run_mode
    )

    message(
      "Feature Importance direction: ",
      pastaa_direction
    )

    message(
      "Number of regions: ",
      pastaa_n_regions
    )

    # =====================================================================
    # 6. SELECT GENES FOR CELL-TYPE ANALYSIS
    # =====================================================================

    get_genes_for_celltype_analysis <- function() {

      genes <- strip_ensembl_version(
        get_genes(viz_obj)
      )

      if (!isTRUE(pastaa_filter_performance)) {
        return(genes)
      }

      corr_df <- viz_obj$performance$correlation |>
        tibble::rownames_to_column("gene")

      if (!pastaa_performance_metric %in% colnames(corr_df)) {
        stop(
          "Performance metric not available: ",
          pastaa_performance_metric
        )
      }

      eligible_names <- corr_df |>
        dplyr::transmute(
          gene = strip_ensembl_version(gene),
          correlation = suppressWarnings(
            as.numeric(
              .data[[pastaa_performance_metric]]
            )
          )
        ) |>
        dplyr::filter(
          is.finite(correlation),
          correlation >= pastaa_min_test_correlation
        ) |>
        dplyr::pull(gene)

      # Performance row names may contain either Ensembl IDs or gene symbols.
      # The app resolves these against the GTF before loading Feature Importance.
      gtf <- as.data.frame(viz_obj$gtf)

      gene_id_col <- intersect(
        c("gene_id", "geneid", "ensembl_gene_id"),
        colnames(gtf)
      )[1]

      gene_name_col <- intersect(
        c("gene_name", "gene_symbol", "external_gene_name"),
        colnames(gtf)
      )[1]

      eligible_ids <- eligible_names[
        eligible_names %in% genes
      ]

      if (
        !is.na(gene_id_col) &&
        !is.na(gene_name_col)
      ) {

        mapping <- data.frame(
          gene_id = strip_ensembl_version(
            gtf[[gene_id_col]]
          ),
          gene_symbol = as.character(
            gtf[[gene_name_col]]
          ),
          stringsAsFactors = FALSE
        ) |>
          dplyr::filter(
            !is.na(gene_id),
            !is.na(gene_symbol)
          ) |>
          dplyr::distinct(
            gene_id,
            .keep_all = TRUE
          )

        symbol_ids <- mapping$gene_id[
          toupper(mapping$gene_symbol) %in%
            toupper(eligible_names)
        ]

        eligible_ids <- unique(
          c(
            eligible_ids,
            symbol_ids
          )
        )
      }

      eligible_ids[
        eligible_ids %in% genes
      ]
    }

    # =====================================================================
    # 7. PARSE GENOMIC REGIONS
    # =====================================================================

    parse_region <- function(x) {

      x <- gsub("[:_-]", " ", x)
      parts <- strsplit(x, "\\\\s+")

      chrom <- vapply(
        parts,
        function(z) z[1],
        character(1)
      )

      start <- suppressWarnings(
        as.integer(
          vapply(parts, function(z) z[2], character(1))
        )
      )

      end <- suppressWarnings(
        as.integer(
          vapply(
            parts,
            function(z) {
              if (length(z) >= 3) z[3] else NA_character_
            },
            character(1)
          )
        )
      )

      end[
        is.na(end) &
          !is.na(start)
      ] <- start[
        is.na(end) &
          !is.na(start)
      ] + 100L

      data.frame(
        chrom = chrom,
        start = start,
        end = end,
        stringsAsFactors = FALSE
      )
    }

    # =====================================================================
    # 8. RANK, PARSE AND DEDUPLICATE REGIONS
    # =====================================================================

    prepare_ranked_regions <- function(
      df,
      direction,
      n_regions
    ) {

      df <- rank_regions(
        df = df,
        direction = direction,
        n_regions = nrow(df)
      )

      coords <- parse_region(
        df$interaction
      )

      df <- dplyr::bind_cols(
        df,
        coords
      ) |>
        dplyr::mutate(
          chrom = ifelse(
            grepl("^chr", chrom),
            chrom,
            paste0("chr", chrom)
          ),
          start2 = pmin(start, end),
          end2 = pmax(start, end),
          region_id = paste0(
            chrom,
            ":",
            start2,
            "-",
            end2
          )
        ) |>
        dplyr::filter(
          !is.na(chrom),
          !is.na(start2),
          !is.na(end2),
          end2 > start2
        ) |>
        dplyr::distinct(
          region_id,
          .keep_all = TRUE
        )

      if (nrow(df) < n_regions) {
        stop(
          "Only ",
          nrow(df),
          " valid unique regions available; requested ",
          n_regions,
          "."
        )
      }

      df |>
        dplyr::slice_head(
          n = n_regions
        )
    }

    # =====================================================================
    # 9. BUILD THE CURRENT REGION SET
    # =====================================================================

    genes_for_analysis <- get_genes_for_celltype_analysis()

    message(
      "Eligible genes for cell-type analysis: ",
      length(genes_for_analysis)
    )

    # For gene mode:
    #
    # shap_df <- prepare_gene_shap(
    #   gene = gene_id,
    #   cell_type = pastaa_cell_types[1]
    # )
    #
    # For cell-type mode:
    #
    # shap_df <- prepare_celltype_shap(
    #   genes = genes_for_analysis,
    #   cell_type = pastaa_cell_types[1]
    # )
    #
    # ranked_regions <- prepare_ranked_regions(
    #   df = shap_df,
    #   direction = "pos",
    #   n_regions = pastaa_n_regions
    # )
    # =====================================================================
    # 10. WRITE BED + RANK FILES
    # =====================================================================

    write_pastaa_inputs <- function(
      ranked_regions,
      prefix = "epivoi_repro"
    ) {

      bed_file <- paste0(prefix, ".bed")
      ranked_file <- paste0(prefix, "_ranked_regions.txt")

      bed <- ranked_regions |>
        dplyr::select(
          chrom,
          start = start2,
          end = end2
        )

      write.table(
        bed,
        bed_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
      )

      ranked <- ranked_regions |>
        dplyr::mutate(
          rank = dplyr::row_number()
        ) |>
        dplyr::select(
          region_id,
          rank
        )

      write.table(
        ranked,
        ranked_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
      )

      list(
        bed_file = bed_file,
        ranked_file = ranked_file
      )
    }

    # =====================================================================
    # 11. EXTRACT FASTA WITH BEDTOOLS
    # =====================================================================

    extract_fasta <- function(
      bed_file,
      genome_fasta,
      prefix = "epivoi_repro"
    ) {

      fasta_file <- paste0(prefix, ".fa")

      status <- system2(
        "bedtools",
        args = c(
          "getfasta",
          "-fi", genome_fasta,
          "-bed", bed_file,
          "-fo", fasta_file
        )
      )

      if (!identical(status, 0L)) {
        stop("bedtools getfasta failed.")
      }

      fasta_file
    }

    # =====================================================================
    # 12. RUN TRAP
    # =====================================================================

    run_trap <- function(
      fasta_file,
      trap_bin,
      energy_matrix,
      prefix = "epivoi_repro"
    ) {

      affinity_file <- paste0(
        prefix,
        "_Affinity.txt"
      )

      status <- system2(
        trap_bin,
        args = c(
          energy_matrix,
          fasta_file
        ),
        stdout = affinity_file
      )

      if (!identical(status, 0L)) {
        stop("TRAP failed.")
      }

      affinity_file
    }

    # =====================================================================
    # 13. RESTORE TRAP MOTIF NAMES
    # =====================================================================

    patch_trap_header <- function(
      affinity_file,
      energy_matrix,
      prefix = "epivoi_repro"
    ) {

      lines <- readLines(
        affinity_file,
        warn = FALSE
      )

      if (length(lines) == 0) {
        stop("TRAP affinity file is empty.")
      }

      energy_lines <- readLines(
        energy_matrix,
        warn = FALSE
      )

      motif_header_lines <- energy_lines[
        grepl("^>", energy_lines)
      ]

      correct_motif_names <- sub(
        "^>",
        "",
        motif_header_lines
      )

      correct_motif_names <- sub(
        "\t.*$",
        "",
        correct_motif_names
      )

      correct_motif_names <- trimws(
        correct_motif_names
      )

      affinity_header <- strsplit(
        lines[1],
        "\t",
        fixed = TRUE
      )[[1]]

      n_affinity_motifs <-
        length(affinity_header) - 1L

      if (
        length(correct_motif_names) !=
          n_affinity_motifs
      ) {
        stop(
          "Motif-count mismatch between energy matrix and TRAP output."
        )
      }

      restored_header <- c(
        affinity_header[1],
        correct_motif_names,
        "DUMMY_MOTIF"
      )

      lines[1] <- paste(
        restored_header,
        collapse = "\t"
      )

      patched_file <- paste0(
        prefix,
        "_Affinity_patched.txt"
      )

      writeLines(
        lines,
        patched_file
      )

      patched_file
    }

    # =====================================================================
    # 14. RUN PASTAA
    # =====================================================================

    run_pastaa <- function(
      affinity_file,
      ranked_file,
      pastaa_bin,
      prefix = "epivoi_repro"
    ) {

      raw_output <- paste0(
        prefix,
        "_pastaa_raw.txt"
      )

      status <- system2(
        pastaa_bin,
        args = c(
          affinity_file,
          ranked_file
        ),
        stdout = raw_output
      )

      if (!identical(status, 0L)) {
        stop("PASTAA failed.")
      }

      result <- read.table(
        raw_output,
        sep = "\t",
        header = FALSE,
        quote = "",
        stringsAsFactors = FALSE,
        fill = TRUE
      )

      result <- result |>
        dplyr::arrange(
          suppressWarnings(
            as.numeric(V2)
          )
        )

      colnames(result)[
        seq_len(min(ncol(result), 7))
      ] <- c(
        "TF",
        "p_value",
        "optimal_targets_in_tissue",
        "optimal_genes_in_tissue",
        "optimal_all_targets",
        "num_genes",
        "num_user_genes"
      )[
        seq_len(min(ncol(result), 7))
      ]

      result |>
        dplyr::filter(
          TF != "DUMMY_MOTIF"
        ) |>
        dplyr::mutate(
          p_value = suppressWarnings(
            as.numeric(p_value)
          ),
          BH_FDR = stats::p.adjust(
            p_value,
            method = "BH"
          )
        ) |>
        dplyr::arrange(
          BH_FDR,
          p_value
        )
    }

    # =====================================================================
    # 15. EXTERNAL RESOURCES USED BY EPIVOI
    # =====================================================================

    genome_fasta <- "', pastaa_hg38_fa, '"
    trap_bin <- "', pastaa_trap_bin, '"
    pastaa_bin <- "', pastaa_bin, '"
    energy_matrix <- "', pastaa_energy_matrix, '"

    # =====================================================================
    # 16. RUN THE CURRENT PASTAA ANALYSIS
    # =====================================================================

    if (length(pastaa_cell_types) == 0) {
      stop(
        "No PASTAA cell type was selected when this reproducibility script was generated."
      )
    }

    directions_to_run <- switch(
      pastaa_direction,
      "pos" = "pos",
      "neg" = "neg",
      "both" = c("pos", "neg"),
      "abs" = "abs",
      stop("Unknown PASTAA direction: ", pastaa_direction)
    )

    all_pastaa_results <- list()

    for (cell_type in pastaa_cell_types) {

      if (pastaa_run_mode == "gene") {

        shap_df <- prepare_gene_shap(
          gene = gene_id,
          cell_type = cell_type
        )

      } else {

        shap_df <- prepare_celltype_shap(
          genes = genes_for_analysis,
          cell_type = cell_type
        )
      }

      for (direction_to_run in directions_to_run) {

        ranked_regions <- prepare_ranked_regions(
          df = shap_df,
          direction = direction_to_run,
          n_regions = pastaa_n_regions
        )

        run_prefix <- paste(
          "epivoi_repro",
          gsub("[^A-Za-z0-9]+", "_", cell_type),
          direction_to_run,
          sep = "_"
        )

        files <- write_pastaa_inputs(
          ranked_regions = ranked_regions,
          prefix = run_prefix
        )

        fasta_file <- extract_fasta(
          bed_file = files$bed_file,
          genome_fasta = genome_fasta,
          prefix = run_prefix
        )

        affinity_file <- run_trap(
          fasta_file = fasta_file,
          trap_bin = trap_bin,
          energy_matrix = energy_matrix,
          prefix = run_prefix
        )

        patched_affinity <- patch_trap_header(
          affinity_file = affinity_file,
          energy_matrix = energy_matrix,
          prefix = run_prefix
        )

        pastaa_result <- run_pastaa(
          affinity_file = patched_affinity,
          ranked_file = files$ranked_file,
          pastaa_bin = pastaa_bin,
          prefix = run_prefix
        ) |>
          dplyr::mutate(
            cell_type = cell_type,
            direction = direction_to_run,
            run_mode = pastaa_run_mode,
            gene_id = ifelse(
              pastaa_run_mode == "celltype",
              "ALL_GENES",
              gene_id
            )
          )

        all_pastaa_results[[
          paste(cell_type, direction_to_run, sep = "__")
        ]] <- pastaa_result
      }
    }

    pastaa_results <- dplyr::bind_rows(
      all_pastaa_results
    )

    significant_tfs <- pastaa_results |>
      dplyr::filter(
        is.finite(BH_FDR),
        BH_FDR <= pastaa_q_cutoff
      ) |>
      dplyr::group_by(
        cell_type,
        direction
      ) |>
      dplyr::arrange(
        BH_FDR,
        p_value,
        .by_group = TRUE
      ) |>
      dplyr::slice_head(
        n = pastaa_save_top_n
      ) |>
      dplyr::ungroup()

    print(significant_tfs)
    # =====================================================================
    # 17. SAVE PASTAA RESULTS
    # =====================================================================

    write.csv(
      pastaa_results,
      file = "epivoi_pastaa_all_results.csv",
      row.names = FALSE
    )

    write.csv(
      significant_tfs,
      file = "epivoi_pastaa_significant_tfs.csv",
      row.names = FALSE
    )

    # =====================================================================
    # 18. PLOT TOP ENRICHED TFS
    # =====================================================================

    plot_df <- significant_tfs |>
      dplyr::mutate(
        minus_log10_q = -log10(
          pmax(BH_FDR, .Machine$double.xmin)
        )
      )

    if (nrow(plot_df) > 0) {

      p <- ggplot(
        plot_df,
        aes(
          x = reorder(TF, minus_log10_q),
          y = minus_log10_q,
          fill = direction
        )
      ) +
        geom_col() +
        coord_flip() +
        facet_wrap(
          ~ cell_type + direction,
          scales = "free_y"
        ) +
        labs(
          x = "TF",
          y = "-log10(q-value / BH FDR)",
          title = "PASTAA TF enrichment"
        ) +
        theme_bw(base_size = 14)

      print(p)

      ggsave(
        filename = "epivoi_pastaa_top_tfs.pdf",
        plot = p,
        width = 9,
        height = 6
      )

    } else {

      message(
        "No TFs passed the selected BH-FDR cutoff."
      )
    }

    message(
      "Reproducibility analysis finished."
    )
    '
    )
  })

  output$repro_code <- renderText({
    repro_code()
  })

  output$download_repro_code <- downloadHandler(
    filename = function() {
      gene_value <- input$igv_gene
      gene_resolved <- resolve_igv_gene_input(gene_value)

      gene_label <- if (!is.null(gene_resolved)) {
        gene_resolved$gene_display
      } else {
        "selected_gene"
      }

      paste0(
        "EpIVoi_reproducibility_",
        make_safe_id(gene_label),
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".R"
      )
    },

    content = function(file) {
      writeLines(repro_code(), con = file)
    }
  )
  
  igv_loaded <- reactiveVal(FALSE)
  igv_current_gene <- reactiveVal(NULL)
  igv_previous_gene <- reactiveVal(NULL)
  current_fimo_export_df <- reactiveVal(NULL)

  igv_gene_history <- reactiveVal(list())
  save_current_igv_state <- function() {
    old_gene <- igv_current_gene()

    if (is.null(old_gene) || !nzchar(old_gene)) {
      return(NULL)
    }

    h <- igv_gene_history()

    h[[old_gene]] <- list(
      celltypes = isolate(input$igv_celltypes),
      fimo_tf_source = isolate(input$fimo_tf_source),
      fimo_manual_tfs = isolate(input$fimo_manual_tfs),
      fimo_top_n = isolate(input$fimo_top_n),
      
      timestamp = Sys.time()
    )

    igv_gene_history(h)

    invisible(NULL)
  }

  restore_igv_state_for_gene <- function(gene) {
    h <- igv_gene_history()

    if (is.null(h[[gene]])) {
      return(NULL)
    }

    saved <- h[[gene]]

    if (!is.null(saved$celltypes)) {
      updateSelectizeInput(
        session,
        "igv_celltypes",
        choices = igv_cell_types,
        selected = intersect(saved$celltypes, igv_cell_types),
        server = FALSE
      )
    }

    if (!is.null(saved$fimo_tf_source)) {
      updateRadioButtons(
        session,
        "fimo_tf_source",
        selected = saved$fimo_tf_source
      )
    }

    if (!is.null(saved$fimo_manual_tfs)) {
      updateTextAreaInput(
        session,
        "fimo_manual_tfs",
        value = saved$fimo_manual_tfs
      )
    }

    if (!is.null(saved$fimo_top_n)) {
      updateNumericInput(
        session,
        "fimo_top_n",
        value = saved$fimo_top_n
      )
    }

    

    invisible(NULL)
  }

  load_igv_for_gene <- function(gene_input, selected_cts = NULL, reason = "Loading IGV") {
    gene_resolved <- resolve_igv_gene_input(gene_input)

    if (is.null(gene_resolved)) {
      igv_loaded(FALSE)
      status_msg(paste("Gene not found:", gene_input))
      return(FALSE)
    }

    gene <- gene_resolved$gene_id
    gene_label <- gene_resolved$gene_display
    current_fimo_export_df(NULL)

    old_gene <- igv_current_gene()

    if (!is.null(old_gene) && !identical(old_gene, gene)) {
      save_current_igv_state()
      igv_previous_gene(old_gene)
    }

    restore_igv_state_for_gene(gene)

    h <- igv_gene_history()

    if (
      (is.null(selected_cts) || length(selected_cts) == 0) &&
      !is.null(h[[gene]]) &&
      !is.null(h[[gene]]$celltypes) &&
      length(h[[gene]]$celltypes) > 0
    ) {
      selected_cts <- h[[gene]]$celltypes
    }

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      selected_cts <- if (length(igv_cell_types) > 0) {
        sort(igv_cell_types)[1]
      } else {
        character(0)
      }
    }

    selected_cts <- intersect(selected_cts, igv_cell_types)

    updateSelectizeInput(
      session,
      "igv_celltypes",
      choices = igv_cell_types,
      selected = selected_cts,
      server = FALSE
    )

    if (length(selected_cts) == 0) {
      igv_loaded(FALSE)
      status_msg("No valid cell types selected for IGV.")
      return(FALSE)
    }

    igv_loaded(FALSE)
    status_msg(paste0(reason, " for ", gene_label, " ..."))

    ok <- tryCatch({
      status_msg(
        paste0(
          "Generating IGV track files",
          "\nGene: ", gene_label,
          "\nCell types: ", paste(selected_cts, collapse = ", ")
        )
      )

      cat("[IGV R] STEP A: ensure_gene_tracks\n")
      flush.console()

      ensure_gene_tracks(gene, selected_cts)

      cat("[IGV R] STEP B: send_cfg\n")
      flush.console()

      send_cfg(gene, selected_cts)

      cat("[IGV R] STEP C: configuration sent\n")
      flush.console()

      igv_current_gene(gene)

      # Keep FALSE until JavaScript confirms successful browser creation.
      igv_loaded(FALSE)

      status_msg(
        paste0(
          "IGV configuration sent; waiting for JavaScript",
          "\nGene: ", gene_label,
          "\nCell types: ", length(selected_cts)
        )
      )

      TRUE
    }, error = function(e) {
      error_message <- conditionMessage(e)

      cat(
        "\n[IGV R] load_igv_for_gene FAILED\n",
        "[IGV R] Error: ", error_message, "\n",
        sep = ""
      )
      flush.console()

      igv_loaded(FALSE)

      status_msg(
        paste0(
          "R-side IGV preparation failed",
          "\nGene: ", gene_label,
          "\nError: ", error_message
        )
      )

      FALSE
    })

    ok
  }

  

  pastaa_status_msg <- reactiveVal("")
  output$pastaa_status <- renderText(pastaa_status_msg())

  # ===== CHANGED: shared SHAP bounds across selected IGV tracks ==============
  calc_track_ylim <- function(gene, selected_cts, min_y = 1e-6) {
    all_vals <- unlist(
      lapply(selected_cts, function(ct) {
        df <- tryCatch(
          load_igv_shap_df(gene, ct),
          error = function(e) NULL
        )

        if (is.null(df) || nrow(df) == 0) {
          return(numeric(0))
        }

        vals <- suppressWarnings(as.numeric(df$score))
        vals[is.finite(vals)]
      }),
      use.names = FALSE
    )

    if (length(all_vals) == 0) {
      return(min_y)
    }

    y <- max(abs(all_vals), na.rm = TRUE)

    if (!is.finite(y) || y <= 0) {
      y <- min_y
    }

    max(as.numeric(y), min_y)
  }

  calc_atac_ylim <- function(gene, selected_cts, min_y = 1e-6) {
    all_vals <- unlist(
      lapply(selected_cts, function(ct) {
        df <- tryCatch(
          load_atac_df(gene, ct),
          error = function(e) NULL
        )

        if (is.null(df) || nrow(df) == 0) {
          return(numeric(0))
        }

        vals <- suppressWarnings(as.numeric(df$score))
        vals[is.finite(vals)]
      }),
      use.names = FALSE
    )

    if (length(all_vals) == 0) {
      return(min_y)
    }

    y <- max(all_vals, na.rm = TRUE)

    if (!is.finite(y) || y <= 0) {
      y <- min_y
    }

    max(as.numeric(y), min_y)
  }
  # ===== END CHANGED =========================================================

  calc_gene_locus_from_shap <- function(gene, pad = 50000) {
    best <- NULL

    for (ct in igv_cell_types) {
      df <- tryCatch(
        load_igv_shap_df(gene, ct),
        error = function(e) NULL
      )

      if (is.null(df) || nrow(df) == 0) next

      df <- df[, c("interaction", "score"), drop = FALSE]

      coords <- parse_interaction_safe(df$interaction)
      coords <- coords[
        !is.na(coords$chrom) &
          !is.na(coords$start) &
          !is.na(coords$end),
        ,
        drop = FALSE
      ]

      if (nrow(coords) == 0) next

      chrom <- coords$chrom[1]
      s <- min(coords$start, na.rm = TRUE)
      e <- max(coords$end, na.rm = TRUE)
      span <- e - s

      if (is.null(best) || span > best$span) {
        best <- list(chrom = chrom, start = s, end = e, span = span)
      }
    }

    if (is.null(best)) return(NULL)

    start <- max(0, best$start - pad)
    end <- best$end + pad

    if (!grepl("^chr", best$chrom)) {
      best$chrom <- paste0("chr", best$chrom)
    }

    sprintf("%s:%d-%d", best$chrom, start, end)
  }

  write_gene_bedgraph <- function(gene, ct, out_bg, max_rows = 200000) {
    df <- tryCatch(
      load_igv_shap_df(gene, ct),
      error = function(e) NULL
    )

    if (is.null(df) || nrow(df) == 0) return(FALSE)

    df <- df[, c("interaction", "score"), drop = FALSE]

    coords <- parse_interaction_safe(df$interaction)
    df2 <- cbind(coords, score = df$score)

    df2 <- df2 |>
      dplyr::mutate(
        chrom = as.character(chrom),
        chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
        start = suppressWarnings(as.integer(start)),
        end = suppressWarnings(as.integer(end)),
        score = suppressWarnings(as.numeric(score))
      ) |>
      dplyr::filter(
        !is.na(chrom),
        !is.na(start),
        !is.na(end),
        is.finite(score)
      ) |>
      dplyr::mutate(
        start2 = pmin(start, end),
        end2 = pmax(start, end)
      ) |>
      dplyr::select(chrom, start = start2, end = end2, score) |>
      dplyr::filter(end > start) |>
      dplyr::arrange(chrom, start, end)

    if (nrow(df2) == 0) return(FALSE)

    if (nrow(df2) > max_rows) {
      ord <- order(abs(df2$score), decreasing = TRUE)
      df2 <- df2[ord[seq_len(max_rows)], , drop = FALSE]
    }

    out <- df2[, c("chrom", "start", "end", "score")]

    write.table(
      out,
      out_bg,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    TRUE
  }

  write_gene_atac_bedgraph <- function(gene, ct, out_bg, max_rows = 200000) {
    df <- tryCatch(
      load_atac_df(gene, ct),
      error = function(e) NULL
    )

    if (is.null(df) || nrow(df) == 0) return(FALSE)

    df <- df[, c("interaction", "score"), drop = FALSE]

    coords <- parse_interaction_safe(df$interaction)
    df2 <- cbind(coords, score = df$score)

    df2 <- df2 |>
      dplyr::mutate(
        chrom = as.character(chrom),
        chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
        start = suppressWarnings(as.integer(start)),
        end = suppressWarnings(as.integer(end)),
        score = suppressWarnings(as.numeric(score))
      ) |>
      dplyr::filter(
        !is.na(chrom),
        !is.na(start),
        !is.na(end),
        is.finite(score)
      ) |>
      dplyr::mutate(
        start2 = pmin(start, end),
        end2 = pmax(start, end)
      ) |>
      dplyr::select(chrom, start = start2, end = end2, score) |>
      dplyr::filter(end > start) |>
      dplyr::arrange(chrom, start, end)

    if (nrow(df2) == 0) return(FALSE)

    if (nrow(df2) > max_rows) {
      ord <- order(abs(df2$score), decreasing = TRUE)
      df2 <- df2[ord[seq_len(max_rows)], , drop = FALSE]
    }

    write.table(
      df2,
      out_bg,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    TRUE
  }

  write_gene_overview_regions_bedgraph <- function(
    gene,
    out_bg,
    top_n_per_ct = 200,
    gap_bp = 1
  ) {
    all_regions <- list()

    for (ct in igv_cell_types) {
      df <- tryCatch(
        load_igv_shap_df(gene, ct),
        error = function(e) NULL
      )

      if (is.null(df) || nrow(df) == 0) next

      df <- df |>
        dplyr::mutate(abs_score = abs(score)) |>
        dplyr::filter(is.finite(abs_score)) |>
        dplyr::arrange(dplyr::desc(abs_score)) |>
        dplyr::slice_head(n = top_n_per_ct)

      if (nrow(df) == 0) next

      coords <- parse_interaction_safe(df$interaction)
      coords$score <- df$abs_score

      coords <- coords[
        !is.na(coords$chrom) &
          !is.na(coords$start) &
          !is.na(coords$end) &
          is.finite(coords$score),
        ,
        drop = FALSE
      ]

      if (nrow(coords) > 0) {
        all_regions[[ct]] <- coords[, c("chrom", "start", "end"), drop = FALSE]
      }
    }

    if (length(all_regions) == 0) {
      return(FALSE)
    }

    bed <- dplyr::bind_rows(all_regions) |>
      dplyr::mutate(
        chrom = as.character(chrom),
        chrom = ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom)),
        start = suppressWarnings(as.integer(start)),
        end = suppressWarnings(as.integer(end))
      ) |>
      dplyr::filter(
        !is.na(chrom),
        !is.na(start),
        !is.na(end)
      ) |>
      dplyr::mutate(
        start2 = pmin(start, end),
        end2 = pmax(start, end)
      ) |>
      dplyr::select(chrom, start = start2, end = end2) |>
      dplyr::filter(end > start) |>
      dplyr::distinct(chrom, start, end) |>
      dplyr::arrange(chrom, start, end)

    if (nrow(bed) == 0) {
      return(FALSE)
    }

    merge_one_chrom <- function(starts, ends, gap_bp = 1) {
      n <- length(starts)

      if (n == 0) {
        return(data.frame(start = integer(0), end = integer(0)))
      }

      out_s <- integer(0)
      out_e <- integer(0)

      cur_s <- starts[1]
      cur_e <- ends[1]

      if (n > 1) {
        for (i in 2:n) {
          if (starts[i] <= (cur_e + gap_bp)) {
            if (ends[i] > cur_e) cur_e <- ends[i]
          } else {
            out_s <- c(out_s, cur_s)
            out_e <- c(out_e, cur_e)
            cur_s <- starts[i]
            cur_e <- ends[i]
          }
        }
      }

      out_s <- c(out_s, cur_s)
      out_e <- c(out_e, cur_e)

      data.frame(start = out_s, end = out_e)
    }

    splitted <- split(bed, bed$chrom)

    merged_list <- lapply(names(splitted), function(chr) {
      dfc <- splitted[[chr]]
      m <- merge_one_chrom(dfc$start, dfc$end, gap_bp = gap_bp)

      if (nrow(m) == 0) return(NULL)

      data.frame(
        chrom = chr,
        start = m$start,
        end = m$end,
        value = 1
      )
    })

    merged <- dplyr::bind_rows(merged_list) |>
      dplyr::arrange(chrom, start, end)

    if (is.null(merged) || nrow(merged) == 0) {
      return(FALSE)
    }

    write.table(
      merged[, c("chrom", "start", "end", "value")],
      out_bg,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    TRUE
  }

  

  build_paired_shap_atac_tracks <- function(
    gene,
    selected_cts,
    cache_bust = as.integer(Sys.time())
  ) {
    out <- list()
    order_i <- 10L

    # ===== CHANGED: one symmetric SHAP range for all selected cell types =====
    shap_y <- calc_track_ylim(
      gene = gene,
      selected_cts = selected_cts
    )
    atac_y <- calc_atac_ylim(
      gene = gene,
      selected_cts = selected_cts
    )

    cat(
      "[IGV R] Shared Feature Importance scale",
      "| gene:", gene,
      "| cell types:", paste(selected_cts, collapse = ", "),
      "| range:", paste0("-", signif(shap_y, 4), " to ", signif(shap_y, 4)),
      "\n"
    )
    # ===== END CHANGED ========================================================

    for (ct in selected_cts) {
      ct_file <- ct_to_file(ct)
      col <- get_ct_color(ct, fallback = "#666666")

      shap_fname <- sprintf("SHAP_%s_%s.bedGraph", gene, ct_file)
      atac_fname <- sprintf("ATAC_%s_%s.bedGraph", gene, ct_file)

      shap_path <- file.path(bw_dir, shap_fname)
      atac_path <- file.path(bw_dir, atac_fname)

      if (file.exists(shap_path)) {
        out[[length(out) + 1L]] <- list(
          id = paste0("shap__", ct_file),
          name = paste0(
            "Feature Importance (",
            input$feature_importance_type,
            ") - ",
            ct
          ),
          type = "wig",
          format = "bedgraph",
          url = sprintf("shap_bw/%s?v=%s", shap_fname, cache_bust),
          height = 150,
          autoscale = TRUE,
          min = -shap_y,
          max = shap_y,
          color = col,
          altColor = col,
          labelColor = col,
          visibilityWindow = 1e9,
          order = order_i
        )

        order_i <- order_i + 1L
      }

      if (file.exists(atac_path)) {
        out[[length(out) + 1L]] <- list(
          id = paste0("atac__", ct_file),
          name = paste0(
            "Epigenetic Signal (",
            input$epigenetic_signal_type,
            ") - ",
            ct
          ),
          type = "wig",
          format = "bedgraph",
          url = sprintf("shap_bw/%s?v=%s", atac_fname, cache_bust),
          height = 115,
          autoscale = TRUE,
          min = 0,
          max = atac_y,
          color = col,
          altColor = col,
          labelColor = col,
          visibilityWindow = 1e9,
          order = order_i
        )

        order_i <- order_i + 1L
      }
    }

    if (length(out) == 0) {
      stop(
        "No existing IGV track files found for gene ",
        gene,
        "."
      )
    }

    unname(out)
  }
  write_empty_regulatory_bed <- function(out_bed) {
    writeLines(character(0), con = out_bed)
  }

  ensure_gene_tracks <- function(gene, selected_cts) {
    dir.create(bw_dir, showWarnings = FALSE, recursive = TRUE)

    for (ct in selected_cts) {
      shap_bg <- file.path(
        bw_dir,
        sprintf("SHAP_%s_%s.bedGraph", gene, ct_to_file(ct))
      )

      status_msg(
        paste0(
          "Generating Feature Importance bedGraph: ",
          gene,
          " / ",
          ct,
          " ..."
        )
      )

      if (file.exists(shap_bg)) {
        unlink(shap_bg)
      }

      write_gene_bedgraph(gene, ct, shap_bg, max_rows = 200000)

      atac_bg <- file.path(
        bw_dir,
        sprintf("ATAC_%s_%s.bedGraph", gene, ct_to_file(ct))
      )

      if (!file.exists(atac_bg)) {
        status_msg(
          paste0(
            "Generating Epigenetic Signal bedGraph: ",
            gene,
            " / ",
            ct,
            " ..."
          )
        )
        write_gene_atac_bedgraph(gene, ct, atac_bg, max_rows = 200000)
      }
    }

    
  }

  normalise_regulatory_build_df <- function(reg) {
    reg <- as.data.frame(reg, stringsAsFactors = FALSE)

    if (nrow(reg) == 0 || ncol(reg) < 3) {
      return(data.frame())
    }

    names_lower <- tolower(colnames(reg))

    find_col <- function(candidates) {
      index <- match(candidates, names_lower, nomatch = 0L)
      index <- index[index > 0L]

      if (length(index) == 0) {
        return(NA_integer_)
      }

      index[1]
    }

    chrom_col <- find_col(c(
      "chrom", "chromosome_name",
      "seq_region_name", "chromosome"
    ))

    start_col <- find_col(c(
      "start", "chromosome_start"
    ))

    end_col <- find_col(c(
      "end", "chromosome_end"
    ))

    name_col <- find_col(c(
      "name",
      "feature_type_name",
      "regulatory_feature_type_name",
      "regulatory_feature_type",
      "feature_type"
    ))

    if (
      is.na(chrom_col) ||
      is.na(start_col) ||
      is.na(end_col)
    ) {
      return(data.frame())
    }

    region_name <- if (!is.na(name_col)) {
      as.character(reg[[name_col]])
    } else {
      rep("regulatory_region", nrow(reg))
    }

    region_name <- sub(";.*$", "", region_name)
    region_name <- trimws(region_name)

    invalid <- is.na(region_name) | !nzchar(region_name) | region_name == "."
    region_name[invalid] <- "regulatory_region"

    out <- data.frame(
      chrom = as.character(reg[[chrom_col]]),
      start = suppressWarnings(as.integer(reg[[start_col]])),
      end = suppressWarnings(as.integer(reg[[end_col]])),
      name = region_name,
      stringsAsFactors = FALSE
    )

    out$chrom <- ifelse(
      grepl("^chr", out$chrom, ignore.case = TRUE),
      out$chrom,
      paste0("chr", out$chrom)
    )

    out |>
      dplyr::filter(
        !is.na(start),
        !is.na(end),
        end > start
      ) |>
      dplyr::distinct(chrom, start, end, name)
  }


  read_local_regulatory_build <- function() {
    if (!file.exists(regulatory_build_bed)) {
      return(data.frame())
    }

    reg <- tryCatch(
      read.table(
        regulatory_build_bed,
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = "",
        fill = TRUE
      ),
      error = function(e) {
        status_msg(
          paste("Could not read local regulatory BED:", e$message)
        )
        data.frame()
      }
    )

    if (nrow(reg) == 0 || ncol(reg) < 3) {
      return(data.frame())
    }

    colnames(reg)[1:3] <- c("chrom", "start", "end")

    if (ncol(reg) >= 4) {
      colnames(reg)[4] <- "name"
    } else {
      reg$name <- "regulatory_region"
    }

    cat(
      "[REGULATORY DEBUG] raw local BED:",
      nrow(reg), "rows |",
      ncol(reg), "cols\n"
    )

    norm_reg <- normalise_regulatory_build_df(reg)

    cat(
      "[REGULATORY DEBUG] after normalisation:",
      nrow(norm_reg), "rows\n"
    )

    if (nrow(norm_reg) > 0) {
      cat(
        "[REGULATORY DEBUG] chromosomes:",
        paste(head(unique(norm_reg$chrom), 10), collapse = ", "),
        "\n"
      )

      cat("[REGULATORY DEBUG] first normalised row:\n")
      print(head(norm_reg, 1))
    }

    norm_reg
  }

  find_regulatory_dataset <- function(species_name) {

    if (
      is.null(species_name) ||
      is.na(species_name) ||
      !nzchar(trimws(species_name))
    ) {
      status_msg(
        "Species information is missing from the EpIVoi object. Regulatory build cannot be queried."
      )
      return(NULL)
    }

    regulation_mart <- tryCatch(
      biomaRt::useEnsembl(
        biomart = "regulation"
      ),
      error = function(e) {
        status_msg(
          paste(
            "Could not connect to Ensembl Regulation BioMart:",
            e$message
          )
        )
        NULL
      }
    )

    if (is.null(regulation_mart)) {
      return(NULL)
    }

    datasets <- tryCatch(
      biomaRt::listDatasets(regulation_mart),
      error = function(e) {
        status_msg(
          paste(
            "Could not retrieve Ensembl Regulation datasets:",
            e$message
          )
        )
        data.frame()
      }
    )

    if (nrow(datasets) == 0) {
      return(NULL)
    }

    species_name <- trimws(as.character(species_name)[1])

    # First try to match the scientific species name in the
    # BioMart dataset description.
    description_match <- grepl(
      species_name,
      datasets$description,
      ignore.case = TRUE,
      fixed = TRUE
    )

    candidates <- datasets[description_match, , drop = FALSE]

    # Fallback: derive the standard Ensembl species prefix.
    # Example:
    # Homo sapiens -> hsapiens
    # Mus musculus -> mmusculus
    # Danio rerio -> drerio
    if (nrow(candidates) == 0) {

      species_parts <- strsplit(
        tolower(species_name),
        "\\s+"
      )[[1]]

      if (length(species_parts) >= 2) {

        ensembl_prefix <- paste0(
          substr(species_parts[1], 1, 1),
          species_parts[2]
        )

        dataset_match <- grepl(
          paste0("^", ensembl_prefix, "_"),
          datasets$dataset,
          ignore.case = TRUE
        )

        candidates <- datasets[
          dataset_match,
          ,
          drop = FALSE
        ]
      }
    }

    if (nrow(candidates) == 0) {
      status_msg(
        paste0(
          "No Ensembl Regulation dataset is available for species: ",
          species_name,
          ". Regulatory build track will not be shown."
        )
      )
      return(NULL)
    }

    # Prefer an actual regulatory-feature dataset if several datasets match.
    regulatory_candidates <- candidates[
      grepl(
        "regulatory",
        candidates$dataset,
        ignore.case = TRUE
      ),
      ,
      drop = FALSE
    ]

    if (nrow(regulatory_candidates) > 0) {
      candidates <- regulatory_candidates
    }

    dataset <- as.character(candidates$dataset[1])

    status_msg(
      paste0(
        "Using Ensembl Regulation dataset ",
        dataset,
        " for ",
        species_name,
        "."
      )
    )

    dataset
  }

  load_regulatory_build_from_ensembl_ftp <- function() {

    ensembl_release <- 116L

    species_slug <- tolower(gsub(" ", "_", species_name))

    assembly <- dplyr::case_when(
      genome_version == "hg38" ~ "GRCh38",
      genome_version == "hg19" ~ "GRCh37",
      genome_version == "mm39" ~ "GRCm39",
      genome_version == "mm10" ~ "GRCm38",
      genome_version == "danRer11" ~ "GRCz11",
      TRUE ~ genome_version
    )

    species_file <- paste0(
      toupper(substr(species_slug, 1, 1)),
      substr(species_slug, 2, nchar(species_slug))
    )

    ftp_url <- paste0(
      "https://ftp.ensembl.org/pub/release-",
      ensembl_release,
      "/regulation/",
      species_slug,
      "/",
      assembly,
      "/annotation/",
      species_file,
      ".",
      assembly,
      ".regulatory_features.v",
      ensembl_release,
      ".gff3.gz"
    )

    status_msg(
      paste0(
        "Loading Ensembl regulatory annotation for ",
        species_name,
        " (",
        assembly,
        ")..."
      )
    )

    message("[REGULATORY] Ensembl FTP URL: ", ftp_url)

    cache_dir <- file.path(
      tempdir(),
      "epivoi_regulatory_cache"
    )

    dir.create(
      cache_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    cache_file <- file.path(
      cache_dir,
      paste0(
        species_slug,
        "_",
        assembly,
        "_regulatory_features_v",
        ensembl_release,
        ".gff3.gz"
      )
    )

    if (!file.exists(cache_file)) {

      download_ok <- tryCatch(
        {
          utils::download.file(
            url = ftp_url,
            destfile = cache_file,
            mode = "wb",
            quiet = TRUE
          )
          TRUE
        },
        error = function(e) {
          status_msg(
            paste(
              "Could not download Ensembl regulatory annotation:",
              e$message
            )
          )
          FALSE
        }
      )

      if (!isTRUE(download_ok)) {
        if (file.exists(cache_file)) {
          unlink(cache_file)
        }
        return(data.frame())
      }
    }

    reg_raw <- tryCatch(
      utils::read.delim(
        gzfile(cache_file),
        header = FALSE,
        sep = "\t",
        comment.char = "#",
        quote = "",
        stringsAsFactors = FALSE
      ),
      error = function(e) {
        status_msg(
          paste(
            "Could not read Ensembl regulatory GFF3:",
            e$message
          )
        )
        data.frame()
      }
    )

    if (nrow(reg_raw) == 0 || ncol(reg_raw) < 5) {
      status_msg(
        "Ensembl regulatory GFF3 contains no usable regions."
      )
      return(data.frame())
    }

    reg <- data.frame(
      chrom = as.character(reg_raw[[1]]),
      start = suppressWarnings(
        as.integer(reg_raw[[4]]) - 1L
      ),
      end = suppressWarnings(
        as.integer(reg_raw[[5]])
      ),
      name = as.character(reg_raw[[3]]),
      stringsAsFactors = FALSE
    )

    reg$chrom <- ifelse(
      grepl("^chr", reg$chrom, ignore.case = TRUE),
      reg$chrom,
      paste0("chr", reg$chrom)
    )

    reg <- reg |>
      dplyr::filter(
        !is.na(start),
        !is.na(end),
        start >= 0,
        end > start
      ) |>
      dplyr::distinct(
        chrom,
        start,
        end,
        name
      )

    status_msg(
      paste0(
        "Ensembl regulatory annotation loaded and cached: ",
        nrow(reg),
        " regions."
      )
    )

    message(
      "[REGULATORY] Loaded from Ensembl FTP: ",
      nrow(reg),
      " regions"
    )

    reg
  }


  get_regulatory_build_for_session <- function() {
    cache_key <- paste(species_name, genome_version, sep = "__")

    cache <- regulatory_build_cache()

    if (
      !is.null(cache[[cache_key]]) &&
      is.data.frame(cache[[cache_key]]) &&
      nrow(cache[[cache_key]]) > 0
    ) {
      return(cache[[cache_key]])
    }

    reg <- data.frame()

    if (nzchar(regulatory_build_bed) && file.exists(regulatory_build_bed)) {
      reg <- read_local_regulatory_build()

      if (nrow(reg) > 0) {
        status_msg(
          paste0(
            "Regulatory build loaded from local BED and cached for this session: ",
            nrow(reg),
            " regions."
          )
        )
      }
    }

    if (nrow(reg) == 0) {
      reg <- load_regulatory_build_from_ensembl_ftp()
    }

    if (is.data.frame(reg) && nrow(reg) > 0) {
      cache[[cache_key]] <- reg
      regulatory_build_cache(cache)
    }

    reg
  }


  load_regulatory_build_for_locus <- function(locus) {
    loc <- parse_locus_string(locus)

    if (is.null(loc)) {
      return(data.frame())
    }

    reg <- get_regulatory_build_for_session()

    if (nrow(reg) == 0) {
      return(data.frame())
    }

    reg |>
      dplyr::filter(
        chrom == loc$chrom,
        end >= loc$start,
        start <= loc$end
      )
  }


  write_regulatory_bed_for_locus <- function(locus, out_bed) {
    cat("\n[REGULATORY] Requested locus:", locus, "\n")
    cat("[REGULATORY] Source BED:", regulatory_build_bed, "\n")
    cat("[REGULATORY] Source exists:", file.exists(regulatory_build_bed), "\n")

    reg_all <- get_regulatory_build_for_session()

    cat(
      "[REGULATORY] Total regulatory regions loaded:",
      nrow(reg_all),
      "\n"
    )

    reg <- load_regulatory_build_for_locus(locus)

    cat(
      "[REGULATORY] Regions overlapping current locus:",
      nrow(reg),
      "\n"
    )

    if (nrow(reg) == 0) {
      cat("[REGULATORY] No overlapping regions -> no IGV track created\n")

      writeLines(character(0), con = out_bed)
      return(FALSE)
    }

    bed <- reg |>
      dplyr::transmute(
        chrom = chrom,
        start = as.integer(start),
        end = as.integer(end),
        name = ifelse(
          !is.na(name) & nzchar(name),
          name,
          "regulatory_region"
        ),
        score = 0,
        strand = "."
      )

    write.table(
      bed,
      out_bed,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )

    cat("[REGULATORY] Output BED:", out_bed, "\n")
    cat("[REGULATORY] Output exists:", file.exists(out_bed), "\n")
    cat("[REGULATORY] Output size:", file.info(out_bed)$size, "bytes\n")
    cat("[REGULATORY] BED rows written:", nrow(bed), "\n")

    if (nrow(bed) > 0) {
      cat("[REGULATORY] First BED row:\n")
      print(utils::head(bed, 1))
    }

    TRUE
  }

  send_cfg <- function(gene, selected_cts) {
    locus <- calc_gene_locus_from_shap(gene, pad = 50000)
    if (is.null(locus) || !nzchar(locus)) locus <- gene

    cache_bust <- as.integer(Sys.time())

    refseq_url <- if (genome_version == "hg19") {
      sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz?v=%d", cache_bust)
    } else {
      sprintf("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz?v=%d", cache_bust)
    }

    refseq_track <- list(
      id = "refseq",
      name = "Genes / RefSeq",
      type = "annotation",
      format = "refgene",
      url = refseq_url,
      displayMode = "EXPANDED",
      height = 130,
      visibilityWindow = 1e9,
      order = 1L
    )
    signal_tracks <- build_paired_shap_atac_tracks(
      gene = gene,
      selected_cts = selected_cts,
      cache_bust = cache_bust
    )
    reg_tracks <- list()

    if (isTRUE(isolate(input$show_regulatory_track))) {
      reg_bed_name <- sprintf("REGULATORY_%s.bed", gene)
      reg_bed_path <- file.path(bw_dir, reg_bed_name)

      ok_reg <- write_regulatory_bed_for_locus(locus, reg_bed_path)

      if (isTRUE(ok_reg)) {
        cat(
          "[REGULATORY] IGV URL:",
          sprintf(
            "shap_bw/%s?v=%s",
            reg_bed_name,
            cache_bust
          ),
          "\n"
        )
        reg_tracks <- list(
          list(
            id = "ensembl_regulatory",
            name = "Ensembl Regulatory Build",
            type = "annotation",
            format = "bed",
            url = sprintf(
              "shap_bw/%s?v=%s",
              reg_bed_name,
              cache_bust
            ),
            displayMode = "EXPANDED",
            height = 80,
            visibilityWindow = 1e9,
            order = 5L
          )
        )

        status_msg("Regulatory build track added.")
      } else {
        reg_tracks <- list()

        status_msg(
          paste(
            "No regulatory regions could be loaded.",
            "The regulatory track will not be shown."
          )
        )
      }
      
    }

    status_msg(
      paste0(
        "DEBUG first IGV track URL: ",
        signal_tracks[[1]]$url
      )
    )

    fimo_files <- list.files(
      fimo_www_dir,
      pattern = paste0("^", gene, "_.*_fimo_sites\\.bed$"),
      full.names = FALSE
    )

    fimo_tracks <- lapply(seq_along(fimo_files), function(i) {
      f <- fimo_files[[i]]
      name <- sub(paste0("^", gene, "_"), "", f)
      name <- sub("_fimo_sites\\.bed$", "", name)

      list(
        id = paste0("fimo__", i),
        name = paste0("FIMO - ", name),
        type = "annotation",
        format = "bed",
        url = paste0("fimo/", f, "?v=", cache_bust),
        displayMode = "EXPANDED",
        height = 80,
        visibilityWindow = 1e9,
        order = 1000L + i
      )
    })
    
    cfg <- list(
      genome = genome_version,
      locus = locus,
      tracks = c(
        list(refseq_track),
        signal_tracks,
        reg_tracks
      )
    )

    cat(
      "[IGV R] Sending igv-create message",
      "| genome:", cfg$genome,
      "| locus:", cfg$locus,
      "| tracks:", length(cfg$tracks),
      "\n"
    )

    session$sendCustomMessage("igv-create", cfg)
    status_msg(sprintf("Loading IGV... gene=%s | %d cell types", gene, length(selected_cts)))
  }

  output$igv_hint <- renderUI({
    gene_resolved <- resolve_igv_gene_input(input$igv_gene)
    if (is.null(gene_resolved)) return(NULL)

    gene <- gene_resolved$gene_id
    gene_label <- gene_resolved$gene_display

    locus_txt <- calc_gene_locus_from_shap(gene, pad = 0)
    if (is.null(locus_txt) || !nzchar(locus_txt)) return(NULL)

    tags$div(
      tags$small(
        style = "color:#555; display:block; margin-bottom:6px;",
        paste("Strong Feature Importance signal for", gene_label, "is typically observed around")
      ),
      actionLink("igv_jump", locus_txt),
      tags$small(
        style = "color:#777; display:block; margin-top:6px; margin-bottom:8px;",
        "(click to jump; zoom out to kb scale if needed)"
      )
    )
  })

  observeEvent(TRUE, {
    default_gene <- if ("ABHD5" %in% igv_gene_list) {
      "ABHD5"
    } else if ("NFATC2" %in% igv_gene_list) {
      "NFATC2"
    } else if (length(igv_gene_list) > 0) {
      igv_gene_list[1]
    } else {
      NULL
    }
    default_cts <- if (length(igv_cell_types) > 0) {
      sort(igv_cell_types)[1]
    } else {
      character(0)
    }
    updateSelectizeInput(
      session,
      "igv_celltypes",
      choices = igv_cell_types,
      selected = default_cts,
      server = FALSE
    )

    updateSelectizeInput(
      session,
      "igv_gene",
      choices = igv_gene_list,
      selected = default_gene,
      server = TRUE
    )

    updateSelectizeInput(
      session,
      "pastaa_celltypes",
      choices = igv_cell_types,
      selected = NULL,
      server = FALSE
    )

    updateSelectizeInput(
      session,
      "pastaa_gene",
      choices = igv_gene_list,
      selected = default_gene,
      server = TRUE
    )
  }, once = TRUE)

  observeEvent(input$igv_jump, {
    gene_resolved <- resolve_igv_gene_input(input$igv_gene)

    if (is.null(gene_resolved)) {
      status_msg(paste("Gene not found:", input$igv_gene))
      return(NULL)
    }

    gene <- gene_resolved$gene_id
    locus <- calc_gene_locus_from_shap(gene, pad = 50000)
    req(locus)

    status_msg(paste("Jumping to", locus))
    session$sendCustomMessage("igv-search", list(locus = locus))
  })





  

  observeEvent(input$igv_back_gene, {
    prev_gene <- igv_previous_gene()

    if (is.null(prev_gene) || !nzchar(prev_gene)) {
      status_msg("No previous gene stored yet.")
      return(NULL)
    }

    gene_label <- id2symbol[[prev_gene]]

    if (is.null(gene_label) || is.na(gene_label) || !nzchar(gene_label)) {
      gene_label <- prev_gene
    }

    h <- igv_gene_history()

    selected_cts <- NULL

    if (!is.null(h[[prev_gene]]) && !is.null(h[[prev_gene]]$celltypes)) {
      selected_cts <- h[[prev_gene]]$celltypes
    }

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      selected_cts <- if (length(igv_cell_types) > 0) {
        sort(igv_cell_types)[1]
      } else {
        character(0)
      }
    }

    updateSelectizeInput(
      session,
      "igv_gene",
      choices = igv_gene_list,
      selected = gene_label,
      server = TRUE
    )

    updateSelectizeInput(
      session,
      "igv_celltypes",
      choices = igv_cell_types,
      selected = selected_cts,
      server = FALSE
    )

    restore_igv_state_for_gene(prev_gene)

    load_igv_for_gene(
      gene_input = gene_label,
      selected_cts = selected_cts,
      reason = "Back to previous gene"
    )
  })

  observeEvent(input$igv_reload, {
    cat("\n[IGV R] Reload button clicked\n")

    gene_input <- isolate(input$igv_gene)
    selected_cts <- isolate(input$igv_celltypes)

    cat("[IGV R] gene input:", gene_input, "\n")
    cat(
      "[IGV R] selected cell types:",
      paste(selected_cts, collapse = ", "),
      "\n"
    )

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      selected_cts <- if (length(igv_cell_types) > 0) {
        sort(igv_cell_types)[1]
      } else {
        character(0)
      }
    }

    result <- load_igv_for_gene(
      gene_input = gene_input,
      selected_cts = selected_cts,
      reason = "Reloading IGV"
    )

    cat("[IGV R] load_igv_for_gene result:", result, "\n")
  }, ignoreInit = TRUE)

  
  

  observeEvent(input$run_pastaa, {

    shinyjs::disable("run_pastaa")
    on.exit(shinyjs::enable("run_pastaa"), add = TRUE)
    if (!is.finite(input$pastaa_n_regions) || input$pastaa_n_regions < 10) {
      pastaa_status_msg("Number of top regions must be at least 10.")
      return(NULL)
    }

    if (
      !is.finite(input$pastaa_save_top_n) ||
      input$pastaa_save_top_n < 1 ||
      input$pastaa_save_top_n > 20
    ) {
      pastaa_status_msg("Save top N TFs must be between 1 and 20.")
      return(NULL)
    }

    if (!is.finite(input$pastaa_q_cutoff) ||
        input$pastaa_q_cutoff < 0 ||
        input$pastaa_q_cutoff > 1) {
      pastaa_status_msg("q-value cutoff must be between 0 and 1.")
      return(NULL)
    }
    run_mode <- input$pastaa_run_mode

    if (identical(run_mode, "gene")) {
      gene_resolved <- resolve_igv_gene_input(input$pastaa_gene)

      if (is.null(gene_resolved)) {
        pastaa_status_msg(paste("Gene not found:", input$pastaa_gene))
        return(NULL)
      }

      gene <- gene_resolved$gene_id
      gene_label <- gene_resolved$gene_display
    } else {
      gene <- NULL
      gene_label <- "all eligible genes"
    }

    filter_performance <- identical(run_mode, "celltype") &&
      isTRUE(input$pastaa_filter_performance)

    min_test_correlation <- suppressWarnings(
      as.numeric(input$pastaa_min_test_correlation)[1]
    )
    performance_metric <- input$pastaa_performance_metric %||% default_test_column

    if (
      filter_performance &&
      (!is.finite(min_test_correlation) ||
       min_test_correlation < -1 ||
       min_test_correlation > 1)
    ) {
      pastaa_status_msg(
        "Minimum test correlation must be between -1 and 1."
      )
      return(NULL)
    }

    selected_cts <- isolate(input$pastaa_celltypes)

    if (is.null(selected_cts) || length(selected_cts) == 0) {
      pastaa_status_msg("Please select at least one cell type.")
      return(NULL)
    }

    directions <- switch(
      input$pastaa_direction,
      "pos" = "pos",
      "neg" = "neg",
      "both" = c("pos", "neg"),
      "abs" = "abs",
      c("pos", "neg")
    )

    n_regions <- input$pastaa_n_regions

    use_abs_shap <- identical(
      input$pastaa_direction,
      "abs"
    )



    
    pastaa_status_msg(
      paste0(
        if (identical(run_mode, "gene")) {
          paste0("Running gene-specific PASTAA for ", gene_label)
        } else {
          paste0(
            "Running cell-type PASTAA across ",
            if (filter_performance) {
              paste0("genes with ", performance_metric, " >= ", min_test_correlation)
            } else {
              "all genes"
            }
          )
        },
        " | ",
        length(selected_cts),
        " cell type(s) ..."
      )
    )

    out <- list()
    errors <- character(0)

    withProgress(message = "Running PASTAA", value = 0, {
      total <- length(selected_cts) * length(directions)
      k <- 0

      for (ct in selected_cts) {
        for (direction in directions) {
          k <- k + 1

          

          

          res <- tryCatch(
            run_pastaa_for_gene_ct_direction(
              gene = gene,
              ct = ct,
              direction = direction,
              n_regions = n_regions,
              run_mode = run_mode,
              use_abs_shap = use_abs_shap,
              filter_performance = filter_performance,
              performance_metric = performance_metric,
              min_test_correlation = min_test_correlation
            ),
            error = function(e) e
          )

          incProgress(
            1 / total,
            detail = paste(gene_label, ct, direction)
          )

          if (inherits(res, "error")) {
            errors <- c(errors, paste(ct, direction, ":", res$message))
          } else {
            out[[res$run_id]] <- res
          }
        }
      }
    })

    previous_results <- pastaa_results()

    if (length(previous_results) == 0) {
      combined_results <- out
    } else {
      combined_results <- previous_results
      combined_results[names(out)] <- out
    }

    pastaa_results(combined_results)

    saved_file <- save_top_pastaa_tfs(
      res_list = out,
      top_n = input$pastaa_save_top_n,
      q_cutoff = input$pastaa_q_cutoff
    )

    if (length(out) == 0) {
      pastaa_status_msg(
        paste(
          "PASTAA finished, but no successful results were created.",
          if (length(errors) > 0) paste("Errors:", paste(errors, collapse = " | ")) else ""
        )
      )
      return(NULL)
    }

    if (length(errors) > 0) {
      pastaa_status_msg(
        paste(
          "PASTAA loaded with warnings/errors.",
          "Successful runs:", length(out),
          "| Failed runs:", length(errors),
          "|",
          paste(errors, collapse = " | "),
          if (!is.null(saved_file)) paste("| Top TFs saved:", saved_file) else "| No TFs passed q-value cutoff."
        )
      )
    } else {
      pastaa_status_msg(
        paste(
          "PASTAA loaded successfully.",
          "Successful runs:", length(out),
          "|",
          if (!is.null(saved_file)) paste("Top TFs saved:", saved_file) else "No TFs passed q-value cutoff."
        )
      )
    }
  })

  output$pastaa_result_tabs <- renderUI({
    res_list <- pastaa_results()

    if (length(res_list) == 0) {
      return(tags$em("No PASTAA results yet. Select a gene/cell type and click Start PASTAA."))
    }

    tabs <- lapply(names(res_list), function(nm) {
      tabPanel(
        title = nm,
        tags$p(
          paste0(
            "Regions used: ",
            res_list[[nm]]$n_regions,
            " | Result file: ",
            res_list[[nm]]$sorted_file
          )
        ),
        DTOutput(paste0("pastaa_table_", nm))
      )
    })

    do.call(
      tabsetPanel,
      c(
        list(id = "pastaa_selected_result_tab"),
        tabs
      )
    )
  })

  pastaa_selected_df <- reactive({
    res_list <- pastaa_results()

    if (length(res_list) == 0) {
      return(NULL)
    }

    selected <- input$pastaa_selected_result_tab

    if (
      is.null(selected) ||
      !nzchar(selected) ||
      !selected %in% names(res_list)
    ) {
      selected <- names(res_list)[1]
    }

    df <- res_list[[selected]]$result

    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    df |>
      dplyr::filter(TF != "DUMMY_MOTIF") |>
      dplyr::mutate(
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR)),
        minus_log10_q = -log10(
          pmax(BH_FDR, .Machine$double.xmin)
        ),
        TF = as.character(TF)
      ) |>
      dplyr::filter(
        is.finite(BH_FDR),
        BH_FDR <= input$pastaa_q_cutoff
      )
  })

  pastaa_selected_direction_pair_df <- reactive({
    res_list <- pastaa_results()

    if (length(res_list) == 0) return(NULL)

    selected <- input$pastaa_selected_result_tab

    if (
      is.null(selected) ||
      !nzchar(selected) ||
      !selected %in% names(res_list)
    ) {
      selected <- names(res_list)[1]
    }

    # Same PASTAA run settings, ignoring only direction
    selected_base <- sub(
      "_(pos|neg|abs)_",
      "_DIRECTION_",
      selected
    )

    matching_runs <- names(res_list)[
      sub(
        "_(pos|neg|abs)_",
        "_DIRECTION_",
        names(res_list)
      ) == selected_base
    ]

    if (length(matching_runs) == 0) return(NULL)

    df <- dplyr::bind_rows(
      lapply(matching_runs, function(nm) {
        res_list[[nm]]$result
      })
    )

    if (nrow(df) == 0) return(NULL)

    df |>
      dplyr::filter(
        TF != "DUMMY_MOTIF",
        direction %in% c("pos", "neg")
      ) |>
      dplyr::mutate(
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR)),
        minus_log10_q = -log10(
          pmax(BH_FDR, .Machine$double.xmin)
        ),
        TF = as.character(TF)
      ) |>
      dplyr::filter(
        is.finite(BH_FDR),
        BH_FDR <= input$pastaa_q_cutoff
      )
  })

  pastaa_summary_df <- reactive({
    res_list <- pastaa_results()
    if (length(res_list) == 0) return(NULL)

    df <- dplyr::bind_rows(
      lapply(names(res_list), function(nm) {
        x <- res_list[[nm]]$result
        x$run_id <- nm
        x
      })
    )

    if (nrow(df) == 0) return(NULL)

    df <- df |>
      dplyr::filter(TF != "DUMMY_MOTIF") |>
      dplyr::mutate(
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR)),
        minus_log10_q = -log10(pmax(BH_FDR, .Machine$double.xmin)),
        TF = as.character(TF),
        cell_type = as.character(cell_type),
        direction = as.character(direction),
        run_mode = as.character(run_mode)
      ) |>
      dplyr::filter(
        is.finite(BH_FDR),
        is.finite(minus_log10_q)
      )

    if (!is.null(input$pastaa_q_cutoff) && is.finite(input$pastaa_q_cutoff)) {
      df <- df |>
        dplyr::filter(BH_FDR <= input$pastaa_q_cutoff)
    }

    if (nrow(df) == 0) return(NULL)

    df
  })

  pastaa_all_table <- reactive({
    res_list <- pastaa_results()

    if (length(res_list) == 0) {
      return(data.frame())
    }

    df <- dplyr::bind_rows(
      lapply(names(res_list), function(nm) {
        x <- res_list[[nm]]$result
        x$run_id <- nm
        x
      })
    )

    if (nrow(df) == 0) {
      return(data.frame())
    }

    df |>
      dplyr::filter(TF != "DUMMY_MOTIF") |>
      dplyr::mutate(
        p_value = suppressWarnings(as.numeric(p_value)),
        BH_FDR = suppressWarnings(as.numeric(BH_FDR))
      ) |>
      dplyr::arrange(BH_FDR, p_value)
  })

  pastaa_top_table <- reactive({
    df <- pastaa_all_table()

    if (nrow(df) == 0) {
      return(df)
    }

    df |>
      dplyr::mutate(
        TF_clean = clean_tf_name(TF)
      ) |>
      dplyr::filter(
        is.finite(BH_FDR),
        BH_FDR <= input$pastaa_q_cutoff,
        !is.na(TF_clean),
        nzchar(TF_clean)
      ) |>
      dplyr::group_by(TF_clean) |>
      dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(BH_FDR, p_value) |>
      dplyr::slice_head(n = min(input$pastaa_save_top_n, 20))
  })

  make_pastaa_top_tf_plot <- function(
    df,
    top_n,
    q_cutoff,
    direction
  ) {
    req(!is.null(df), nrow(df) > 0)

    plot_df <- df |>
      dplyr::filter(.data$direction == .env$direction)

    shiny::validate(
      shiny::need(
        nrow(plot_df) > 0,
        paste0(
          "No ",
          switch(
            direction,
            pos = "positive",
            neg = "negative",
            abs = "absolute",
            direction
          ),
          " PASTAA results pass the selected q-value cutoff."
        )
      )
    )

    top_df <- plot_df |>
      dplyr::group_by(TF) |>
      dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(BH_FDR, p_value) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::mutate(
        TF = factor(TF, levels = rev(TF))
      )

    direction_label <- switch(
      direction,
      pos = "Positive Feature Importance",
      neg = "Negative Feature Importance",
      abs = "Absolute Feature Importance",
      direction
    )

    ggplot(
      top_df,
      aes(
        x = TF,
        y = minus_log10_q
      )
    ) +
      geom_col() +
      coord_flip() +
      labs(
        x = "TF",
        y = "-log10(q-value / BH FDR)",
        title = paste0(
          direction_label,
          ": top enriched TFs, q-value cutoff ≤ ",
          q_cutoff
        )
      ) +
      theme_bw(base_size = 14)
  }

  make_pastaa_posneg_combined_plot <- function(
      df,
      top_n,
      q_cutoff,
      max_groups = 10
    ) {
      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          paste0(
            "No TF enrichment results pass the selected q-value cutoff (≤ ",
            q_cutoff,
            ")."
          )
        )
      )

      df <- df |>
        dplyr::filter(direction %in% c("pos", "neg")) |>
        dplyr::mutate(cell_type = as.character(cell_type))

      shiny::validate(
        shiny::need(
          nrow(df) > 0,
          "No positive/negative TF enrichment results are available."
        )
      )

      all_groups <- sort(unique(df$cell_type))
      n_groups_total <- length(all_groups)
      groups_used <- all_groups[seq_len(min(max_groups, n_groups_total))]
      groups_dropped <- setdiff(all_groups, groups_used)

      df <- df |>
        dplyr::filter(cell_type %in% groups_used)

      rank_top_tfs <- function(direction_filter) {
        df |>
          dplyr::filter(direction == direction_filter) |>
          dplyr::group_by(TF) |>
          dplyr::summarise(
            best_fdr = min(BH_FDR, na.rm = TRUE),
            best_p = min(p_value, na.rm = TRUE),
            .groups = "drop"
          ) |>
          dplyr::arrange(best_fdr, best_p) |>
          dplyr::slice_head(n = top_n) |>
          dplyr::pull(TF)
      }

      top_pos <- rank_top_tfs("pos")
      top_neg <- rank_top_tfs("neg")
      top_tfs <- union(top_pos, top_neg)

      shiny::validate(
        shiny::need(
          length(top_tfs) > 0,
          "No TFs pass the selected q-value cutoff in either direction."
        )
      )

      n_available_pos <- dplyr::n_distinct(df$TF[df$direction == "pos"])
      n_available_neg <- dplyr::n_distinct(df$TF[df$direction == "neg"])
      n_shown_pos <- length(top_pos)
      n_shown_neg <- length(top_neg)

      plot_df <- df |>
        dplyr::filter(TF %in% top_tfs) |>
        dplyr::group_by(TF, cell_type, direction) |>
        dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          direction = factor(direction, levels = c("pos", "neg")),
          cell_type = factor(cell_type, levels = groups_used)
        )

      tf_order <- plot_df |>
        dplyr::group_by(TF) |>
        dplyr::summarise(
          max_sig = max(minus_log10_q),
          .groups = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(max_sig)) |>
        dplyr::pull(TF)

      plot_df$TF <- factor(
        plot_df$TF,
        levels = rev(tf_order)
      )

      subtitle_lines <- c(
        paste0(
          "Positive: showing ", n_shown_pos, " of ", n_available_pos,
          " eligible TFs | Negative: showing ", n_shown_neg, " of ",
          n_available_neg, " eligible TFs"
        ),
        if (length(groups_dropped) > 0) {
          paste0(
            n_groups_total,
            " metadata groups loaded, showing first ",
            length(groups_used),
            " (",
            paste(groups_used, collapse = ", "),
            "). Not shown: ",
            paste(groups_dropped, collapse = ", "),
            ". Re-run PASTAA with fewer groups selected to compare others."
          )
        } else {
          NULL
        }
      )

      ggplot(
        plot_df,
        aes(
          x = cell_type,
          y = TF,
          size = minus_log10_q,
          colour = direction,
          group = direction
        )
      ) +
        geom_point(
          alpha = 0.9,
          position = position_dodge(width = 0.3)
        ) +
        scale_colour_manual(
          values = c(
            "pos" = "#D55E00",
            "neg" = "#0072B2"
          )
        ) +
        labs(
          x = "Metadata group",
          y = "TF",
          size = "-log10(q-value / BH FDR)",
          colour = "Feature Importance direction",
          title = paste0(
            "Top enriched TFs by group, top ",
            top_n,
            " per direction (ranked across all shown groups), q-value cutoff ≤ ",
            q_cutoff
          ),
          subtitle = paste(
            subtitle_lines,
            collapse = "\n"
          )
        ) +
        theme_bw(base_size = 14) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          plot.subtitle = element_text(
            size = 10,
            colour = "grey40"
          )
        )
    }

  make_pastaa_abs_dotplot <- function(
      df,
      top_n,
      q_cutoff,
      max_groups = 10
    ) {
      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          paste0(
            "No TF enrichment results pass the selected q-value cutoff (≤ ",
            q_cutoff,
            ")."
          )
        )
      )

      df <- df |>
        dplyr::filter(direction == "abs") |>
        dplyr::mutate(cell_type = as.character(cell_type))

      shiny::validate(
        shiny::need(
          nrow(df) > 0,
          "No absolute-direction TF enrichment results are available."
        )
      )

      all_groups <- sort(unique(df$cell_type))
      n_groups_total <- length(all_groups)
      groups_used <- all_groups[seq_len(min(max_groups, n_groups_total))]
      groups_dropped <- setdiff(all_groups, groups_used)

      df <- df |>
        dplyr::filter(cell_type %in% groups_used)

      top_tfs <- df |>
        dplyr::group_by(TF) |>
        dplyr::summarise(
          best_fdr = min(BH_FDR, na.rm = TRUE),
          best_p = min(p_value, na.rm = TRUE),
          .groups = "drop"
        ) |>
        dplyr::arrange(best_fdr, best_p) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::pull(TF)

      shiny::validate(
        shiny::need(
          length(top_tfs) > 0,
          "No TFs pass the selected q-value cutoff."
        )
      )

      n_available <- dplyr::n_distinct(df$TF)
      n_shown <- length(top_tfs)

      plot_df <- df |>
        dplyr::filter(TF %in% top_tfs) |>
        dplyr::group_by(TF, cell_type) |>
        dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::mutate(
          cell_type = factor(cell_type, levels = groups_used)
        )

      tf_order <- plot_df |>
        dplyr::group_by(TF) |>
        dplyr::summarise(
          max_sig = max(minus_log10_q),
          .groups = "drop"
        ) |>
        dplyr::arrange(dplyr::desc(max_sig)) |>
        dplyr::pull(TF)

      plot_df$TF <- factor(
        plot_df$TF,
        levels = rev(tf_order)
      )

      subtitle_lines <- c(
        paste0(
          "Showing ",
          n_shown,
          " of ",
          n_available,
          " eligible TFs"
        ),
        if (length(groups_dropped) > 0) {
          paste0(
            n_groups_total,
            " metadata groups loaded, showing first ",
            length(groups_used),
            " (",
            paste(groups_used, collapse = ", "),
            "). Not shown: ",
            paste(groups_dropped, collapse = ", "),
            "."
          )
        } else {
          NULL
        }
      )

      ggplot(
        plot_df,
        aes(
          x = cell_type,
          y = TF,
          size = minus_log10_q
        )
      ) +
        geom_point(
          alpha = 0.9,
          colour = "#5B2C6F"
        ) +
        labs(
          x = "Metadata group",
          y = "TF",
          size = "-log10(q-value / BH FDR)",
          title = paste0(
            "Top enriched TFs by group, top ",
            top_n,
            " (ranked across all shown groups) -- Absolute Feature Importance, ",
            "q-value cutoff ≤ ",
            q_cutoff
          ),
          subtitle = paste(
            subtitle_lines,
            collapse = "\n"
          )
        ) +
        theme_bw(base_size = 14) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
          ),
          plot.subtitle = element_text(
            size = 10,
            colour = "grey40"
          )
        )
    }


  make_pastaa_dotplot <- function(df, top_n, q_cutoff) {
    shiny::validate(
      shiny::need(
        !is.null(df) && nrow(df) > 0,
        paste0(
          "No TF enrichment results pass the selected q-value cutoff (≤ ",
          q_cutoff,
          "). Increase the q-value cutoff to display the dotplot."
        )
      )
    )

    top_tfs <- df |>
      dplyr::group_by(TF) |>
      dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::arrange(BH_FDR, p_value) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::pull(TF)

    plot_df <- df |>
      dplyr::filter(
        TF %in% top_tfs,
        direction %in% c("pos", "neg", "abs")
      ) |>
      dplyr::group_by(TF, cell_type, direction) |>
      dplyr::arrange(BH_FDR, p_value, .by_group = TRUE) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        TF = factor(TF, levels = rev(top_tfs)),
        direction = factor(
          direction,
          levels = c("pos", "neg", "abs")
        )
      )

    shiny::validate(
      shiny::need(
        nrow(plot_df) > 0,
        if (identical(input$pastaa_direction, "abs")) {
          paste0(
            "The dotplot could not be created for absolute Feature Importance because ",
            "no valid TF enrichment results are available for the selected settings. ",
            "Absolute Feature Importance combines positive and negative Feature Importance ",
            "values by ranking regions according to their absolute magnitude."
          )
        } else {
          "No PASTAA results are available for the selected settings."
        }
      )
    )

    ggplot(plot_df, aes(
      x = cell_type,
      y = TF,
      size = minus_log10_q,
      colour = direction,
      group = direction
    )) +
      geom_point(
        alpha = 0.9,
        shape = 16,
        position = position_dodge(width = 0.22)
      ) +
      labs(
        x = "Cell type",
        y = "TF",
        size = "-log10(q-value / BH FDR)",
        colour = "Feature Importance direction",
        title = paste0(
          "TF enrichment across cell types, q-value cutoff ≤ ",
          q_cutoff
        )
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
  }

  pastaa_selected_direction <- reactive({
    res_list <- pastaa_results()

    if (length(res_list) == 0) {
      return(NULL)
    }

    selected <- input$pastaa_selected_result_tab

    if (
      is.null(selected) ||
      !nzchar(selected) ||
      !selected %in% names(res_list)
    ) {
      selected <- names(res_list)[1]
    }

    df <- res_list[[selected]]$result

    if (
      is.null(df) ||
      nrow(df) == 0 ||
      !"direction" %in% names(df)
    ) {
      return(NULL)
    }

    unique(as.character(df$direction))[1]
  })


  output$pastaa_top_tf_barplot_title <- renderUI({
    direction <- pastaa_selected_direction()

    if (is.null(direction) || length(direction) != 1L || !nzchar(direction)) {
      return(NULL)
    }

    label <- switch(
      direction,
      pos = "Positive Feature Importance",
      neg = "Negative Feature Importance",
      abs = "Absolute Feature Importance",
      "Feature Importance"
    )

    h4(paste0(label, " (current selection)"))
  })


  output$pastaa_top_tf_barplot_current <- renderPlot({

    direction <- pastaa_selected_direction()

    shiny::validate(
      shiny::need(
        !is.null(direction),
        "No PASTAA results yet. Select a gene/biological group and click Start PASTAA."
      )
    )

    df <- pastaa_selected_direction_pair_df()

    make_pastaa_top_tf_plot(
      df = df,
      top_n = input$pastaa_save_top_n,
      q_cutoff = input$pastaa_q_cutoff,
      direction = direction
    )
  })


  output$pastaa_abs_dotplot <- renderPlot({
    df <- pastaa_summary_df()

    shiny::validate(
      shiny::need(
        !is.null(df) &&
          nrow(df) > 0 &&
          any(df$direction == "abs", na.rm = TRUE),
        "Run TF enrichment with Absolute Feature Importance to display this plot."
      )
    )

    make_pastaa_abs_dotplot(
      df = df,
      top_n = input$pastaa_save_top_n,
      q_cutoff = input$pastaa_q_cutoff
    )
  })


  output$pastaa_top_tf_barplot_combined <- renderPlot({

    df <- pastaa_summary_df()

    make_pastaa_posneg_combined_plot(
      df = df,
      top_n = input$pastaa_save_top_n,
      q_cutoff = input$pastaa_q_cutoff
    )
  })

  observe({
    res_list <- pastaa_results()
    if (length(res_list) == 0) return(NULL)

    for (nm in names(res_list)) {
      local({
        key <- nm
        output[[paste0("pastaa_table_", key)]] <- renderDT({
          df <- pastaa_results()[[key]]$result |>
            dplyr::filter(TF != "DUMMY_MOTIF") |>
            dplyr::mutate(
              p_value = suppressWarnings(as.numeric(p_value)),
              BH_FDR = suppressWarnings(as.numeric(BH_FDR))
            ) |>
            dplyr::arrange(BH_FDR, p_value) |>
            dplyr::mutate(
              p_value = signif(p_value, 4),
              BH_FDR = signif(BH_FDR, 4)
            )

          datatable(
            df,
            options = list(
              pageLength = 20,
              scrollX = TRUE,
              displayStart = 0
            ),
            rownames = FALSE
          )
        })
      })
    }
  })

  observeEvent(input$run_fimo, {

    # ------------------------------------------------------------
    # Safety checks before running FIMO
    # ------------------------------------------------------------

    if (!isTRUE(igv_loaded())) {
      status_msg("Please load IGV first before running FIMO.")
      return(NULL)
    }

    if (
      !is.finite(input$fimo_top_n) ||
        input$fimo_top_n < 1 ||
        input$fimo_top_n > 20
    ) {
      status_msg("FIMO top N must be between 1 and 20.")
      return(NULL)
    }

    selected_backgrounds <- c("encode")

    gene_resolved <- resolve_igv_gene_input(input$igv_gene)

    if (is.null(gene_resolved)) {
      status_msg("Cannot run FIMO: selected gene was not found.")
      return(NULL)
    }

    gene <- gene_resolved$gene_id
    gene_label <- gene_resolved$gene_display

    if (!identical(gene, igv_current_gene())) {
      status_msg("Please reload IGV for the selected gene before running FIMO.")
      return(NULL)
    }

    # Only disable the button after all quick checks passed.
    # This prevents the button from briefly disabling for simple input errors.
    shinyjs::disable("run_fimo")
    on.exit(shinyjs::enable("run_fimo"), add = TRUE)

    # ------------------------------------------------------------
    # Prepare TFs
    # ------------------------------------------------------------

    status_msg("Preparing TFs for FIMO...")

    tf_prep <- tryCatch(
      get_fimo_motif_ids_from_input(
        top_n = input$fimo_top_n,
        manual_text = input$fimo_manual_tfs,
        source = input$fimo_tf_source,
        current_gene = gene
      ),
      error = function(e) e
    )

    if (inherits(tf_prep, "error")) {
      status_msg(paste("Cannot prepare TFs for FIMO:", tf_prep$message))
      return(NULL)
    }

    motif_ids_for_fimo <- tf_prep$motif_ids
    tf_mapping_table <- tf_prep$mapping_table
    tf_source_table <- tf_prep$source_table
    missing_tfs_for_fimo <- tf_prep$missing_tfs %||% character(0)
    motif_label_map_for_fimo <- tf_prep$motif_label_map %||% character(0)

    tf_names_selected <- if ("TF_clean" %in% colnames(tf_source_table)) {
      unique(tf_display_name(tf_source_table$TF_clean))
    } else {
      unique(tf_display_name(tf_source_table$TF))
    }

    tf_names_selected <- tf_names_selected[
      !is.na(tf_names_selected) & nzchar(tf_names_selected)
    ]

    status_msg(
      paste0(
        "Prepared ",
        length(tf_names_selected),
        " TF(s), mapped to ",
        length(motif_ids_for_fimo),
        " motif ID(s)."
      )
    )

    # ------------------------------------------------------------
    # Decide which PASTAA combinations / regions should be used
    # ------------------------------------------------------------

    if (input$fimo_tf_source == "saved") {

      combos <- tryCatch(
        get_latest_saved_pastaa_combinations(current_gene = gene),
        error = function(e) e
      )

      if (inherits(combos, "error")) {
        status_msg(paste("Cannot prepare PASTAA combinations for FIMO:", combos$message))
        return(NULL)
      }

      if (nrow(combos) == 0) {
        status_msg("No saved PASTAA combinations found for current IGV gene.")
        return(NULL)
      }

    } else {

      selected_cts <- input$igv_celltypes

      if (is.null(selected_cts) || length(selected_cts) == 0) {
        selected_cts <- if (length(igv_cell_types) > 0) sort(igv_cell_types)[1] else character(0)
      }

      selected_cts <- intersect(selected_cts, igv_cell_types)

      if (length(selected_cts) == 0) {
        status_msg("Cannot run FIMO: no valid IGV cell types selected.")
        return(NULL)
      }

      combos <- data.frame(
        gene_id = gene,
        cell_type = selected_cts,
        stringsAsFactors = FALSE
      )
    }

    if (is.null(combos) || nrow(combos) == 0) {
      status_msg("Cannot run FIMO: no gene/cell-type/direction combinations available.")
      return(NULL)
    }

    status_msg(
      paste0(
        "Running FIMO for ",
        gene_label,
        " using ",
        length(tf_names_selected),
        " selected TF(s), ",
        length(motif_ids_for_fimo),
        " mapped motif ID(s), ",
        nrow(combos),
        " region set(s), and ENCODE background..."
      )
    )

    fimo_bed_files_by_bg <- setNames(
      vector("list", length(selected_backgrounds)),
      selected_backgrounds
    )

    failed <- character(0)
    total_sites <- 0

    # ------------------------------------------------------------
    # Run FIMO for selected backgrounds
    # ------------------------------------------------------------

    withProgress(message = "Running FIMO", value = 0, {

      total <- nrow(combos) * length(selected_backgrounds)

      if (!is.finite(total) || total <= 0) {
        stop("No valid FIMO runs available.")
      }

      for (bg_type in selected_backgrounds) {

        bg_file <- tryCatch(
          ensure_fimo_bg_file(bg_type),
          error = function(e) e
        )

        if (inherits(bg_file, "error")) {
          failed <<- c(failed, paste0(bg_type, ": ", bg_file$message))
          next
        }

        for (i in seq_len(nrow(combos))) {

          incProgress(1 / total)

          combo_gene <- as.character(combos$gene_id[i])
          ct <- as.character(combos$cell_type[i])
          

          run_gene_for_file <- if (combo_gene == "ALL_GENES") {
            "ALL_GENES"
          } else {
            gene
          }

          run_id <- paste(
            make_safe_id(run_gene_for_file),
            make_safe_id(ct),
            "abs",
            sep = "_"
          )

          message("Running FIMO: ", run_id, " | background=", bg_type)

          run_mode_for_fimo <- if (combo_gene == "ALL_GENES") {
            "celltype"
          } else {
            "gene"
          }

          region_direction_for_fimo <- "abs"

          res <- tryCatch({

            ensure_fimo_fasta(
              gene = gene,
              ct = ct,
              direction = region_direction_for_fimo,
              run_id = run_id,
              n_regions = input$pastaa_n_regions,
              run_mode = run_mode_for_fimo,
              use_abs_shap = TRUE
            )

            run_fimo_for_pastaa(
              run_id = run_id,
              top_n = input$fimo_top_n,
              motif_ids_override = motif_ids_for_fimo,
              motif_label_map = motif_label_map_for_fimo,
              bg_type = bg_type,
              bg_file = bg_file
            )

          }, error = function(e) {
            failed <<- c(
              failed,
              paste0(run_id, " / ", bg_type, ": ", e$message)
            )
            NULL
          })

          if (!is.null(res) && !is.null(res$bed_file) && file.exists(res$bed_file)) {
            fimo_bed_files_by_bg[[bg_type]] <- c(
              fimo_bed_files_by_bg[[bg_type]],
              res$bed_file
            )

            total_sites <- total_sites + res$n_sites
          }
        }
      }
    })

    # ------------------------------------------------------------
    # Stop if no FIMO sites were created
    # ------------------------------------------------------------

    if (sum(lengths(fimo_bed_files_by_bg)) == 0) {

      if (length(failed) > 0) {
        status_msg(
          paste(
            "FIMO finished, but no binding sites were added:",
            paste(failed, collapse = " | ")
          )
        )
      } else {
        status_msg("FIMO finished, but no binding sites were found.")
      }

      return(NULL)
    }

    # ------------------------------------------------------------
    # Merge BED files and create one IGV track per TF/motif
    # ------------------------------------------------------------

    cache_bust <- as.integer(Sys.time())
    fimo_tracks <- list()
    sites_by_bg <- c()
    sites_by_tf <- c()
    fimo_export_rows <- list()

    for (bg_type in names(fimo_bed_files_by_bg)) {

      bg_files <- fimo_bed_files_by_bg[[bg_type]]

      if (length(bg_files) == 0) next

      all_bed <- dplyr::bind_rows(lapply(bg_files, function(f) {
        read.table(
          f,
          sep = "\t",
          stringsAsFactors = FALSE,
          quote = "",
          col.names = c("chrom", "start", "end", "name", "score", "strand")
        )
      }))

      all_bed <- all_bed |>
        dplyr::mutate(
          chrom = as.character(chrom),
          start = suppressWarnings(as.integer(start)),
          end = suppressWarnings(as.integer(end)),
          name = as.character(name),
          score = suppressWarnings(as.integer(score)),
          strand = as.character(strand),

          # name looks like: TFNAME|+ or TFNAME|-
          tf_label = sub("\\|.*$", "", name),
          tf_label = ifelse(
            !is.na(tf_label) & nzchar(tf_label),
            tf_label,
            name
          )
        ) |>
        dplyr::filter(
          !is.na(chrom),
          !is.na(start),
          !is.na(end),
          end > start,
          !is.na(name),
          nzchar(name),
          !is.na(tf_label),
          nzchar(tf_label)
        ) |>
        dplyr::distinct(
          chrom,
          start,
          end,
          name,
          score,
          strand,
          tf_label,
          .keep_all = TRUE
        ) |>
        dplyr::arrange(chrom, start, end)

      if (nrow(all_bed) == 0) next

      sites_by_bg[bg_type] <- nrow(all_bed)

      tf_labels <- sort(unique(all_bed$tf_label))

      if (length(tf_labels) > 20) {
        tf_labels <- tf_labels[seq_len(20)]
      }

      track_color <- "#7B1FA2"

      for (tf in tf_labels) {
        tf_bed <- all_bed |>
          dplyr::filter(tf_label == tf) |>
          dplyr::select(
            chrom,
            start,
            end,
            name,
            score,
            strand
          ) |>
          dplyr::arrange(chrom, start, end)

        if (nrow(tf_bed) == 0) next
        fimo_export_rows[[length(fimo_export_rows) + 1]] <- tf_bed |>
        dplyr::mutate(
          track = paste0("FIMO - ", tf),
          track_type = "FIMO",
          fimo_label = gsub("\\|", " ", name)
        )

        tf_safe <- make_safe_id(tf)

        merged_bed_name <- paste0(
          gene,
          "_fimo_",
          tf_safe,
          ".bed"
        )

        merged_bed_file <- file.path(fimo_www_dir, merged_bed_name)

        write.table(
          tf_bed,
          merged_bed_file,
          sep = "\t",
          quote = FALSE,
          row.names = FALSE,
          col.names = FALSE
        )

        sites_by_tf[tf] <- nrow(tf_bed)

        fimo_tracks[[length(fimo_tracks) + 1]] <- list(
          id = paste0(
            "fimo__",
            tf_safe,
            "__",
            cache_bust
          ),
          name = paste0("FIMO - ", tf),
          type = "annotation",
          format = "bed",
          url = paste0("fimo/", merged_bed_name, "?v=", cache_bust),

          displayMode = "SQUISHED",
          height = 40,

          color = track_color,
          labelColor = track_color,

          showDataRange = FALSE,
          autoscale = FALSE,

          visibilityWindow = 1e9,
          order = 1000L + length(fimo_tracks)
        )
      }
    }
    if (length(fimo_export_rows) > 0) {
      current_fimo_export_df(dplyr::bind_rows(fimo_export_rows))
    } else {
      current_fimo_export_df(NULL)
    }

    if (length(fimo_tracks) == 0) {
      status_msg("FIMO finished, but no valid merged BED tracks were created.")
      return(NULL)
    }

    # ------------------------------------------------------------
    # Add FIMO tracks to IGV
    # ------------------------------------------------------------

    if (isTRUE(input$fimo_replace_tracks)) {
      session$sendCustomMessage("igv-remove-fimo-tracks", list())
    }

    session$sendCustomMessage(
      "igv-add-tracks",
      list(tracks = fimo_tracks)
    )

    # ------------------------------------------------------------
    # Status message
    # ------------------------------------------------------------

    tf_names_for_status <- unique(tf_display_name(tf_names_selected))
    tracks_with_hits <- unique(names(sites_by_tf))

    missing_motif_display <- unique(tf_display_name(missing_tfs_for_fimo))
    missing_motif_display <- missing_motif_display[
      !is.na(missing_motif_display) & nzchar(missing_motif_display)
    ]

    no_hit_tfs <- setdiff(tf_names_for_status, tracks_with_hits)
    no_hit_tfs <- setdiff(no_hit_tfs, missing_motif_display)

    no_track_reasons <- c()

    if (length(missing_motif_display) > 0) {
      no_track_reasons <- c(
        no_track_reasons,
        paste0(
          missing_motif_display,
          ": no matching JASPAR motif found"
        )
      )
    }

    if (length(no_hit_tfs) > 0) {
      no_track_reasons <- c(
        no_track_reasons,
        paste0(
          no_hit_tfs,
          ": no FIMO hits in selected regions"
        )
      )
    }

    if (length(failed) > 0) {
      no_track_reasons <- c(
        no_track_reasons,
        paste0(
          "Skipped FIMO run(s): ",
          paste(failed, collapse = " | ")
        )
      )
    }

    status_msg(
      paste0(
        "FIMO tracks added: ",
        length(fimo_tracks),
        " track(s)",
        "\nGene: ",
        gene_label,
        "\nSelected TFs from PASTAA/manual input: ",
        paste(tf_names_for_status, collapse = ", "),
        "\nTracks with FIMO hits: ",
        if (length(tracks_with_hits) > 0) {
          paste(tracks_with_hits, collapse = ", ")
        } else {
          "none"
        },
        "\nENCODE background sites: ",
        paste(sites_by_bg, collapse = ", "),
        "\nTFs without track / reason: ",
        if (length(no_track_reasons) > 0) {
          paste(no_track_reasons, collapse = " | ")
        } else {
          "none"
        }
      )
    )
  })
  # ======================================================================
  # DOWNLOAD HANDLERS
  # ======================================================================
  output$download_pastaa_top_tf_pos_pdf <- downloadHandler(
    filename = function() {
      gene_label <- input$pastaa_gene

      if (is.null(gene_label) || !nzchar(gene_label)) {
        gene_label <- "PASTAA"
      }

      paste0(
        make_safe_id(gene_label),
        "_TF_enrichment_positive_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      df <- pastaa_selected_direction_pair_df()

      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          "No PASTAA results available."
        )
      )

      p <- make_pastaa_top_tf_plot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff,
        direction = "pos"
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 8.5,
        height = 5.5,
        onefile = TRUE
      )

      print(p)

      grDevices::dev.off()
    }
  )
  pastaa_filename_label <- function(df) {
    if (is.null(df) || nrow(df) == 0 || !"run_mode" %in% names(df)) {
      return("PASTAA")
    }

    modes <- unique(df$run_mode)

    if (length(modes) == 1L && identical(modes, "celltype")) {
      groups <- unique(as.character(df$cell_type))
      paste(make_safe_id(groups), collapse = "_vs_")

    } else if (length(modes) == 1L && identical(modes, "gene")) {
      genes <- unique(as.character(df$gene_id))
      make_safe_id(paste(genes, collapse = "_"))

    } else {
      "PASTAA"
    }
  }
  output$download_pastaa_top_tf_current_pdf <- downloadHandler(
    filename = function() {
      label <- pastaa_filename_label(
        pastaa_selected_direction_pair_df()
      )

      direction <- pastaa_selected_direction()

      direction_label <- switch(
        direction,
        pos = "positive",
        neg = "negative",
        abs = "absolute",
        "unknown"
      )

      paste0(
        label,
        "_TF_enrichment_",
        direction_label,
        "_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      direction <- pastaa_selected_direction()

      shiny::validate(
        shiny::need(
          !is.null(direction),
          "No PASTAA results available."
        )
      )

      df <- pastaa_selected_direction_pair_df()

      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          "No PASTAA results available."
        )
      )

      p <- make_pastaa_top_tf_plot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff,
        direction = direction
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 8.5,
        height = 5.5,
        onefile = TRUE
      )

      print(p)
      grDevices::dev.off()
    }
  )


  output$download_pastaa_abs_dotplot_pdf <- downloadHandler(
    filename = function() {
      df <- pastaa_summary_df()

      label <- if (
        is.null(df) ||
        nrow(df) == 0 ||
        !"cell_type" %in% names(df)
      ) {
        "PASTAA"

      } else {
        groups <- unique(
          as.character(
            df$cell_type[df$direction == "abs"]
          )
        )

        safe_groups <- make_safe_id(groups)

        if (length(safe_groups) == 0L) {
          "PASTAA"

        } else if (length(safe_groups) <= 3L) {
          paste(
            safe_groups,
            collapse = "_vs_"
          )

        } else {
          paste0(
            paste(
              safe_groups[1:3],
              collapse = "_vs_"
            ),
            "_and_",
            length(safe_groups) - 3L,
            "more"
          )
        }
      }

      paste0(
        label,
        "_TF_enrichment_absolute_dotplot_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      df <- pastaa_summary_df()

      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          "No PASTAA results available."
        )
      )

      p <- make_pastaa_abs_dotplot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 9.5,
        height = 6.5,
        onefile = TRUE
      )

      print(p)
      grDevices::dev.off()
    }
  )


  output$download_pastaa_top_tf_combined_pdf <- downloadHandler(
    filename = function() {
      df <- pastaa_summary_df()

      label <- if (
        is.null(df) ||
        nrow(df) == 0 ||
        !"cell_type" %in% names(df)
      ) {
        "PASTAA"

      } else {
        groups <- unique(
          as.character(df$cell_type)
        )

        safe_groups <- make_safe_id(groups)

        if (length(safe_groups) <= 3L) {
          paste(
            safe_groups,
            collapse = "_vs_"
          )

        } else {
          paste0(
            paste(
              safe_groups[1:3],
              collapse = "_vs_"
            ),
            "_and_",
            length(safe_groups) - 3L,
            "more"
          )
        }
      }

      paste0(
        label,
        "_TF_enrichment_combined_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      df <- pastaa_summary_df()

      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          "No PASTAA results available."
        )
      )

      p <- make_pastaa_posneg_combined_plot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 9.5,
        height = 6.5,
        onefile = TRUE
      )

      print(p)
      grDevices::dev.off()
    }
  )

  output$download_pastaa_top_tf_neg_pdf <- downloadHandler(
    filename = function() {
      gene_label <- input$pastaa_gene

      if (is.null(gene_label) || !nzchar(gene_label)) {
        gene_label <- "PASTAA"
      }

      paste0(
        make_safe_id(gene_label),
        "_TF_enrichment_negative_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },

    content = function(file) {
      df <- pastaa_selected_direction_pair_df()

      shiny::validate(
        shiny::need(
          !is.null(df) && nrow(df) > 0,
          "No PASTAA results available."
        )
      )

      p <- make_pastaa_top_tf_plot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff,
        direction = "neg"
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 8.5,
        height = 5.5,
        onefile = TRUE
      )

      print(p)

      grDevices::dev.off()
    }
  )


  output$download_pastaa_dotplot_pdf <- downloadHandler(
    filename = function() {
      gene_label <- input$pastaa_gene
      if (is.null(gene_label) || !nzchar(gene_label)) {
        gene_label <- "PASTAA"
      }

      paste0(
        make_safe_id(gene_label),
        "_TF_enrichment_dotplot_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".pdf"
      )
    },
    content = function(file) {
      df <- pastaa_summary_df()

      shiny::validate(
        shiny::need(!is.null(df) && nrow(df) > 0, "No PASTAA results available.")
      )

      p <- make_pastaa_dotplot(
        df = df,
        top_n = input$pastaa_save_top_n,
        q_cutoff = input$pastaa_q_cutoff
      )

      grDevices::cairo_pdf(
        filename = file,
        width = 9.5,
        height = 6.5,
        onefile = TRUE
      )

      print(p)

      grDevices::dev.off()
    }
  )

  make_table_downloads <- function(
    id_prefix,
    all_reactive,
    file_prefix,
    gene_reactive = NULL,
    ct_reactive = NULL
  ) {
    output[[paste0(id_prefix, "_all_csv")]] <- downloadHandler(
      filename = function() {
        make_download_name(
          paste0(file_prefix, "_all"),
          "csv",
          gene = if (!is.null(gene_reactive)) gene_reactive() else NULL,
          ct = if (!is.null(ct_reactive)) ct_reactive() else NULL
        )
      },
      content = function(file) {
        write_download_table(all_reactive(), file, sep = ",")
      }
    )

    output[[paste0(id_prefix, "_all_tsv")]] <- downloadHandler(
      filename = function() {
        make_download_name(
          paste0(file_prefix, "_all"),
          "tsv",
          gene = if (!is.null(gene_reactive)) gene_reactive() else NULL,
          ct = if (!is.null(ct_reactive)) ct_reactive() else NULL
        )
      },
      content = function(file) {
        write_download_table(all_reactive(), file, sep = "\t")
      }
    )
  }

  make_table_downloads(
    id_prefix = "overview",
    all_reactive = overview_table_all,
    file_prefix = "performance"
  )

  make_table_downloads(
    id_prefix = "train_test",
    all_reactive = train_test_table_all,
    file_prefix = "train_test_genes"
  )

  make_table_downloads(
    id_prefix = "shap",
    all_reactive = shap_all_current,
    file_prefix = "shap_regions",
    gene_reactive = reactive(input$shap_gene),
    ct_reactive = reactive(input$shap_cell_type)
  )

  
  output$pastaa_all_csv <- downloadHandler(

    filename = function() {

      run_mode <- input$pastaa_run_mode %||% "gene"
      n_regions <- as.integer(input$pastaa_n_regions)

      if (identical(run_mode, "gene")) {

        gene_resolved <- resolve_igv_gene_input(input$pastaa_gene)

        gene_id <- if (!is.null(gene_resolved)) {
          gene_resolved$gene_id
        } else {
          make_safe_id(input$pastaa_gene %||% "gene")
        }

        paste0(
          "pastaa_",
          make_safe_id(gene_id),
          "_n",
          n_regions,
          ".csv"
        )

      } else {

        groups <- input$pastaa_celltypes %||% character(0)

        group_label <- if (length(groups) == 1) {
          make_safe_id(groups)
        } else {
          paste0(length(groups), "_biological_groups")
        }

        parts <- c(
          "pastaa",
          group_label,
          paste0("n", n_regions)
        )

        if (isTRUE(input$pastaa_filter_performance)) {

          threshold <- format(
            as.numeric(input$pastaa_min_test_correlation),
            scientific = FALSE,
            trim = TRUE
          )

          threshold <- gsub("\\.", "p", threshold)

          parts <- c(
            parts,
            paste0("corr", threshold)
          )
        }

        paste0(
          paste(parts, collapse = "_"),
          ".csv"
        )
      }
    },

    content = function(file) {
      write_download_table(
        pastaa_all_table(),
        file,
        sep = ","
      )
    }
  )


  output$pastaa_all_tsv <- downloadHandler(

    filename = function() {

      run_mode <- input$pastaa_run_mode %||% "gene"
      n_regions <- as.integer(input$pastaa_n_regions)

      if (identical(run_mode, "gene")) {

        gene_resolved <- resolve_igv_gene_input(input$pastaa_gene)

        gene_id <- if (!is.null(gene_resolved)) {
          gene_resolved$gene_id
        } else {
          make_safe_id(input$pastaa_gene %||% "gene")
        }

        paste0(
          "pastaa_",
          make_safe_id(gene_id),
          "_n",
          n_regions,
          ".tsv"
        )

      } else {

        groups <- input$pastaa_celltypes %||% character(0)

        group_label <- if (length(groups) == 1) {
          make_safe_id(groups)
        } else {
          paste0(length(groups), "_biological_groups")
        }

        parts <- c(
          "pastaa",
          group_label,
          paste0("n", n_regions)
        )

        if (isTRUE(input$pastaa_filter_performance)) {

          threshold <- format(
            as.numeric(input$pastaa_min_test_correlation),
            scientific = FALSE,
            trim = TRUE
          )

          threshold <- gsub("\\.", "p", threshold)

          parts <- c(
            parts,
            paste0("corr", threshold)
          )
        }

        paste0(
          paste(parts, collapse = "_"),
          ".tsv"
        )
      }
    },

    content = function(file) {
      write_download_table(
        pastaa_all_table(),
        file,
        sep = "\t"
      )
    }
  )

}

shinyApp(ui, server)