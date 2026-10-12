# Prepare the existing chromosome 21 EpIVoi object
# without modifying the original RDS.

data_dir <- "epivoi_object/data"

original <- file.path(
  data_dir,
  "epivoi_obj_cre_rf_IHEC_1MB_test_performance_chr21.rds"
)

output <- file.path(data_dir, "epivoi_chr21_local.rds")

if (!file.exists(original)) {
  stop("Original chromosome 21 RDS not found: ", original)
}

obj <- readRDS(original)

for (slot in c("feature_importance", "epigenetic_signal")) {
  idx <- obj[[slot]]$file_index

  idx$file_path <- file.path(
    slot,
    basename(idx$file_path)
  )

  missing <- idx$file_path[
    !file.exists(file.path(data_dir, idx$file_path))
  ]

  if (length(missing) > 0) {
    stop(
      "Missing files in ", slot, ": ",
      paste(head(missing, 5), collapse = ", ")
    )
  }

  obj[[slot]]$file_index <- idx

  cat(slot, ":", nrow(idx), "files OK\n")
}

saveRDS(obj, output)

cat("\nSUCCESS: local chromosome 21 object prepared\n")
cat("Output:", output, "\n")
