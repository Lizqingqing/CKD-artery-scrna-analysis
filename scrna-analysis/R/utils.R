assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
need <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(!length(missing), paste("Install required packages:", paste(missing, collapse = ", ")))
}
write_table <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}
counts_of <- function(obj) {
  # JoinLayers must be called after merging and before extracting raw counts.
  SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
}
validate_counts <- function(x) {
  assert(inherits(x, "sparseMatrix"), "Expected a sparse gene-by-cell count matrix.")
  assert(nrow(x) > 2 && ncol(x) > 2, "The count matrix is too small.")
  assert(!anyDuplicated(rownames(x)) && !anyDuplicated(colnames(x)), "Duplicate matrix names.")
  assert(all(is.finite(x@x)) && all(x@x >= 0) && all(x@x == floor(x@x)),
         "Input must contain finite, nonnegative integer raw counts.")
  assert(all(Matrix::colSums(x) > 0), "Input contains empty cells; use a filtered 10x matrix.")
}
validate_config <- function(cfg) {
  assert(cfg$species %in% c("human", "mouse"), "species must be human or mouse.")
  q <- cfg$qc
  assert(q$min_features >= 0 && q$max_features > q$min_features && q$max_counts > 0 &&
           q$max_percent_mt >= 0 && q$max_percent_mt <= 100, "Invalid QC thresholds.")
  assert(cfg$analysis$pcs >= 2 && cfg$analysis$variable_features >= 3 &&
           cfg$analysis$resolution > 0, "Invalid analysis parameters.")
  assert(cfg$metabolism$max_dense_gb > 0, "max_dense_gb must be positive.")
}
plot_save <- function(plot, path, width = 8, height = 6) {
  ggplot2::ggsave(path, plot, width = width, height = height, limitsize = TRUE)
}
