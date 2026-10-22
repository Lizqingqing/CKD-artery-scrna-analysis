annotate_cells <- function(obj, cfg, out) {
  map <- utils::read.csv(cfg$annotation_file, colClasses = "character", stringsAsFactors = FALSE)
  assert(all(c("cluster", "cell_type") %in% names(map)), "Annotation needs cluster and cell_type columns.")
  assert(!anyNA(map) && !anyDuplicated(map$cluster), "Missing or duplicated annotation entries.")
  clusters <- as.character(obj$seurat_clusters)
  assert(all(map$cluster %in% clusters), "Annotation contains unknown clusters.")
  assert(all(nzchar(trimws(map$cell_type))), "Annotation labels cannot be blank.")
  labels <- map$cell_type[match(clusters, map$cluster)]
  labels[is.na(labels)] <- "Unassigned"
  obj$cell_type <- labels
  write_table(data.frame(cell_id = colnames(obj), obj[[]]), file.path(out, "cell_metadata.csv"))
  plot_save(Seurat::DimPlot(obj, group.by = "cell_type", label = TRUE), file.path(out, "annotation.pdf"))
  saveRDS(obj, file.path(out, "annotated.rds"))
  obj
}
