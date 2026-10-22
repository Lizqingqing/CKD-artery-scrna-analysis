run_trajectory <- function(obj, cfg, out) {
  need(c("monocle3", "SingleCellExperiment"))
  tr <- cfg$trajectory
  md <- obj[[]]
  assert(length(tr$subset_values) > 0 && tr$subset_column %in% names(md),
         "Select a biologically coherent lineage with trajectory subset_column and subset_values.")
  cells <- rownames(md)[md[[tr$subset_column]] %in% unlist(tr$subset_values)]
  assert(length(cells) >= 30, "Trajectory needs at least 30 selected cells.")
  x <- counts_of(obj)[, cells, drop = FALSE]
  x <- x[Matrix::rowSums(x) > 0, , drop = FALSE]
  md <- md[colnames(x), , drop = FALSE]
  assert(tr$root_column %in% names(md) && length(tr$root_values) > 0,
         "Specify an evidence-based root_column and root_values; roots are never guessed.")
  roots <- rownames(md)[md[[tr$root_column]] %in% unlist(tr$root_values)]
  assert(length(roots) > 0, "No selected cells match the root criteria.")
  genes <- data.frame(gene_short_name = rownames(x), row.names = rownames(x))
  cds <- monocle3::new_cell_data_set(x, cell_metadata = md, gene_metadata = genes)
  nd <- min(tr$num_dim, nrow(x) - 1L, ncol(x) - 1L)
  assert(nd >= 2, "Too few genes or cells for trajectory PCA.")
  set.seed(cfg$seed)
  cds <- monocle3::preprocess_cds(cds, num_dim = nd)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", cores = 1,
                                    umap.fast_sgd = FALSE, umap.n_neighbors = min(15L, ncol(cds) - 1L))
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", random_seed = cfg$seed,
                               k = min(20L, ncol(cds) - 1L))
  cds <- monocle3::learn_graph(cds, use_partition = tr$use_partition, close_loop = FALSE)
  saveRDS(cds, file.path(out, "trajectory_unordered.rds"))
  parts <- monocle3::partitions(cds)
  write_table(data.frame(cell_id = names(parts), partition = as.character(parts)),
              file.path(out, "trajectory_partitions.csv"))
  if (isTRUE(tr$use_partition)) {
    assert(all(unique(parts) %in% parts[roots]),
           "At least one partition lacks a root. Review lineage selection or supply justified roots for each partition.")
  }
  cds <- monocle3::order_cells(cds, root_cells = roots)
  pt <- monocle3::pseudotime(cds)
  write_table(data.frame(cell_id = names(pt), pseudotime = as.numeric(pt), finite = is.finite(pt)),
              file.path(out, "pseudotime.csv"))
  assert(all(is.finite(pt)), "Unreachable cells have infinite pseudotime. Inspect graph connectivity and roots.")
  saveRDS(cds, file.path(out, "trajectory.rds"))
  plot_save(monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE),
             file.path(out, "trajectory_pseudotime.pdf"))
  if (isTRUE(tr$graph_test)) {
    test <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
    write_table(data.frame(gene = rownames(test), test), file.path(out, "trajectory_graph_test.csv"))
  }
  invisible(cds)
}
