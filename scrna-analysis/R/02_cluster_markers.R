cluster_markers <- function(obj, cfg, out) {
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$analysis$variable_features)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(obj))
  npcs <- min(cfg$analysis$pcs, ncol(obj) - 1L, length(Seurat::VariableFeatures(obj)) - 1L)
  assert(npcs >= 2, "Insufficient variable features for PCA.")
  obj <- Seurat::RunPCA(obj, npcs = npcs, seed.use = cfg$seed)
  obj <- Seurat::FindNeighbors(obj, dims = seq_len(npcs), k.param = min(20L, ncol(obj) - 1L))
  obj <- Seurat::FindClusters(obj, resolution = cfg$analysis$resolution, random.seed = cfg$seed)
  obj <- Seurat::RunUMAP(obj, dims = seq_len(npcs), n.neighbors = min(30L, ncol(obj) - 1L), seed.use = cfg$seed)
  markers <- Seurat::FindAllMarkers(obj, assay = "RNA", slot = "data", test.use = "wilcox",
                                   only.pos = TRUE, min.pct = cfg$analysis$min_pct,
                                   logfc.threshold = cfg$analysis$logfc_threshold)
  write_table(markers, file.path(out, "markers.csv"))
  plot_save(Seurat::DimPlot(obj, group.by = "seurat_clusters", label = TRUE), file.path(out, "clusters.pdf"))
  plot_save(Seurat::DimPlot(obj, group.by = "sample_id"), file.path(out, "samples.pdf"))
  plot_save(Seurat::ElbowPlot(obj, ndims = npcs), file.path(out, "pca_elbow.pdf"))
  write_table(data.frame(cluster = sort(unique(as.character(obj$seurat_clusters))), cell_type = ""),
              file.path(out, "annotation_template.csv"))
  saveRDS(obj, file.path(out, "clustered.rds"))
  obj
}
