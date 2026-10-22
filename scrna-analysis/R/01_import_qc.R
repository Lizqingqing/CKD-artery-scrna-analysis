import_qc <- function(cfg, out) {
  samples <- utils::read.csv(cfg$samples_file, stringsAsFactors = FALSE, colClasses = "character")
  required <- c("sample_id", "condition", "donor_id", "matrix_dir")
  assert(all(required %in% names(samples)) && nrow(samples) > 0, "Invalid sample manifest.")
  assert(!anyNA(samples[required]) && all(nzchar(as.matrix(samples[required]))), "Missing sample metadata.")
  assert(!anyDuplicated(samples$sample_id) && all(grepl("^[A-Za-z0-9-]+$", samples$sample_id)),
         "Sample IDs must be unique and contain only letters, numbers, and hyphens.")
  objects <- vector("list", nrow(samples))
  for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    assert(dir.exists(s$matrix_dir), paste("Missing 10x directory:", s$matrix_dir))
    x <- Seurat::Read10X(s$matrix_dir, gene.column = 2, unique.features = TRUE)
    if (is.list(x)) {
      assert("Gene Expression" %in% names(x), "No Gene Expression matrix in multimodal input.")
      x <- x[["Gene Expression"]]
    }
    validate_counts(x)
    obj <- Seurat::CreateSeuratObject(x, project = s$sample_id, min.cells = 0, min.features = 0)
    obj <- Seurat::RenameCells(obj, add.cell.id = s$sample_id)
    obj$sample_id <- s$sample_id
    obj$condition <- s$condition
    obj$donor_id <- s$donor_id
    mt <- grep(cfg$qc$mitochondrial_pattern, rownames(obj), value = TRUE)
    assert(length(mt) > 0, "No mitochondrial genes matched. Check species, symbols, and regex.")
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, features = mt)
    objects[[i]] <- obj
  }
  obj <- if (length(objects) == 1) objects[[1]] else merge(objects[[1]], y = objects[-1])
  obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  md <- obj[[]]
  q <- cfg$qc
  md$qc_pass <- with(md, nFeature_RNA >= q$min_features & nFeature_RNA <= q$max_features &
                      nCount_RNA <= q$max_counts & percent.mt <= q$max_percent_mt)
  write_table(data.frame(cell_id = rownames(md), md), file.path(out, "qc_all_cells.csv"))
  write_table(as.data.frame(table(sample_id = md$sample_id, retained = md$qc_pass)),
              file.path(out, "qc_summary.csv"))
  plot_save(Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          group.by = "sample_id", pt.size = 0, ncol = 3), file.path(out, "qc_before.pdf"), 12, 5)
  assert(sum(md$qc_pass) >= 20, "Fewer than 20 cells passed QC.")
  assert(all(samples$sample_id %in% md$sample_id[md$qc_pass]), "A sample lost all cells during QC.")
  obj <- subset(obj, cells = rownames(md)[md$qc_pass])
  plot_save(Seurat::VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                          group.by = "sample_id", pt.size = 0, ncol = 3), file.path(out, "qc_after.pdf"), 12, 5)
  saveRDS(obj, file.path(out, "qc.rds"))
  obj
}
