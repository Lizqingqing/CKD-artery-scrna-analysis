#!/usr/bin/env Rscript
# Run from the repository root. Each stage writes to a new immutable run directory.
source("R/utils.R")
need(c("yaml", "Seurat", "SeuratObject", "Matrix", "ggplot2"))
assert(utils::packageVersion("Seurat") >= "5.0.0", "Seurat >= 5.0.0 is required.")
args <- commandArgs(trailingOnly = TRUE)
assert(length(args) >= 2, "Usage: Rscript scripts/run_pipeline.R config.yml stage [input.rds]")
cfg <- yaml::read_yaml(args[1])
validate_config(cfg)
stage <- args[2]
assert(stage %in% c("standard", "annotate", "metabolism", "trajectory"), "Unknown stage.")
for (file in list.files("R", pattern = "^[0-9].*\\.R$", full.names = TRUE)) source(file)
out <- file.path(cfg$output_dir, stage)
assert(!dir.exists(out), paste("Output already exists; choose a new output_dir:", out))
dir.create(out, recursive = TRUE)
file.copy(args[1], file.path(out, "config.yml"))
for (f in c(cfg$samples_file, cfg$annotation_file)) {
  if (file.exists(f)) file.copy(f, file.path(out, basename(f)))
}
set.seed(cfg$seed)
options(Seurat.object.assay.version = "v5")
status <- "FAILED"
tryCatch({
  if (stage == "standard") {
    obj <- import_qc(cfg, out)
    obj <- cluster_markers(obj, cfg, out)
  } else {
    assert(length(args) == 3 && file.exists(args[3]), "Supply an existing input RDS file.")
    writeLines(paste(args[3], tools::md5sum(args[3])), file.path(out, "input_checksum.txt"))
    obj <- readRDS(args[3])
    assert(inherits(obj, "Seurat"), "Input must be a Seurat object.")
    obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
    if (stage == "annotate") annotate_cells(obj, cfg, out)
    if (stage == "metabolism") {
      assert("cell_type" %in% names(obj[[]]), "Run annotation first.")
      score_metabolism(obj, cfg, out)
    }
    if (stage == "trajectory") run_trajectory(obj, cfg, out)
  }
  status <- "SUCCESS"
}, finally = {
  writeLines(status, file.path(out, "STATUS.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  packages <- as.data.frame(utils::installed.packages()[, c("Package", "Version")])
  write_table(packages, file.path(out, "package_versions.csv"))
  remote_info <- lapply(intersect(c("Seurat", "scMetabolism", "monocle3", "AUCell", "VISION"),
                                  rownames(utils::installed.packages())), function(pkg) {
    d <- utils::packageDescription(pkg)
    value <- function(key) if (is.null(d[[key]])) "" else as.character(d[[key]])
    data.frame(package = pkg, version = value("Version"), remote_repo = value("RemoteRepo"),
               remote_ref = value("RemoteRef"), remote_sha = value("RemoteSha"))
  })
  if (length(remote_info)) write_table(do.call(rbind, remote_info), file.path(out, "package_sources.csv"))
  sources <- c(list.files("R", full.names = TRUE), "scripts/run_pipeline.R")
  write_table(data.frame(file = sources, md5 = unname(tools::md5sum(sources))), file.path(out, "source_checksums.csv"))
})
