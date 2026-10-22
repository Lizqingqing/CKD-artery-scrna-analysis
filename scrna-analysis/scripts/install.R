#!/usr/bin/env Rscript
# Supply reviewed immutable refs; installation is deliberately separate from analysis.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || any(!nzchar(args))) {
  stop("Usage: Rscript scripts/install.R SC_METABOLISM_REF MONOCLE3_REF VISION_REF")
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("remotes", "BiocManager", "renv", "yaml", "Seurat", "ggplot2"))
BiocManager::install(c("AUCell", "GSEABase", "GSVA", "SingleCellExperiment",
                       "SummarizedExperiment", "BiocGenerics", "DelayedArray",
                       "DelayedMatrixStats", "batchelor", "HDF5Array"),
                     ask = FALSE, update = FALSE)
remotes::install_github(paste0("YosefLab/VISION@", args[3]), upgrade = "never")
remotes::install_github(paste0("wu-yc/scMetabolism@", args[1]), upgrade = "never")
remotes::install_github(paste0("cole-trapnell-lab/monocle3@", args[2]), upgrade = "never")
print(sessionInfo())
