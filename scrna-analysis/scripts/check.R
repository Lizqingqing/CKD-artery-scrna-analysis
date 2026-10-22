#!/usr/bin/env Rscript
# These dependency-free checks exercise parsing and configuration guard behavior.
files <- c(list.files("R", pattern = "\\.R$", full.names = TRUE),
           list.files("scripts", pattern = "\\.R$", full.names = TRUE))
for (f in files) parse(f)
source("R/utils.R")
expect_error <- function(expr) {
  caught <- inherits(tryCatch({ force(expr); NULL }, error = identity), "error")
  stopifnot(caught)
}
assert(TRUE, "Valid assertion unexpectedly failed")
expect_error(assert(FALSE, "Expected failure"))
expect_error(assert(NA, "Expected failure"))
cfg <- list(species = "human", qc = list(min_features = 200, max_features = 6000,
             max_counts = 50000, max_percent_mt = 20),
            analysis = list(pcs = 30, variable_features = 2000, resolution = 0.5),
            metabolism = list(max_dense_gb = 4))
validate_config(cfg)
bad <- cfg; bad$qc$max_features <- 100
expect_error(validate_config(bad))
bad <- cfg; bad$qc$max_percent_mt <- 101
expect_error(validate_config(bad))
bad <- cfg; bad$species <- "unknown"
expect_error(validate_config(bad))
cat("PASS:", length(files), "R files parsed; assertion and invalid-QC/species guards checked.\n")
