score_metabolism <- function(obj, cfg, out) {
  need(c("scMetabolism", "AUCell", "GSEABase"))
  actual <- as.character(utils::packageVersion("scMetabolism"))
  assert(utils::compareVersion(actual, cfg$metabolism$expected_version) == 0,
         paste("scMetabolism version mismatch. Requested", cfg$metabolism$expected_version,
               "but installed", actual, "See docs/ENVIRONMENT.md; do not silently relabel versions."))
  assert(cfg$species == "human", "Bundled KEGG sets use human symbols. Mouse requires curated ortholog mapping first.")
  x <- counts_of(obj)
  # The upstream matrix API densifies counts; budget for several temporary matrices.
  estimated_gb <- as.double(nrow(x)) * ncol(x) * 8 * 4 / 1024^3
  assert(estimated_gb <= cfg$metabolism$max_dense_gb,
         paste("Estimated dense working memory exceeds budget:", round(estimated_gb, 2), "GiB."))
  gmt <- system.file("data", "KEGG_metabolism_nc.gmt", package = "scMetabolism")
  assert(nzchar(gmt) && file.exists(gmt), "Bundled KEGG GMT was not found.")
  sets <- GSEABase::getGmt(gmt)
  ids <- GSEABase::geneIds(sets)
  coverage <- data.frame(pathway = names(ids), size = lengths(ids),
                         matched = vapply(ids, function(g) sum(g %in% rownames(x)), integer(1)))
  coverage$fraction <- coverage$matched / coverage$size
  write_table(coverage, file.path(out, "kegg_gene_coverage.csv"))
  assert(all(coverage$matched >= cfg$metabolism$min_pathway_genes &
               coverage$fraction >= cfg$metabolism$min_pathway_fraction),
         "At least one KEGG pathway has low coverage. Inspect coverage before explicitly adjusting thresholds.")
  writeLines(paste("KEGG GMT MD5:", tools::md5sum(gmt)), file.path(out, "kegg_checksum.txt"))
  set.seed(cfg$seed)
  # Use the matrix API to avoid legacy Seurat slot access in the upstream wrapper.
  scores <- as.matrix(scMetabolism::sc.metabolism(countexp = x, method = "AUCell",
                         imputation = FALSE, ncores = 1, metabolism.type = "KEGG"))
  assert(!anyDuplicated(colnames(scores)) && setequal(colnames(scores), colnames(obj)),
         "Metabolism output cell IDs do not match the input.")
  scores <- scores[, colnames(obj), drop = FALSE]
  assert(all(is.finite(scores)), "Non-finite metabolic scores.")
  saveRDS(scores, file.path(out, "kegg_aucell_scores.rds"))
  write_table(data.frame(pathway = rownames(scores), scores, check.names = FALSE),
              file.path(out, "kegg_aucell_scores.csv"))
  md <- obj[[]]
  groups <- interaction(md$sample_id, md$cell_type, drop = TRUE)
  summaries <- lapply(split(seq_len(ncol(obj)), groups), function(i) {
    data.frame(sample_id = md$sample_id[i[1]], cell_type = md$cell_type[i[1]], n_cells = length(i),
               pathway = rownames(scores), mean_auc = rowMeans(scores[, i, drop = FALSE]))
  })
  write_table(do.call(rbind, summaries), file.path(out, "kegg_sample_celltype_means.csv"))
  invisible(scores)
}
