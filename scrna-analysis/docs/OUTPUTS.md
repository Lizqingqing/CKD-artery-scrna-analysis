# Output contract

Outputs are written beneath `<output_dir>/<stage>/`. Paths in the default configuration start with `results/run01/`.

| Stage | Main artifacts | Meaning |
|---|---|---|
| standard | qc_all_cells.csv, qc_summary.csv | QC metrics and retention audit |
| standard | qc_before.pdf, qc_after.pdf | QC distributions by sample |
| standard | qc.rds, clustered.rds | Raw-count-preserving Seurat checkpoints |
| standard | markers.csv | Positive cluster markers; may be empty if none pass |
| standard | clusters.pdf, samples.pdf, pca_elbow.pdf | Clustering and diagnostic plots |
| standard | annotation_template.csv | Cluster labels to review and fill |
| annotate | annotated.rds, cell_metadata.csv, annotation.pdf | Reviewed annotations, metadata, and UMAP |
| metabolism | kegg_gene_coverage.csv, kegg_checksum.txt | Identifier coverage and signature provenance |
| metabolism | kegg_aucell_scores.rds/.csv | Pathway-by-cell numeric AUC scores |
| metabolism | kegg_sample_celltype_means.csv | Sample/cell-type descriptive means |
| trajectory | trajectory_unordered.rds, trajectory_partitions.csv | Graph checkpoint and partition membership |
| trajectory | trajectory.rds, pseudotime.csv, trajectory_pseudotime.pdf | Rooted graph, cell ordering, and plot |
| trajectory | trajectory_graph_test.csv | Optional graph-associated expression test |
| all | config.yml, manifests, STATUS.txt, sessionInfo.txt, package_versions.csv, package_sources.csv, source_checksums.csv | Run provenance and status |
| downstream | input_checksum.txt | Input RDS path and MD5 |

RDS objects preserve internal data structures for downstream R analysis. MD5 is used for accidental-change detection, not security. Raw matrix checksums must be archived separately. Cell identifiers include the sample prefix; all downstream score alignment is by cell name, never assumed column order. Non-finite pseudotime is exported with a `finite` flag before failure so disconnected cells can be inspected.

The repository ships an empty results directory because no input dataset was provided. Results files may contain sensitive sample/cell metadata and are excluded from Git by default.
