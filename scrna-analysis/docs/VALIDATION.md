# Validation report

## Performed during packaging

- Parsed all R files with the available R interpreter.
- Exercised dependency-free assertions and QC/species validation failures.
- Checked that the ZIP includes all repository files, including hidden configuration files.

The exact local check output is recorded in `validation-checks.txt`.

## Not performed

No 10x dataset was provided. Seurat, scMetabolism, Monocle3, and yaml were not installed in the packaging environment. Dependency installation, analysis execution, plotted output rendering, biological validation, and compatibility with a genuine scMetabolism v1.0 release were therefore not tested. The lightweight GitHub workflow also does not run the biological pipeline.

## Acceptance checks before scientific use

1. Restore/install dependencies in a clean environment; record exact source refs.
2. Run a representative 10x dataset through the standard stage; compare counts, sample IDs, QC thresholds, and retained cells with independent expectations.
3. Inspect marker evidence and annotate; ensure condition is not confounded with donor or batch.
4. Resolve the scMetabolism version discrepancy, check identifier coverage, and verify score dimensions/cell alignment on a small input before scaling up.
5. Choose a coherent lineage and justified root populations; inspect every partition and confirm finite pseudotime and sensible branching.
6. Repeat with relevant QC, dimension, and clustering choices to assess robustness. Run biological replicate-aware statistical analyses separately when required.
7. Save and restore a real renv lockfile and rerun before calling the environment reproducible.

Known implementation pitfalls addressed include direct access to obsolete Seurat count slots, merged layers left split, fixed PCA dimensions exceeding matrix rank, duplicated cross-sample barcodes, unreviewed cell-type assignment, normalized counts passed as raw counts, automatic trajectory roots, unrooted partitions, non-finite pseudotime, dense-memory surprises, species/gene-set mismatch, and silently misreported package versions.
