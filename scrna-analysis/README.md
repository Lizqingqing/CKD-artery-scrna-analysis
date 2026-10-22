# Single-cell RNA-seq analysis in R

A modular 10x → Seurat QC/clustering/markers → manual annotation → scMetabolism AUCell/KEGG → Monocle3 trajectory workflow. All code comments and documentation are in English.

## Provenance and version caveat

This repository was reconstructed from the requested workflow. The referenced conversation returned placeholders instead of the preceding code, and no detailed methods or input data were available. It therefore does not claim exact reproduction of an unseen manuscript method or line-by-line repair of that code. Thresholds are editable examples, not study-specific decisions.

**The requested scMetabolism v1.0 could not be verified.** The upstream DESCRIPTION inspected on 2026-10-22 reports **0.2.2**. The default config requests `1.0` and stops if the installed version differs. Supply the intended v1.0 source and validate its API, or deliberately change `expected_version` to the actual installed version and report that version in your methods. Do not relabel 0.2.2 as 1.0. See [environment guidance](docs/ENVIRONMENT.md).

## Layout

```text
R/                         Analysis functions, one module per stage
scripts/run_pipeline.R     Command-line entry point
scripts/install.R          Dependency bootstrap with explicit Git refs
config/config.yml          Parameters and run output directory
config/samples.csv         Sample manifest
config/annotations.csv     Reviewed cluster-to-cell-type mapping
scripts/check.R            Dependency-free syntax and guard checks
.github/workflows/check.yml Automated lightweight checks
results/                   Generated outputs, excluded from Git
data/raw/                 Local input, excluded from Git
```

## Input

Use filtered 10x gene-expression matrices: `matrix.mtx[.gz]`, `barcodes.tsv[.gz]`, and `features.tsv[.gz]` (legacy `genes.tsv` also supported by Read10X) in each manifest directory. Paths are relative to the repository root. Column 2 must contain gene symbols. Multimodal input selects `Gene Expression`. Raw integer counts are required; do not provide normalized expression.

Use unique sample IDs containing letters, digits, and hyphens. `donor_id` identifies the biological replicate; several libraries from one donor are not independent replicates. Barcodes are prefixed by sample ID. Read10X makes duplicate symbols unique; review symbol aliases and duplicated symbols before pathway scoring. No automatic Ensembl conversion or cross-species conversion is performed.

## Quick start

1. Install dependencies following [ENVIRONMENT.md](docs/ENVIRONMENT.md).
2. Put data in `data/raw/`; edit the sample manifest and configuration. Use `^MT-` for human or `^mt-` for mouse when appropriate to your annotation.
3. Run from the repository root:

```sh
Rscript scripts/check.R
Rscript scripts/run_pipeline.R config/config.yml standard
```

4. Inspect QC, clusters, sample mixing, and marker tables. Copy `results/run01/standard/annotation_template.csv` to `config/annotations.csv` and fill `cell_type` using multiple canonical markers and biological knowledge. Blank mappings are rejected; a header-only map is allowed and yields `Unassigned` for all clusters. Partial mappings yield `Unassigned` for missing clusters.

```sh
Rscript scripts/run_pipeline.R config/config.yml annotate results/run01/standard/clustered.rds
Rscript scripts/run_pipeline.R config/config.yml metabolism results/run01/annotate/annotated.rds
```

5. Configure a coherent lineage and biologically justified starting population before running trajectory. For example, after defining these labels in your annotation:

```yaml
trajectory:
  subset_column: cell_type
  subset_values: [Progenitor, Intermediate, Mature]
  root_column: cell_type
  root_values: [Progenitor]
  num_dim: 30
  use_partition: true
  graph_test: false
```

```sh
Rscript scripts/run_pipeline.R config/config.yml trajectory results/run01/annotate/annotated.rds
```

Each stage refuses to overwrite its output directory. For a revised stage or failed run, set a new `output_dir` (for example `results/run02`) and pass the desired previous-stage RDS explicitly. Scripts do not silently resume stale results. Preserve the printed console log if needed using your shell's standard redirection.

## Analysis choices and safeguards

- Sparse counts are retained through import and Seurat/Monocle3 setup. Seurat v5 layers are joined before extraction and marker testing.
- QC uses inclusive minimum/maximum gene thresholds, total-count ceiling, and mitochondrial fraction. Pre/post plots and a complete cell-level QC audit are saved; empty samples stop the run.
- LogNormalize (10,000), variable-feature selection, variable-gene scaling, PCA, neighbors, Louvain clustering, UMAP, and positive Wilcoxon markers form the standard workflow. PCA dimensions and neighbor counts are capped by the data size. Marker results are descriptive cluster markers, not donor-aware condition differential expression.
- Samples are merged without batch correction. Inspect sample UMAPs and donor composition. Integration, ambient RNA correction, doublet detection, and donor-level differential expression are not implemented; decide whether they are needed before scientific interpretation. No covariates are silently regressed out.
- scMetabolism uses its matrix API with `AUCell`, no imputation, and KEGG. This avoids assumptions about legacy Seurat assay slots. Raw gene-by-cell counts are supplied, and returned cells are explicitly matched. Scoring uses the package's AUCell defaults; record the installed AUCell version. The matrix API may require more memory than the rough four-matrix estimate; it is not a guaranteed memory limit.
- Bundled KEGG sets use human symbols. Mouse standard/trajectory analysis is supported with an appropriate regex, but metabolism stops until a separately reviewed human-ortholog workflow is implemented. Uppercasing mouse symbols is not ortholog mapping.
- Every pathway's gene overlap is exported and checked. Low coverage stops scoring; investigate identifiers before reducing thresholds. AUC scores describe relative signature enrichment, not measured metabolic flux.
- Metabolic summaries are calculated within sample and cell type, with cell counts. No cell-level hypothesis tests pretending cells are independent donors are included.
- Monocle3 builds its own PCA/UMAP, clusters, partitions, and graph from raw counts of the selected lineage. Seurat UMAP is not copied into a new, inconsistent graph. Roots are explicit; partitions without roots and unreachable cells stop the run. Multiple justified roots can be supplied. Pseudotime is relative ordering, not elapsed time or proof of causality.
- Optional `graph_test` measures expression association with graph structure (Moran's I), not necessarily monotonic change along pseudotime.
- Seeds, effective config, package versions, input-RDS checksum, source checksums, and session information are saved. Stages write `SUCCESS` only on completion; partial files in a `FAILED` directory are not final results.

## Outputs and validation

See [OUTPUTS.md](docs/OUTPUTS.md) and [VALIDATION.md](docs/VALIDATION.md). No fabricated biological outputs or datasets are included. An exact lockfile must be generated after a successful installation and representative full run; this archive does not pretend to contain a validated dependency lock.

## References

- [Seurat guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [Seurat v5 layers and accessors](https://satijalab.org/seurat/articles/seurat5_essential_commands.html)
- [scMetabolism repository and citation](https://github.com/wu-yc/scMetabolism)
- [scMetabolism package version declaration](https://github.com/wu-yc/scMetabolism/blob/main/DESCRIPTION)
- [Monocle3 trajectory documentation](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/)
- [AUCell documentation](https://bioconductor.org/packages/AUCell/)

Before public release, add your authorship, study citation, and a license appropriate for your contributions. No third-party source code or KEGG gene-set content is redistributed here; dependencies retain their own licenses.
