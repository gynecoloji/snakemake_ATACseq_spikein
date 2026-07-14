# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.1](https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.1.0...v1.1.1) (2026-07-14)


### Documentation

* finalize citation metadata (author, ORCID, Zenodo DOI) ([d01753e](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/d01753e9003fe84e338a06020ad748c5c8ccdbf5))

## [Unreleased]

## [1.1.0] - 2026-07-14

Adds the Snakemake Workflow Catalog's auto-rendering features (tube map +
parameter table), config validation, and citation metadata.

### Added

- `.test/` executable test case (tiny placeholder inputs) so the catalog renders
  the workflow **tube map** — the rule graph, via snakevision — from
  `snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph`.
- `workflow/schemas/config.schema.yaml` — a JSON schema documenting every config
  parameter (type, default, description). The catalog renders it as a
  **parameter table**, and the workflow now runs `validate(config, ...)` on every
  invocation (fail-fast on bad/missing parameters; fills in defaults).
- `CITATION.cff` and a Zenodo DOI badge / citation section in the README.

### Changed

- Documentation now points at `config.schema.yaml` as the single source of truth
  for parameters; removed the duplicated parameter tables from `README.md` and
  `config/README.md`.
- Replaced the static workflow-diagram files (`atacseq_pipeline_flowchart.mermaid`,
  `ATACseq_workflow.svg`, `ATACseq_workflow_interactive.html`) with the
  catalog's auto-rendered tube map.

### Removed

- The three static workflow-diagram artifacts listed above (recoverable from git
  history).

## [1.0.0] - 2026-07-13

First release published to the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/).
This release restructures the project into the catalog's standardized
`workflow/` + `config/` layout and consolidates the two Snakefiles into a single
unified-DAG entry point, so the pipeline can be deployed into other projects with
`snakedeploy deploy-workflow` or imported via Snakemake's `module` / `use rule`
system.

It contains the complete ATAC-seq spike-in pipeline — concatenated (human +
spike-in) Bowtie2 alignment, filtering/dedup/blacklist, MACS2 peaks, Active
Motif–style spike-in normalization, reproducible fixed-width consensus peaks with
a fragment-count matrix, and a comprehensive QC stage (deepTools QC, FRiP, IDR,
library complexity, TSS enrichment, spike-in QC, and a self-contained interactive
HTML report). The changes below describe the layout migration this release
introduces.

### Added

- `workflow/Snakefile` — single standardized entry point. A unified DAG builds the
  primary stage and the QC report in dependency order; run subsets with the new
  `atacseq_all` and `qc_all` target rules.
- `.snakemake-workflow-catalog.yml` — Workflow Catalog metadata (mandatory
  `--use-conda`; conda / apptainer / apptainer+conda deployment) enabling the
  "standardized usage" tier.
- `config/README.md` — configuration reference (sample sheet, parameters,
  required reference data).
- `workflow/rules/common.smk` — shared config, sample/group tables, output
  directory constants, and helper functions used by both stages.
- Deployment via `snakedeploy deploy-workflow` and the `module` / `use rule`
  module system (see the README).
- A `conda:` environment for the `peak_summary` rule (previously ran in the base
  environment).

### Changed

- Restructured into the standardized Snakemake Workflow Catalog layout:
  - `envs/` → `workflow/envs/`
  - `ref/*.py`, `ref/*.R` → `workflow/scripts/`
  - `ref/config.yaml` → `config/config.yaml`; `ref/samples.csv` →
    `config/samples.csv` (with `samples_table` updated accordingly)
  - `snakefile_ATACseq` / `snakefile_ATAC_QC` rule bodies → `workflow/rules/atacseq.smk`
    and `workflow/rules/qc.smk`; the shared preamble was extracted to `common.smk`.
- Docker entrypoint, `docker-compose.yml`, and `run_pipeline.sh` now target
  `-s workflow/Snakefile`; `.dockerignore` updated for the new layout. The
  pre-baked per-rule conda environments continue to be reused (env file contents
  are unchanged).
- `README.md` and `DOCKER.md` updated with the new run commands and a
  "Deploying with snakedeploy" section; test import paths point at
  `workflow/scripts/`.
- `ref/` now holds reference **data** only (genomes, annotations, BEDs).

### Removed

- Root `snakefile_ATACseq` and `snakefile_ATAC_QC` — superseded by
  `workflow/Snakefile` together with the `atacseq_all` / `qc_all` targets.

### Migration notes

- **Breaking:** invocations using `-s snakefile_ATACseq` / `-s snakefile_ATAC_QC`
  no longer work. Use `snakemake --use-conda -s workflow/Snakefile` (optionally
  with the `atacseq_all` or `qc_all` target). The two stages, previously run as
  two separate commands, now build in one dependency-ordered run by default.
- Configuration moved from `ref/` to `config/`; point your edits at
  `config/config.yaml` and `config/samples.csv`.

[Unreleased]: https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/gynecoloji/snakemake_ATACseq_spikein/releases/tag/v1.0.0
