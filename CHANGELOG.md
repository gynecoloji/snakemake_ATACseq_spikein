# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.1](https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.4.0...v1.4.1) (2026-07-23)


### Documentation

* **diffopen:** add diffopen_rna_table format + example table ([9bca79d](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/9bca79d2386e43ae31a9793f2e7ff81389b28efa))
* **diffopen:** add diffopen_rna_table format + example table ([bda5df1](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/bda5df1e8a116dda099b4a17f78023e25ddd620b))

## [1.4.0](https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.3.0...v1.4.0) (2026-07-23)


### Added

* **diffopen:** add rnastable normalization mode ([359874a](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/359874a44d0f4f343b9bffd37d10a3fc96160b1a))
* **diffopen:** add rnastable wildcard + params wiring to diffopen.smk ([bcc8982](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/bcc8982685053b428f30556909847dae8d5c0dc1))
* **diffopen:** add size_factors_rnastable orchestrator + anchor floor ([c33ace2](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/c33ace21df7b7f560906271432f058baa329ddf3))
* **diffopen:** add stable_genes_from_de for rnastable mode ([ed7a4ad](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/ed7a4ad729c6462b4e38a13316ae41da92aaa89b))
* **diffopen:** add TSS-window + anchor selection for rnastable ([03680e5](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/03680e5aa41995fa7fe4173159851ec795e61250))
* **diffopen:** config schema + fast-fail guard for rnastable mode ([808ba76](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/808ba76f1d3b25f307614e05e8ad4c8ff0a0a2fe))
* **diffopen:** include rnastable in the normalization comparison report ([a424939](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/a4249391f3324371991a3235e6988400a2af562f))
* **diffopen:** wire rnastable mode into diffopen.R main() ([febf1b8](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/febf1b8ee71eac1750f0e163808e29a0004275f3))


### Changed

* **diffopen:** final-review cleanups (test path, stale comments, friendly rna_table error) ([a437ca5](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/a437ca5f4aa4cd9c6f5628407b7d21e529a0f267))


### Documentation

* document the rnastable diffopen normalization mode ([86d015d](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/86d015d4e136ecef515c6bf418c4a58d6a51a5e7))

## [1.3.0](https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.2.0...v1.3.0) (2026-07-21)


### Added

* add downstream gene/enrichment/track stages and split anchor_shape by class ([73bf28b](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/73bf28bb0f0ea5dcf0b9ca78bfda62873308add0))
* mode-scaled bigWigs for tracks, and pairwise IDR in the QC report ([a469d49](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/a469d49be1af59f6171d28d9896f26f0f8e29eaa))
* split differential openness by promoter vs enhancer class ([074c182](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/074c182d57d4aee8d221c2741e69c0d73349f34e))


### Fixed

* declare diffopen scripts as rule inputs and repair valueless-flag parsing ([cd648d5](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/cd648d52addc02f73017670d3d1f91ce7f23b46e))
* declare shell-invoked scripts as inputs in qc.smk and atacseq.smk ([5cbd0fc](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/5cbd0fcd79b5de192c7d405b9909fa6799f156d0))
* report unshrunk MLE effect sizes alongside the shrunk log2FC ([c1d6444](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/c1d64440c91c8db15864cb8730e17a6f5b0931c9))
* second bare {input} in samtools_sort_filter_index ([d541a24](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/d541a2474c8b7f0b166d54a7a9dc3f3847b63ddc))


### Documentation

* state that Gviz track heights are not the effect size ([7923587](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/7923587a47f548d90e4c03367866eec2a245d0c9))

## [1.2.0](https://github.com/gynecoloji/snakemake_ATACseq_spikein/compare/v1.1.1...v1.2.0) (2026-07-19)


### Added

* add differential-openness stage with selectable normalization ([f24bee8](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/f24bee8429f24e26d64492981a9dcd714881c7ba))
* add native Apptainer build definition for the sif ([f93d510](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/f93d510066c1fb7d9982b88d3e29c396f23eb9f0))


### Documentation

* embed tube-map diagram in README ([34864d4](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/34864d4366dd30d635a05eb546c3377757702bd0))
* fix Code of Conduct contact email ([899e385](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/899e385b15027a9cbd1a2c6869910c5014bfd3bb))
* replace rule graph with snakevision tube map ([c63b04b](https://github.com/gynecoloji/snakemake_ATACseq_spikein/commit/c63b04bf9691666794821582809843ff25fccf87))

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
