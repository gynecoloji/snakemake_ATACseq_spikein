# Technical documentation

Step-by-step documentation of the `snakemake_ATACseq_spikein` workflow. For
installation, container usage, and the full narrative, see the top-level
[`README.md`](../README.md); for every configuration parameter, see
[`config/README.md`](../config/README.md) and the schema
[`workflow/schemas/config.schema.yaml`](schemas/config.schema.yaml).

The rule graph is shown in [`images/rulegraph.svg`](../images/rulegraph.svg) and
rendered as a "tube map" on the workflow's Snakemake Workflow Catalog page.

## Overview

A single `snakemake -s workflow/Snakefile --use-conda` run builds a **unified
DAG** covering two stages in dependency order:

1. **Primary** (`atacseq_all` target) — alignment → filtering → peak calling →
   spike-in normalization → reproducible consensus peaks + fragment counts.
2. **QC** (`qc_all` target) — deepTools QC, FRiP, IDR, library complexity,
   spike-in QC, TSS enrichment, and an interactive HTML QC report.

## Inputs

| Input | Location | Notes |
|---|---|---|
| Paired-end reads | `data/<sample_id>_R1_001.fastq.gz`, `_R2_001.fastq.gz` | one pair per sample |
| Sample sheet | `config/samples.csv` | columns `sample_id, type, group` |
| Human genome FASTA | `ref/hg38.fa` | chr-prefixed UCSC |
| Spike-in genome FASTA | `ref/dm6.fa` | any species |
| Blacklist BED | `ref/hg38_blacklist_regions.bed` | ENCODE, chr-prefixed |
| GTF / 2bit / promoter+enhancer BEDs | `ref/…` | QC references |
| Picard | `ref/picard.jar` | duplicate marking |

Configuration is read from `config/config.yaml` and validated against the schema
at parse time (missing/invalid parameters fail fast).

## Steps (primary stage)

1. **`fastqc`** — raw-read quality.
2. **`fastp`** — adapter trimming + quality filtering (auto-detects adapters).
3. **`build_combined_genome`** — prefix the spike-in chromosomes, concatenate
   with the (optionally chromosome-subset) human genome, and build one Bowtie2
   index.
4. **`bowtie2_align`** — single alignment pass to the combined genome; each read
   is assigned to exactly one genome by chromosome prefix.
5. **`samtools_sort_filter_index`** — keep uniquely-mapped, properly-paired human
   reads; record mitochondrial-% QC; restrict to the analysis chromosomes.
6. **`remove_duplicates`** — Picard MarkDuplicates.
7. **`filter_blacklist`** — fragment-level ENCODE blacklist removal.
8. **`call_peaks`** — MACS2 (BAMPE, `-q 0.05`).
9. **Spike-in normalization** (`spikein_extract_count`, `compute_spikein_factors`,
   `create_bigwig`, `create_spikein_bigwig`) — count spike-in reads, derive the
   normalization factor `NF = min(spike-in reads) / spike-in reads`, and emit
   depth-normalized and spike-in-scaled bigWigs.
10. **Consensus peaks** (`relaxed_peaks`, `reproducible_idr`, `consensus_peaks`,
    `count_fragments_consensus`) — per-group reproducibility (majority vote for
    ≥3 replicates, IDR for exactly 2), a fixed-width consensus set, and a
    featureCounts fragment matrix.

## Steps (QC stage)

deepTools coverage/fragment-size/fingerprint/correlation/PCA/GC/TSS, a numeric
TSS-enrichment score, FRiP, IDR on relaxed peaks, library complexity
(NRF/PBC1/PBC2), spike-in fraction QC, reads-in-annotation and peak summaries, a
FastQC-only MultiQC report, and a self-contained interactive HTML QC report
(`results/qc/atacseq_qc_report.html`).

## Outputs

All outputs are written under `results/` (peaks, bigWigs, consensus matrix, QC
tables and reports); per-rule logs under `logs/`. See the README's "Output
Structure" section for the full tree.

## Running the tests

```bash
python -m pytest tests/ -q                               # unit tests
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph   # DAG/tube map
```
