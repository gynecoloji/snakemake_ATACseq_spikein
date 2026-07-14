# Configuration

This workflow is configured through two files in this directory:

- `config.yaml` — all workflow parameters (see below)
- `samples.csv` — the sample sheet

plus reference data you download into `ref/` (not tracked in git; see
[Reference data](#reference-data)).

## Sample sheet (`config/samples.csv`)

CSV with one row per sample and these columns:

| column      | description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| `sample_id` | Sample name. Raw reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`. |
| `type`      | Free-text condition label (e.g. `Control`, `NICD3`).                         |
| `group`     | Replicate group. Reproducibility handling is chosen from group size (below). |

Example:

```csv
sample_id,type,group
GSF4007-Control_1_S11,Control,group1
GSF4007-Control_2_S13,Control,group1
GSF4007-Control_3_S15,Control,group1
GSF4007-NICD3-V5_1_S12,NICD3,group2
GSF4007-NICD3-V5_2_S14,NICD3,group2
GSF4007-NICD3-V5_3_S16,NICD3,group2
```

**Per-group reproducibility** is derived automatically from the number of
replicates in each `group`:

- **≥ 3 replicates** → majority vote (a peak is kept if it recurs in ≥
  `consensus_min_replicates` replicates).
- **exactly 2 replicates** → IDR (`idr_threshold`).
- **1 replicate** → the sample's own peaks are used as-is.

## Parameters (`config/config.yaml`)

### Inputs / alignment

| key              | meaning                                                                                  |
|------------------|------------------------------------------------------------------------------------------|
| `samples_table`  | Path to the sample sheet (`config/samples.csv`).                                          |
| `adapter_r1` / `adapter_r2` | Optional. Leave commented to auto-detect adapters (`--detect_adapter_for_pe`); set to force specific sequences. |
| `human_fasta`    | Human genome FASTA. **Must be chr-prefixed UCSC** (chr1..chrX) to match the blacklist.    |
| `spikein_fasta`  | Spike-in genome FASTA (any species: Drosophila dm6, yeast, E. coli, …).                   |
| `spikein_prefix` | Prepended to spike-in chrom names before concatenation (avoids collisions, enables split).|
| `combined_index` | Bowtie2 index prefix for the built human+spike-in reference.                              |
| `align_chroms`   | Human chromosomes kept when building the index (`[]` = all).                              |
| `keep_chroms`    | Analysis keep-set for the final human BAM (must be a subset of `align_chroms`).           |
| `blacklist`      | ENCODE-style blacklist BED (chr-prefixed).                                                |

### Peaks / spike-in / consensus

| key                        | meaning                                                              |
|----------------------------|----------------------------------------------------------------------|
| `peak_types`               | Peak types to analyze (`["narrowPeak"]`; add `broadPeak` if needed). |
| `effective_genome_size`    | For deepTools RPGC normalization (hg38 default provided).            |
| `bin_size`                 | bigWig bin size (bp).                                                |
| `consensus_window`         | Fixed consensus peak width around each summit (bp).                  |
| `consensus_min_replicates` | Majority-vote threshold for ≥3-replicate groups.                     |
| `idr_threshold`            | IDR threshold for 2-replicate groups.                               |
| `idr_relaxed_pvalue`       | MACS2 `-p` for the relaxed calls used as IDR input.                 |
| `idr_top_n_peaks`          | Top-N relaxed peaks retained per replicate for IDR.                 |
| `keep_chroms_regex`        | Consensus chrom filter; keep consistent with `keep_chroms`.        |
| `macs2_genome`             | MACS2 `-g` (`hs` for human).                                        |

### QC

| key                              | meaning                                                        |
|----------------------------------|----------------------------------------------------------------|
| `gtf`                            | GENCODE GTF (chr-prefixed) for TSS enrichment.                 |
| `promoter_bed` / `enhancer_bed`  | BEDs for the reads-in-annotation QC (shipped in `ref/`).       |
| `spikein_pct_min` / `spikein_pct_max` | Expected spike-in read-fraction band (%) for the QC status flag. |

## Reference data

Genomes, indexes and large annotations are **not** shipped in the repo (they are
`.gitignore`d). Download / place them under `ref/` before running, matching the
paths in `config.yaml`:

- `ref/hg38.fa` — chr-prefixed UCSC human genome
- `ref/dm6.fa` (or another spike-in genome)
- `ref/hg38_blacklist_regions.bed` — ENCODE hg38 blacklist (shipped)
- `ref/gencode.v36.annotation.gtf` — GENCODE annotation (for TSS QC)
- `ref/hg38.2bit` — for `computeGCBias`
- `ref/picard.jar` — Picard (used by MarkDuplicates)

The combined Bowtie2 index (`ref/COMBINED/`) is built automatically by the
`build_combined_genome` rule from `human_fasta` + `spikein_fasta`.

See the top-level `README.md` for full setup and run instructions.
