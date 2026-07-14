## Citation

If you use this pipeline in your research, please cite:

```
https://github.com/gynecoloji/snakemake_ATACseq_spikein
```

**Please also cite the individual tools used:**
- **Snakemake**: Köster, J. and Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.
- **MACS2**: Zhang, Y. et al. (2008). Model-based analysis of ChIP-Seq (MACS). Genome Biology, 9, R137.
- **deepTools**: Ramírez, F. et al. (2016). deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research, 44(W1), W160-W165.
- **Bowtie2**: Langmead, B. and Salzberg, S.L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-359.
- **SAMtools**: Li, H. et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079.
- **IDR**: Li, Q. et al. (2011). Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, 5(3), 1752-1779.
- **featureCounts (Subread)**: Liao, Y. et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.
- **Consensus peaks (fixed-width / score-per-million)**: Corces, M.R. et al. (2018). The chromatin accessibility landscape of primary human cancers. Science, 362(6413), eaav1898.
- **BEDTools**: Quinlan, A.R. and Hall, I.M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

**Author**: gynecoloji  
**Project Repository**: [https://github.com/gynecoloji/snakemake_ATACseq_spikein](https://github.com/gynecoloji/snakemake_ATACseq_spikein)

For questions, issues, or feature requests, please:
1. Check the existing [Issues](https://github.com/gynecoloji/snakemake_ATACseq_spikein/issues) on GitHub
2. Submit a new issue with detailed information about your problem
3. Include relevant log files and system information for troubleshooting

## Acknowledgments

This pipeline was developed based on best practices from the ENCODE consortium and incorporates methodologies from multiple published ATAC-seq analysis workflows. Special thanks to the developers of all the integrated tools that make this comprehensive analysis possible.

---

**Note**: This pipeline is optimized for human genome analysis (hg38) but can be adapted for other organisms by updating reference files and parameters accordingly. It supports Active Motif–style spike-in normalization (e.g. Drosophila) via a concatenated human + spike-in Bowtie2 alignment.# ATAC-seq Analysis Pipeline

A comprehensive Snakemake workflow for processing and analyzing ATAC-seq data from raw reads to peak calling with extensive quality control metrics and differential binding analysis.

## Overview

This pipeline integrates three complementary components for complete ATAC-seq analysis:

1. **Primary ATAC-seq workflow** (`snakefile_ATACseq`) - Raw FASTQ → concatenated (human + spike-in) Bowtie2 alignment → filtering → MACS2 peak calling → **spike-in normalization** (scaled bigWigs) → a reproducible, fixed-width **consensus peak set** with a fragment-count matrix
2. **ATAC-seq QC workflow** (`snakefile_ATAC_QC`) - deepTools QC, FRiP, IDR, library complexity, spike-in QC, TSS enrichment score, a self-contained **interactive HTML QC report** (all QC except FastQC), plus a FastQC-only MultiQC
3. **Differential Analysis Notebook** (`ATACseq_Dx.ipynb`, R / Bioconductor) - **DESeq2** differential binding (NICD3 vs Ctrl) on the consensus count matrix — split into **promoter vs distal** peaks with a **paired** design and run under **both spike-in and default (median-of-ratios) normalization** for comparison — plus **Gviz genome-browser tracks** at Notch target loci as a positive control; runs on the `atacseq-diffbind` container

## Workflow Diagram

The complete workflow is shown below:

![ATAC-seq Workflow Diagram](ATACseq_workflow.svg)

## Features

- **Complete end-to-end processing** of paired-end ATAC-seq data
- **Spike-in normalization** (Active Motif–style) via a concatenated human + spike-in genome, with spike-in-scaled bigWig tracks
- **Concatenated Bowtie2 alignment** with per-read genome assignment (no cross-mapping double-counting)
- **Chromosome control**: align to a configurable set (default chr1–22, chrX, chrM), record mitochondrial-% QC, then keep only the analysis chromosomes
- **Reproducible consensus peaks** (Corces-2018 fixed-width / score-per-million; majority-vote or IDR reproducibility) + a **featureCounts** fragment matrix for differential analysis
- **Blacklist region filtering** for removal of technical artifacts
- **Extensive QC metrics** (fragment size, TSS heatmap + numeric enrichment score, fingerprint, correlation/PCA, GC bias, FRiP, NRF/PBC1/PBC2, spike-in %, reads in promoters/enhancers)
- **Interactive HTML QC report** (self-contained, theme-aware; all QC except FastQC) + a FastQC-only MultiQC
- **Signal track generation** (bigWig, bedGraph) for visualization
- **Conda environment management** (one env per rule) plus a ready-to-run **Docker / Apptainer** image

## Pipeline Components

### 1. Primary Processing Pipeline (`snakefile_ATACseq`)

**Processing Steps:**
```
Raw FASTQ → FastQC → fastp
  → build combined (human + spike-in) Bowtie2 index
  → Bowtie2 alignment (one pass to the combined genome)
     ├─ human reads   → unique + properly-paired filter → mito-% QC → keep chr1–22/chrX
     │                  → Picard dedup → blacklist filter → MACS2 peaks
     │                  → spike-in-scaled + depth-normalized bigWigs
     │                  → consensus peaks → featureCounts fragment matrix
     └─ spike-in reads → dedup → count → normalization factor
```

**Key Features:**
- Quality assessment with FastQC; read trimming/filtering with fastp
- **Concatenated alignment** to human + spike-in with Bowtie2; reads are split by genome (spike-in chroms are prefixed), so each read is assigned to exactly one genome
- Uniquely-mapped, properly-paired filtering; filtering-orphaned mates removed; **mitochondrial-% recorded** then non-primary chromosomes dropped
- Picard duplicate removal + fragment-level ENCODE blacklist filtering
- MACS2 peak calling (BAMPE)
- **Spike-in normalization**: `NF = min(spike-in reads) / spike-in reads`, applied as a bigWig scale factor
- **Consensus peaks** (fixed-width, SPM-ranked; majority-vote/IDR reproducibility) + **featureCounts** matrix

### 2. Quality Control Pipeline (`snakefile_ATAC_QC`)

**Comprehensive QC Metrics** (run after the primary pipeline; consumes its `results/`):
- **Fragment Size Analysis** - Insert size distribution and nucleosomal patterns
- **TSS Enrichment** - Heatmap/profile **and a numeric per-sample enrichment score**
- **Signal Quality** - Fingerprint plots for signal-to-noise assessment
- **Sample Correlation** - Multi-sample correlation heatmap/scatter and PCA
- **GC Bias Assessment** - Sequence composition bias evaluation
- **FRiP Scores** - Fraction of Reads in Peaks
- **Library Complexity** - PCR bottlenecking assessment (NRF, PBC1, PBC2)
- **IDR Analysis** - Irreproducible Discovery Rate on relaxed peak calls, per replicate pair
- **Spike-in QC** - Spike-in % per sample vs the Active Motif 2–10% target
- **Peak + annotation summary** - Peak counts/widths and reads in promoters vs enhancers
- **Interactive HTML QC report** (`atacseq_qc_report.html`) - one self-contained, theme-aware page for all QC **except FastQC** (alignment rate, mito %, duplication, blacklist, spike-in %/norm factors, peaks/FRiP, TSS enrichment, NRF/PBC1/PBC2, fragment size/GC/correlation/PCA/fingerprint, reads-in-annotation, consensus), with ENCODE-threshold pass/warn/fail flags; numeric metrics render as interactive tables/bar charts and the deepTools QC (fragment size, GC bias, correlation/PCA, fingerprint, TSS profile + per-region heatmap) renders as **interactive SVG/canvas charts drawn client-side** — nothing is embedded as a static image. Mitochondrial % comes from the primary pipeline's `idxstats`
- **FastQC MultiQC report** (`multiqc_fastqc.html`) - MultiQC scoped to FastQC (raw-read quality) only

### 3. Differential Binding Analysis (`ATACseq_Dx.ipynb`)

An **R / Bioconductor** notebook (runs on the `atacseq-diffbind` container, `ir` kernel) that
takes the pipeline's consensus fragment matrix (`results/consensus/consensus_counts.txt`) and
performs the differential test end-to-end:

**Differential binding (DESeq2):**
- **NICD3 vs Ctrl** on the consensus peaks, split into **promoter vs distal** sets (overlap with
  `ref/promoter_chr1-22X.bed`) and tested separately
- **Paired design** `~ pair + condition` — each replicate index blocks pair-to-pair variation
- Run under **two normalizations, compared side by side**: the pipeline's **spike-in** scale
  factors (`sizeFactor = 1/NF`, centered) and DESeq2's **default** median-of-ratios (`sf = NULL`)
- Per group × normalization: full results + significant subset (`padj < 0.05 & |log2FC| > 1`),
  MA/volcano plots, sample **PCA** (spike-in-normalized VST), and ChIPseeker nearest-gene
  annotation → `results/diff_region/`

**Positive-control genome-browser tracks (Gviz):**
- Signal at canonical Notch/NICD targets (HES1, HEY1, HEYL, NRARP, DTX1) across all samples, from
  **both** the spike-in-scaled and RPGC bigWigs (**two figures per gene**), with all-transcript
  gene models from the GTF → `results/browser_tracks/` (PNG + PDF)

The executed notebook is generated from `ref/build_diffbind_notebook.py` (helpers in
`ref/diffbind_helpers.R`).

## Requirements

The pipeline requires the following dependencies:

- [Snakemake](https://snakemake.readthedocs.io/) ≥7.0.0
- [Conda](https://docs.conda.io/en/latest/) / [Mamba](https://github.com/mamba-org/mamba) (recommended)
- [Python](https://www.python.org/) ≥3.8
- UNIX-based system (Linux/MacOS)

### Software Dependencies
(automatically installed via conda environments):
- **FastQC** (quality control)
- **fastp** (read trimming)
- **Bowtie2** (alignment to the concatenated human + spike-in genome)
- **SAMtools** (BAM processing, read splitting, idxstats)
- **Picard** (duplicate removal)
- **MACS2** (peak calling)
- **deepTools** (bigWig/bedGraph, QC and visualization)
- **bedtools** + **bc** (genomic interval operations; FRiP / complexity / annotation)
- **IDR** (reproducibility analysis)
- **featureCounts / Subread** (fragment quantification over the consensus set)

### Differential-analysis environment (R / Bioconductor)
`ATACseq_Dx.ipynb` is an **R** notebook (`ir` kernel), not Python. Run it on the
**`atacseq-diffbind`** container, which bundles R 4.6 / Bioconductor with everything it needs:
- **DESeq2**, **apeglm** — differential binding + log2FC shrinkage
- **GenomicRanges** / **IRanges**, **rtracklayer** — region handling, bigWig I/O
- **ChIPseeker**, **TxDb.Hsapiens.UCSC.hg38.knownGene**, **org.Hs.eg.db** — annotation
- **Gviz** — genome-browser tracks; **ggplot2** — MA / volcano / PCA

(edgeR / limma / DiffBind are also installed for alternative differential frameworks.)

### Reference Files (`ref/`)

Files fall into three groups: **shipped** with the repo (present after clone),
**downloaded** from public sources, and **generated** locally from the downloads.
Exact commands are in [Reference files: download & generate](#reference-files-download--generate).

```
ref/
├── config.yaml                          # pipeline configuration                 (shipped)
├── samples.csv                          # sample sheet: sample_id,type,group      (shipped; you edit)
│
│  ── downloaded (public sources; not in the repo) ──
├── hg38.fa                              # human genome FASTA (UCSC, chr-prefixed)
├── dm6.fa                               # spike-in genome FASTA (e.g. Drosophila)
├── hg38.2bit                            # human genome 2bit (QC: GC bias)
├── gencode.v36.annotation.gtf           # GENCODE annotation (TSS, gene models)
├── FANTOM5_CAGE_peaks_hg38.bed.gz       # FANTOM5 CAGE peaks (optional; only to rebuild the CAGE set)
├── picard.jar                           # Picard (MarkDuplicates)
│
│  ── shipped, or generated from the downloads ──
├── hg38_blacklist_regions.bed           # ENCODE hg38 blacklist v2                (shipped)
├── hg38.fa.fai                          # faidx of hg38.fa                        (generated)
├── promoter_chr1-22X.bed                # Ensembl Regulatory Build promoters      (shipped; QC reads-in-annotation)
├── enhancer_chr1-22X.bed                # Ensembl Regulatory Build enhancers      (shipped; QC reads-in-annotation)
├── Promoter_uniqueTSS_hg38_{3000,5000}bp_chr1-22X.bed      # TSS±N, unique TSS          (shipped)
├── Promoter_MANEcanonical_hg38_{3000,5000}bp_chr1-22X.bed  # TSS±N, one per gene (MANE) (shipped)
├── Promoter_FANTOM5CAGE_hg38_{3000,5000}bp_chr1-22X.bed    # TSS±N, empirical CAGE      (shipped)
│      # per-transcript set (Promoter_UCSC_hg38_*bp) not shipped — `build_promoter_beds.py transcript`
├── COMBINED/                            # combined human+spike-in Bowtie2 index   (built by the pipeline)
│
│  ── shipped scripts ──
├── build_promoter_beds.py               # regenerates the promoter TSS BEDs (4 definitions)
├── build_diffbind_notebook.py           # regenerates ATACseq_Dx.ipynb
├── build_qc_report.py  tss_score.py  consensus_peaks.py  process_sam.py
├── compute_spikein_factors.py  blacklist-stats-script.py  downsample_tss_matrix.py
└── diffbind_helpers.R                   # R helpers for the differential-binding notebook
```

### Reference files: download & generate

Only the **downloaded** files are absent after a clone. Fetch them into `ref/`, then
build the small derived files. Everything else ships with the repo.

**1. Download** (run from the repo root):

```bash
cd ref

# Human + spike-in genomes and the human 2bit (UCSC goldenPath)
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz   && gunzip dm6.fa.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# GENCODE v36 gene annotation
curl -O http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
gunzip gencode.v36.annotation.gtf.gz

# FANTOM5 robust CAGE peaks (hg38) — OPTIONAL: only needed to rebuild the CAGE
# promoter set (the built Promoter_FANTOM5CAGE_*.bed already ships with the repo)
curl -L -o FANTOM5_CAGE_peaks_hg38.bed.gz \
  "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"

# Picard (MarkDuplicates)
curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/latest/download/picard.jar

cd ..
```

The ENCODE blacklist ships with the repo; to refresh it (Boyle Lab v2):

```bash
curl -L https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
  | gunzip > ref/hg38_blacklist_regions.bed
```

**2. Generate** (needs the downloads above; run from the repo root):

```bash
# faidx of the human genome (used to clamp promoter windows to chromosome ends)
samtools faidx ref/hg38.fa

# The promoter TSS BEDs (unique / MANE / CAGE) already ship with the repo. Rebuild
# them only after updating the GTF, or to also produce the per-transcript set:
python ref/build_promoter_beds.py all          # modes: unique | mane | cage | transcript | all
# (needs ref/hg38.fa.fai + the GTF; the `cage` mode also needs the FANTOM5 download above)
# Kept chromosomes default to `keep_chroms` in ref/config.yaml (currently chr1-22,X);
# override per run with e.g. --chroms chr1,chr2,chrX  or  --chroms all  (the filename
# scope token, e.g. chr1-22X, auto-adjusts; set it explicitly with --label).

# The combined human+spike-in Bowtie2 index (ref/COMBINED/) is built automatically
# by the pipeline's build_combined_genome rule from hg38.fa + dm6.fa — no manual step.
```

The QC annotation BEDs (`promoter_chr1-22X.bed`, `enhancer_chr1-22X.bed`) ship with the
repo. They are **Ensembl Regulatory Build** features (`feature_type` Promoter / Enhancer,
`ENSR*` IDs), chr-prefixed and restricted to chr1–22,X. To rebuild from a newer release,
download the regulatory GFF from
`http://ftp.ensembl.org/pub/current_regulation/homo_sapiens/`, keep the Promoter /
Enhancer features, and emit `chrom start end ENSR_id`.

### Promoter / TSS definitions

`build_promoter_beds.py` writes promoter windows (TSS ±3kb and ±5kb, strand-aware,
chr1–22,X) under four TSS definitions so you can match the definition to the analysis.
All windows are symmetric width 2N; GTF-derived sets keep a 10-column layout
(`tx_id chrom start end score strand tx_id.ver transcript_type gene_name gene_id`),
CAGE is BED6.

| `mode` | TSS = | Windows (chr1–22,X) | Use when |
|---|---|---|---|
| `transcript` | every transcript's 5′ end (GENCODE v36) | 231,104 | you want all annotated alternative promoters (one per isoform; ChIPseeker-style) |
| `unique` | distinct (chrom, TSS, strand) tuples | 205,696 | as `transcript`, without isoform double-counting |
| `mane` | one canonical TSS per gene (MANE Select, else 5′-most) | 60,058 | one row per gene; ENCODE-comparable TSS-enrichment reference |
| `cage` | FANTOM5 robust CAGE peak position (empirical) | 209,146 | data-defined start sites, incl. promoters annotation misses |

The repo ships the **`unique`, `mane`, and `cage`** sets (each at ±3kb and ±5kb); the
per-transcript set is generate-on-demand (`build_promoter_beds.py transcript`).

Cross-check: **90.6%** of protein-coding MANE-canonical TSSs have a FANTOM5 CAGE peak
within 500 bp (strong annotation↔empirical agreement); non-coding genes only **24.9%**
(lncRNAs etc. are largely CAGE-silent). Note that each transcript has exactly one TSS,
whereas a gene can have several — the pipeline's deepTools TSS-enrichment QC and the
`ChIPseeker` peak annotation both use **transcript-level** TSS by default.

## Installation

```bash
# Clone the repository
git clone https://github.com/gynecoloji/snakemake_ATACseq_spikein.git
cd snakemake_ATACseq_spikein

# You need Snakemake + conda/mamba as the driver. The per-rule tool environments
# (envs/*.yaml) are created automatically on the first `--use-conda` run — you do
# not build them by hand.
mamba create -n atacseq -c conda-forge -c bioconda snakemake-minimal pandas
conda activate atacseq
```

> **No local install?** Skip all of the above and use the container image instead —
> see [Container Execution (Docker / Apptainer)](#container-execution-docker--apptainer).

## Configuration

Edit `ref/config.yaml` to match your experimental setup and reference files:

```yaml
# ── Samples ──
samples_table: "ref/samples.csv"          # columns: sample_id, type, group

# ── Adapter trimming (fastp) ──
# Omit or leave empty  ->  fastp AUTO-DETECTS adapters for PE reads (--detect_adapter_for_pe).
# Provide both to OVERRIDE auto-detection with explicit sequences (Illumina shown):
# adapter_r1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
# adapter_r2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ── Alignment: concatenated (human + spike-in) reference ──
human_fasta:    "ref/hg38.fa"             # chr-prefixed UCSC hg38 (you provide)
spikein_fasta:  "ref/dm6.fa"              # spike-in genome, ANY species (you provide)
spikein_prefix: "spikein_"                # prepended to spike-in chrom names (to split reads)
combined_index: "ref/COMBINED/genome"     # built by the pipeline
align_chroms: [chr1, ..., chr22, chrX, chrM]    # human chroms that enter the index
keep_chroms:  [chr1, ..., chr22, chrX]          # analysis chroms in the final BAM (drop chrM/chrY)

# ── Peaks / blacklist / bigWig ──
blacklist: "ref/hg38_blacklist_regions.bed"
macs2_genome: "hs"
effective_genome_size: 2913022398
bin_size: 25

# ── Consensus peaks & IDR ──
consensus_window: 500                     # fixed peak width (bp)
consensus_min_replicates: 2               # majority-vote threshold (>=3-replicate conditions)
idr_threshold: 0.05
idr_relaxed_pvalue: 0.1
idr_top_n_peaks: 150000
keep_chroms_regex: "^chr([1-9]|1[0-9]|2[0-2]|X)$"

# ── QC pipeline ──
gtf: "ref/gencode.v36.annotation.gtf"
promoter_bed: "ref/promoter_chr1-22X.bed"
enhancer_bed: "ref/enhancer_chr1-22X.bed"
spikein_pct_min: 2                        # Active Motif spike-in target range (%)
spikein_pct_max: 10
```
See `ref/config.yaml` for the complete, commented file. To switch spike-in species,
change `spikein_fasta` (the combined index is rebuilt automatically).

## Data Preparation

### Input Files
Please pay attention to the **common suffix** of the fastq related raw files (_R1/2_001.fastq.gz)
Accepted file format should look like: **sample_id**+**common suffix** (e.g. GSF4007-Control_1_S11_R1_001.fastq.gz)

Place paired-end FASTQ files in the `data/` directory following this naming convention:
```
data/{sample}_R1_001.fastq.gz
data/{sample}_R2_001.fastq.gz
```

### Sample Information

Create `ref/samples.csv`:
```csv
sample_id,type,group
GSF4007-Control_1_S11,Control,group1
GSF4007-Control_2_S13,Control,group1
GSF4007-Control_3_S15,Control,group1
GSF4007-NICD3-V5_1_S12,NICD3,group2
GSF4007-NICD3-V5_2_S14,NICD3,group2
GSF4007-NICD3-V5_3_S16,NICD3,group2
```

The `group` column defines conditions / replicate sets. It drives **consensus-peak
reproducibility** (≥2-of-N majority vote for ≥3-replicate conditions, or IDR for
2-replicate conditions) and pairwise IDR in the QC pipeline.

## Running the Pipeline

### Dry Run

To check the workflow without executing any commands:
```bash
# Check primary processing pipeline
snakemake -s snakefile_ATACseq -n

# Check QC pipeline
snakemake -s snakefile_ATAC_QC -n
```

### Local Execution

```bash
# Run primary processing pipeline
snakemake -s snakefile_ATACseq --use-conda --cores 20

# Run QC pipeline (after primary processing completes)
snakemake -s snakefile_ATAC_QC --use-conda --cores 20
```

### Differential Analysis
`ATACseq_Dx.ipynb` is an **R** notebook (`ir` kernel); run it in the **`atacseq-diffbind`**
environment. Interactively:
```bash
jupyter notebook ATACseq_Dx.ipynb
```
Or execute it headless (writes `results/diff_region/` + `results/browser_tracks/`):
```bash
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=ir ATACseq_Dx.ipynb
```
No local R/Bioconductor install? Run it in the **`atacseq-diffbind`** container (pull it
from [Container Execution](#container-execution-docker--apptainer)); Apptainer auto-mounts
the current directory:
```bash
apptainer exec atacseq-diffbind.sif \
    jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=ir ATACseq_Dx.ipynb
```

### Cluster Execution

For execution on a SLURM cluster: (Not tested)
```bash
snakemake -s snakefile_ATACseq --use-conda \
  --cluster "sbatch -p {params.partition} -c {threads} -t {params.time}" \
  --jobs 100
```

### Container Execution (Docker / Apptainer)

Two prebuilt images cover the whole workflow — you install nothing except Docker or
Apptainer:

| Image | Contents | Used for |
|---|---|---|
| **`gynecoloji/atacseq-spikein`** | Snakemake + all five per-rule conda envs | primary pipeline + QC pipeline |
| **`gynecoloji/atacseq-diffbind`** | R / Bioconductor (DESeq2, ChIPseeker, Gviz) + `ir` Jupyter kernel | `ATACseq_Dx.ipynb` differential analysis |

**Download** — pull with Docker, or convert to a local `.sif` once for Apptainer /
Singularity (HPC):

```bash
# Docker
docker pull gynecoloji/atacseq-spikein:latest
docker pull gynecoloji/atacseq-diffbind:latest

# Apptainer / Singularity  (writes ./atacseq-*.sif in the current directory)
apptainer pull atacseq-spikein.sif  docker://gynecoloji/atacseq-spikein:latest
apptainer pull atacseq-diffbind.sif docker://gynecoloji/atacseq-diffbind:latest
```

Genomes/FASTQs are **not** baked into the images; you mount your project directory at run
time (see [`DOCKER.md`](DOCKER.md) for the exact `ref/` and `data/` files the container
expects). Run the primary pipeline first, then QC; the differential notebook runs in the
diffbind image (see [Differential Analysis](#differential-analysis)).

The image's entrypoint is
`snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda`, so anything
you pass after the image name goes straight to `snakemake` (e.g. `-s <snakefile> --cores N`,
or `-n` for a dry run).

#### Docker

```bash
# Pull the published image (or run `docker compose build` to build it locally)
docker pull gynecoloji/atacseq-spikein:latest

# Run from your project directory (which holds snakefile_*, ref/, data/):
# 1) primary processing pipeline
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    gynecoloji/atacseq-spikein:latest -s snakefile_ATACseq --cores 16

# 2) QC pipeline (after the primary run finishes)
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    gynecoloji/atacseq-spikein:latest -s snakefile_ATAC_QC --cores 16

# dry run: append  -n
```

Convenience wrappers `docker compose` and `./run_pipeline.sh` are also provided:

```bash
./run_pipeline.sh snakefile_ATACseq --cores 16
docker compose run --rm atacseq -s snakefile_ATAC_QC --cores 16
```

#### Apptainer / Singularity (HPC)

On clusters without Docker, convert the image to a SIF once and run it with Apptainer.
Apptainer auto-mounts `$HOME`, `/tmp`, and the current directory, and runs as you (no
`--user` needed):

```bash
# One-time: build a local .sif from the Docker Hub image
apptainer pull atacseq-spikein.sif docker://gynecoloji/atacseq-spikein:latest

# Run from your project directory:
# 1) primary processing pipeline
apptainer exec atacseq-spikein.sif \
    snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda \
    -s snakefile_ATACseq --cores 16

# 2) QC pipeline (after the primary run)
apptainer exec atacseq-spikein.sif \
    snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda \
    -s snakefile_ATAC_QC --cores 16
```

Notes:
- **References/data outside the project dir:** if `ref/` genomes live elsewhere (e.g. on
  scratch), bind them in — `--bind /scratch/genomes:/scratch/genomes` — and point
  `ref/config.yaml` at the bound paths.
- **Pre-built envs:** the five conda environments are baked at `/opt/wf-conda`
  (read-only in the SIF) and reused via `--conda-prefix`. If Apptainer reports a
  read-only error writing there, add `--writable-tmpfs` to the `apptainer exec` command.
- `apptainer run atacseq-spikein.sif -s snakefile_ATACseq --cores 16` also works — it
  invokes the same entrypoint.

## Pipeline Details

### 1. Quality Control and Preprocessing

- **FastQC** - Quality assessment of raw reads
- **Fastp** - Adapter trimming and quality filtering with the following parameters:
  - Minimum read length: 30bp
  - Adapter handling: **auto-detects** adapters for paired-end reads by default
    (`--detect_adapter_for_pe`); set `adapter_r1`/`adapter_r2` in `ref/config.yaml`
    to override with explicit sequences
  - Polyg tail trimming
  - Quality trimming: sliding window of 4 with mean quality 20

### 2. Combined-genome Alignment and Read Splitting

- **Combined index** (`build_combined_genome`) - The spike-in FASTA's chromosome names are
  prefixed (e.g. `>chr2L` → `>spikein_chr2L`), the human genome is subset to `align_chroms`,
  the two are concatenated, and a single Bowtie2 index is built.
- **Bowtie2** - One alignment pass to the combined genome:
  - `-X 3000 -I 0 --no-discordant --no-mixed`
- **Read splitting** (`samtools view` by RNAME prefix):
  - **Human reads** → properly-paired (`-f 2 -F 2316`), uniquely-mapped (`grep -v XS:i:`),
    filtering-orphaned mates removed (`process_sam.py`)
  - **Spike-in reads** (`^spikein_`) → same properly-paired + uniquely-mapped +
    orphan-removed (`process_sam.py`) filter as the human reads, then deduplicated
    and counted for normalization
- **Mitochondrial-% QC** - `samtools idxstats` on the human BAM (incl. chrM) is recorded,
  then reads are restricted to `keep_chroms` (chr1–22, chrX; chrM/chrY/non-primary dropped).

### 3. Post-processing

- **Picard MarkDuplicates** - PCR duplicate removal with `REMOVE_DUPLICATES=true`
- **Fragment-level blacklist filtering** - ENCODE blacklist region removal using bedtools with proper paired-end handling

### 4. Peak Calling

- **MACS2** - Peak calling in paired-end mode:
  - Format: BAMPE, genome `hs`, `--nomodel`, q-value cutoff 0.05

### 5. Spike-in Normalization and Consensus Peaks

- **Spike-in normalization** - `compute_spikein_factors` sets `NF = min(spike-in reads) / spike-in reads`
  (the sample with the fewest spike-in reads = 1.0). `create_spikein_bigwig` applies `NF` as a
  deepTools `--scaleFactor`; `create_bigwig` also emits a depth-normalized (RPGC) track for
  before/after comparison.
- **Consensus peaks** (`consensus_peaks`) - Corces-2018 fixed-width (`consensus_window`),
  score-per-million (SPM) iterative-overlap peaks. Per-condition reproducibility is a
  ≥2-of-N majority vote (≥3 replicates) or IDR (2 replicates, via `relaxed_peaks` +
  `reproducible_idr`), unioned across conditions, with blacklist / chrM / chrY / non-primary
  excluded.
- **Fragment counting** (`count_fragments_consensus`) - `featureCounts` (paired-end) over the
  consensus set → `results/consensus/consensus_counts.txt` (a regions × samples matrix).

### 6. Comprehensive QC Metrics (`snakefile_ATAC_QC`)

- **Fragment Size Analysis** - Nucleosomal pattern assessment using `bamPEFragmentSize`
- **TSS Enrichment** - Heatmap/profile plus a **numeric enrichment score** per sample
- **Fingerprint Analysis** - Signal-to-noise using `plotFingerprint`
- **Sample Correlation** - `multiBamSummary` → correlation heatmap/scatter and PCA
- **GC Bias Assessment** - Sequence composition bias (`computeGCBias`)
- **Library Complexity**: **NRF** (Nd/Total), **PBC1** (N1/Nd), **PBC2** (N1/N2)
- **FRiP Score** - Fraction of reads in peaks
- **Spike-in QC** - Spike-in % per sample vs the 2–10% target
- **Peak + annotation summary** - Peak counts/widths; reads in promoters vs enhancers
- **IDR Analysis** - IDR on relaxed peak calls, per replicate pair
- **Interactive QC report** - `atacseq_qc_report.html` (self-contained; all QC except FastQC) plus `multiqc_fastqc.html` (FastQC only)
  (mitochondrial % comes from the primary pipeline's `idxstats`)

### 7. Differential Binding Analysis (`ATACseq_Dx.ipynb`)

An R / Bioconductor notebook (details in [Component 3](#3-differential-binding-analysis-atacseq_dxipynb)):
- **DESeq2** differential binding **NICD3 vs Ctrl** on the consensus matrix — **promoter vs distal**,
  **paired** design, under **spike-in and default (median-of-ratios) normalization** (compared)
- Sample **PCA**, MA/volcano plots, ChIPseeker nearest-gene annotation → `results/diff_region/`
- **Gviz** genome-browser tracks at Notch target loci (spike-in + RPGC, two figures per gene) →
  `results/browser_tracks/`

## Output Files

The pipeline generates the following output directories:

```
results/
├── fastqc/                 # FastQC reports
├── fastp/                  # Trimmed reads and reports
├── aligned/                # Bowtie2 log ({sample}.bowtie2.log); combined SAM is temporary
├── filtered/               # Human filtered BAM + {sample}.idxstats.txt (mito QC) + summaries
├── dedup/                  # Deduplicated BAM + Picard metrics
├── blacklist_filtered/     # Analysis-ready BAM ({sample}.nobl.bam)
├── peaks/                  # MACS2 peak calls (*_peaks.narrowPeak)
├── spikein/
│   ├── aligned/            # Spike-in BAM ({sample}.spikein.bam) + dedup metrics
│   ├── counts/             # {sample}.spikein_count.txt + flagstat
│   └── normalization_factors.tsv   # sample, spikein_reads, norm_factor
├── bigwig/                 # Depth-normalized (RPGC) bigWigs ({sample}.bw)
├── spikein_bigwig/         # Spike-in-scaled bigWigs ({sample}.spikein.bw)
├── consensus/
│   ├── consensus_peaks.bed # Fixed-width, non-overlapping consensus set
│   ├── consensus_peaks.saf # featureCounts input
│   └── consensus_counts.txt # Fragment count matrix (regions × samples)
├── peaks_relaxed/, qc_relaxed_peaks/  # Relaxed MACS2 calls (IDR input; 2-rep conditions / QC)
├── deeptools/              # Fragment size (+ raw lengths), fingerprint, correlation/PCA (+ matrix),
│                           #   GC bias, TSS matrix/plots + downsampled-heatmap JSON
├── bedgraph/               # RPGC bedGraphs (QC)
├── FRiP/                   # FRiP scores (*.frip.txt)
├── library_complexity/     # NRF / PBC1 / PBC2 (*_complexity.txt)
├── idr/                    # IDR peak calls between replicate pairs
├── spikein_qc/             # Spike-in % table (vs 2–10% target)
├── peak_annotation/        # Reads in promoters / enhancers
├── qc/
│   ├── atacseq_qc_report.html      # Interactive QC report (all QC except FastQC)
│   ├── multiqc_fastqc.html         # FastQC-only MultiQC
│   ├── tss_enrichment_scores.tsv, peak_summary.tsv
│   └── blacklist_filtering_stats.txt
# ── Differential-analysis notebook (ATACseq_Dx.ipynb) outputs ──
├── diff_region/            # DESeq2 DB (promoter/distal · spike-in & default), PCA, MA/volcano, annotation
└── browser_tracks/         # Gviz tracks at Notch loci (spike-in + RPGC, PNG + PDF)
```

### Directory Structure
```
ATAC-seq-Pipeline/
├── data/                   # Raw FASTQ files
├── ref/                    # Reference files, config, helper scripts
│   └── COMBINED/           # Combined human+spike-in Bowtie2 index (built by the pipeline)
├── envs/                   # Per-rule conda environment files
├── snakefile_ATACseq, snakefile_ATAC_QC, create_envs.smk
├── Dockerfile, docker-compose.yml, run_pipeline.sh, DOCKER.md
├── results/                # All pipeline outputs (detailed above)
│   └── tmp/                # Temporary processing files
└── logs/                   # Per-rule logs (see below)
```

## Troubleshooting

### Common Issues

1. **Low alignment rate** 
   - Check reference genome compatibility and version
   - Verify adapter trimming parameters in fastp step
   - Check for sample contamination or incorrect library prep
   - Examine FastQC reports for quality issues

2. **High duplication rate**
   - Indicates low library complexity issue
   - Consider optimizing chromatin extraction or Tn5 tagmentation conditions
   - May require increased sequencing depth for better coverage
   - Check PBC1 and PBC2 metrics for severity assessment

3. **Poor TSS enrichment**
   - Issues with sample quality or chromatin accessibility
   - Check protocol for nuclei isolation and tagmentation efficiency
   - Verify blacklist filtering is not removing genuine signal
   - Consider cell type-specific accessibility patterns

4. **High mitochondrial content**
   - Common in certain cell types but may indicate technical issues
   - Consider optimizing nuclear isolation protocols
   - May require additional sequencing depth to compensate for lost reads
   - Check for cytoplasmic contamination during sample prep

5. **Failed IDR analysis**
   - Low reproducibility between biological replicates
   - Check experimental conditions for consistency across replicates
   - May indicate technical issues with sample preparation or processing
   - Consider peak calling parameters if IDR consistently fails

6. **Fragment size distribution issues**
   - Should show nucleosomal ladder pattern (147bp, 294bp, etc.)
   - Flat distribution may indicate DNA degradation or poor tagmentation
   - Very short fragments may indicate over-tagmentation

### Performance Optimization
- Adjust thread counts based on available resources
- Use SSD storage for temporary files when possible
- Monitor memory usage during peak calling with large datasets
- Consider using `--cluster` mode for large-scale analyses

### Log Files

Comprehensive log files for each step are stored in the `logs/` directory:
```
logs/
# Primary pipeline (snakefile_ATACseq)
├── fastqc/, fastp/, bowtie2/, samtools/, dedup/, blacklist_filter/, macs2/
├── build_combined_genome/, spikein_extract_count/, spikein_factors/
├── bigwig/, spikein_bigwig/
├── relaxed_peaks/, reproducible/, consensus/, consensus_counts/
├── blacklist_stats/
# QC pipeline (snakefile_ATAC_QC)
├── deeptools_bedgraph/, deeptools_fragmentsize/, deeptools_plotfingerprint/
├── deeptools_correlation/, deeptools_gc_bias/, deeptools_tss/, tss_enrichment_score/
├── FRiP/, qc_relaxed_peaks/, idr/, library_complexity/
└── spikein_qc/, reads_in_annotations/, peak_summary/, multiqc_fastqc/
```

## Citation

If you use this pipeline, please cite the relevant tools:
- Snakemake: Köster & Rahmann (2012)
- MACS2: Zhang et al. (2008)
- deepTools: Ramírez et al. (2016)
- Bowtie2: Langmead & Salzberg (2012)
- IDR: Li et al. (2011); featureCounts: Liao et al. (2014); Consensus peaks: Corces et al. (2018)

## Contact

For questions or issues with this pipeline, please refer to the individual tool documentation or submit an issue to the repository.

---

**Note**: This pipeline is optimized for human genome analysis (hg38) but can be adapted for other organisms by updating reference files and parameters.