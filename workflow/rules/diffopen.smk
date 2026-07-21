# ── Differential chromatin openness (OPT-IN stage) ──────────────────────
# Deliberately NOT part of `rule all`: a differential test needs at least two
# conditions in config/samples.csv (ideally with replicates), which is not true
# of every dataset this workflow is deployed on. Request it explicitly:
#
#   snakemake --use-conda --cores N diffopen_all          # none + spikein + ctcf
#   snakemake --use-conda --cores N diffopen_anchor_shape # hybrid (Method 6)
#
# All modes consume the primary stage's consensus count matrix, so the primary
# pipeline runs first automatically.
#
# Shared config, directory constants and _diffopen_extra_input() are in common.smk.


# Three selectable normalizations, one output directory each. See
# workflow/scripts/diffopen.R for what each mode does and its trade-off.
rule diffopen:
    wildcard_constraints:
        mode = "none|spikein|ctcf"
    input:
        unpack(_diffopen_extra_input),        # positional: must precede keywords
        # declared so edits to the script invalidate its outputs
        script   = "workflow/scripts/diffopen.R",
        counts   = f"{CONSENSUS_DIR}/consensus_counts.txt",
        samples  = config["samples_table"],
        promoter = config["promoter_bed"],
        enhancer = config["enhancer_bed"],
    output:
        table    = f"{DIFFOPEN_DIR}/{{mode}}/differential_openness.tsv",
        promoter = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_promoter.tsv",
        enhancer = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_enhancer.tsv",
        # pre-filtered nominal-significance subsets (n=3 -> FDR is very
        # conservative; read the direction balance of these, not their size)
        all_p05  = f"{DIFFOPEN_DIR}/{{mode}}/differential_openness_nominal_p05.tsv",
        all_p01  = f"{DIFFOPEN_DIR}/{{mode}}/differential_openness_nominal_p01.tsv",
        prom_p05 = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_promoter_nominal_p05.tsv",
        prom_p01 = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_promoter_nominal_p01.tsv",
        enh_p05  = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_enhancer_nominal_p05.tsv",
        enh_p01  = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_enhancer_nominal_p01.tsv",
        factors  = f"{DIFFOPEN_DIR}/{{mode}}/size_factors.tsv",
        summary  = f"{DIFFOPEN_DIR}/{{mode}}/run_summary.txt",
        ma       = f"{DIFFOPEN_DIR}/{{mode}}/MA_plot.png",
    params:
        outdir    = lambda w: f"{DIFFOPEN_DIR}/{w.mode}",
        ref_label = config.get("diffopen_ref_label", "Control"),
        trim_k    = config.get("ctcf_trim_k", 2.5),
        trim_iter = config.get("ctcf_trim_iter", 2),
        # mode-specific flag, built from whichever extra input was supplied
        extra     = lambda w, input: (
            f"--spikein {input.spikein}" if w.mode == "spikein"
            else f"--ctcf {input.ctcf}"  if w.mode == "ctcf"
            else ""
        ),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/{mode}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/diffopen
        Rscript workflow/scripts/diffopen.R \
            --mode {wildcards.mode} \
            --counts {input.counts} \
            --samples {input.samples} \
            --outdir {params.outdir} \
            --ref-label '{params.ref_label}' \
            --trim-k {params.trim_k} --trim-iter {params.trim_iter} \
            --promoter-bed {input.promoter} --enhancer-bed {input.enhancer} \
            {params.extra} > {log} 2>&1
        """


# Hybrid "anchor + shape" (Method 6): LEVEL from the spike-in, intensity-dependent
# SHAPE from constitutive CTCF anchors (anchors that move between conditions are
# trimmed), combined into per-region DESeq2 normalizationFactors.
rule diffopen_anchor_shape:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/spikein_anchor_shape.R",
        counts  = f"{CONSENSUS_DIR}/consensus_counts.txt",
        spikein = f"{SPIKEIN_DIR}/normalization_factors.tsv",
        samples = config["samples_table"],
        ctcf    = config.get("ctcf_bed", "ref/constitutive_ctcf_hg38.bed"),
        promoter_bed = config["promoter_bed"],
        enhancer_bed = config["enhancer_bed"],
    output:
        table   = f"{DIFFOPEN_DIR}/anchor_shape/differential_openness.tsv",
        # same promoter/enhancer split as the other modes, so the hybrid feeds
        # the identical downstream annotate/enrich/tracks rules
        promoter = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_promoter.tsv",
        enhancer = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_enhancer.tsv",
        all_p05  = f"{DIFFOPEN_DIR}/anchor_shape/differential_openness_nominal_p05.tsv",
        all_p01  = f"{DIFFOPEN_DIR}/anchor_shape/differential_openness_nominal_p01.tsv",
        prom_p05 = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_promoter_nominal_p05.tsv",
        prom_p01 = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_promoter_nominal_p01.tsv",
        enh_p05  = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_enhancer_nominal_p05.tsv",
        enh_p01  = f"{DIFFOPEN_DIR}/anchor_shape/diffopen_enhancer_nominal_p01.tsv",
        anchors = f"{DIFFOPEN_DIR}/anchor_shape/invariant_ctcf_anchors.txt",
        level   = f"{DIFFOPEN_DIR}/anchor_shape/spikein_level.tsv",
        summary = f"{DIFFOPEN_DIR}/anchor_shape/run_summary.txt",
        ma      = f"{DIFFOPEN_DIR}/anchor_shape/anchored_MA.png",
        ecdf    = f"{DIFFOPEN_DIR}/anchor_shape/absolute_ecdf.png",
        shape   = f"{DIFFOPEN_DIR}/anchor_shape/shape_curves.png",
    params:
        outdir    = f"{DIFFOPEN_DIR}/anchor_shape",
        ref_label = config.get("diffopen_ref_label", "Control"),
        span      = config.get("anchor_shape_span", 0.6),
        trim_k    = config.get("anchor_shape_trim_k", 2.5),
        iter      = config.get("anchor_shape_iter", 2),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/anchor_shape.log"
    shell:
        """
        mkdir -p {params.outdir} logs/diffopen
        Rscript workflow/scripts/spikein_anchor_shape.R \
            --counts {input.counts} \
            --spikein {input.spikein} \
            --samples {input.samples} \
            --ctcf {input.ctcf} \
            --outdir {params.outdir} \
            --promoter-bed {input.promoter_bed} \
            --enhancer-bed {input.enhancer_bed} \
            --span {params.span} \
            --trim-k {params.trim_k} \
            --iter {params.iter} \
            --ref-label '{params.ref_label}' > {log} 2>&1
        """


# ── Downstream annotation / enrichment / tracks (per normalization mode) ──
# Nearest TRANSCRIPT TSS (not gene-level 5' ends -- those misassign genes whose
# span is long, e.g. PGK1 by 193 kb). Also caches the parsed GENCODE models as
# RDS so the Gviz rule does not re-read a 1.3 GB GTF.
# Parse the 1.3 GB GENCODE GTF ONCE into a compact RDS shared by every mode.
# Its own wildcard-free rule: as an output of the per-mode rule, three jobs
# would race to write the same file.
rule diffopen_gene_models:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/diffopen_annotate.R",
        gtf = config["gtf"],
    output:
        models = f"{DIFFOPEN_DIR}/gene_models.rds",
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/gene_models.log",
    shell:
        """
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_annotate.R --models-only \
            --gtf {input.gtf} --models {output.models} > {log} 2>&1
        """


rule diffopen_annotate:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/diffopen_annotate.R",
        promoter = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_promoter.tsv",
        enhancer = f"{DIFFOPEN_DIR}/{{mode}}/diffopen_enhancer.tsv",
        gtf      = config["gtf"],
        models   = f"{DIFFOPEN_DIR}/gene_models.rds",
    output:
        summary  = f"{DIFFOPEN_DIR}/{{mode}}/genes/annotation_summary.tsv",
        universe = f"{DIFFOPEN_DIR}/{{mode}}/genes/universe_genes.txt",
    params:
        indir  = lambda w: f"{DIFFOPEN_DIR}/{w.mode}",
        outdir = lambda w: f"{DIFFOPEN_DIR}/{w.mode}/genes",
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/annotate_{mode}.log"
    shell:
        """
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_annotate.R \
            --indir {params.indir} --gtf {input.gtf} \
            --outdir {params.outdir} --models {input.models} > {log} 2>&1
        """


# Offline GO enrichment (clusterProfiler + org.Hs.eg.db -- no network call).
# Gated: sets with <= diffopen_min_genes genes are skipped, so the padj<0.05
# tier is normally not tested at small n. Universe = coding genes reachable
# from any tested peak, NOT the whole genome (which would inflate significance).
rule diffopen_enrich:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/diffopen_enrich.R",
        universe = f"{DIFFOPEN_DIR}/{{mode}}/genes/universe_genes.txt",
    output:
        summary = f"{DIFFOPEN_DIR}/{{mode}}/enrichment/enrichment_summary.tsv",
    params:
        genedir   = lambda w: f"{DIFFOPEN_DIR}/{w.mode}/genes",
        outdir    = lambda w: f"{DIFFOPEN_DIR}/{w.mode}/enrichment",
        min_genes = config.get("diffopen_min_genes", 10),
        ont       = config.get("diffopen_go_ont", "BP"),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/enrich_{mode}.log"
    shell:
        """
        mkdir -p logs/diffopen
        Rscript workflow/scripts/diffopen_enrich.R \
            --genedir {params.genedir} --outdir {params.outdir} \
            --min-genes {params.min_genes} --ont {params.ont} > {log} 2>&1
        """


# Per-mode bigWigs scaled by that mode's own DESeq2 size factors, so the browser
# tracks are on the same footing as the differential test that selected the
# regions. Built from the SAME .nobl.bam featureCounts counted, so the track and
# the table describe the same data.
#
# NOT RPGC: the size factors already carry the depth correction, and layering
# --normalizeUsing on top would normalize twice. Raw coverage x 1/sf only.
#
# anchor_shape is deliberately absent (wildcard_constraints below). Its
# normalization is a per-region G x n offset matrix, not a per-sample scalar, so
# no single --scaleFactor can express it; its tracks stay on the shared RPGC set.
rule diffopen_bigwig:
    wildcard_constraints:
        mode = "none|spikein|ctcf"
    input:
        bam     = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        bai     = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai",
        factors = f"{DIFFOPEN_DIR}/{{mode}}/size_factors.tsv",
    output:
        bw = f"{DIFFOPEN_DIR}/{{mode}}/bigwig/{{sample}}.bw",
    params:
        bin_size  = config["bin_size"],
        blacklist = config["blacklist"],
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/diffopen/bigwig_{mode}_{sample}.log"
    shell:
        """
        mkdir -p $(dirname {output.bw}) logs/diffopen
        # DESeq2 DIVIDES counts by the size factor; bamCoverage MULTIPLIES by
        # --scaleFactor. The track factor is therefore 1/sf, not sf.
        SF=$(awk -F'\t' -v s="{wildcards.sample}" \
               'NR>1 && $1==s && $2+0>0 {{printf "%.10f", 1/$2}}' {input.factors})
        if [ -z "$SF" ]; then
            echo "no usable size factor for {wildcards.sample} in {input.factors}" >&2
            exit 1
        fi
        echo "scaling {wildcards.sample} ({wildcards.mode}) by 1/sf = $SF" > {log}
        bamCoverage --bam {input.bam} \
            --scaleFactor $SF \
            --binSize {params.bin_size} \
            --numberOfProcessors {threads} \
            --extendReads \
            --blackListFileName {params.blacklist} \
            --outFileName {output.bw} >> {log} 2>&1
        """


# Gviz browser tracks for the top up/down regions, one figure per gene
# (PNG + PDF). Same >min-genes gate as the enrichment.
rule diffopen_tracks:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/diffopen_tracks.R",
        summary = f"{DIFFOPEN_DIR}/{{mode}}/genes/annotation_summary.tsv",
        models  = f"{DIFFOPEN_DIR}/gene_models.rds",
        # mode's own size-factor-scaled tracks; RPGC fallback for anchor_shape
        bigwigs = _diffopen_track_bigwigs,
    output:
        done = touch(f"{DIFFOPEN_DIR}/{{mode}}/tracks/.tracks_done"),
    params:
        genedir   = lambda w: f"{DIFFOPEN_DIR}/{w.mode}/genes",
        outdir    = lambda w: f"{DIFFOPEN_DIR}/{w.mode}/tracks",
        bwdir     = _diffopen_track_bwdir,
        tier      = config.get("diffopen_track_tier", "p01"),
        top       = config.get("diffopen_track_top", 5),
        min_genes = config.get("diffopen_min_genes", 10),
    conda:
        "../envs/r-diffopen.yaml"
    log:
        "logs/diffopen/tracks_{mode}.log"
    shell:
        """
        mkdir -p logs/diffopen {params.outdir}
        Rscript workflow/scripts/diffopen_tracks.R \
            --genedir {params.genedir} --bigwigdir {params.bwdir} \
            --models {input.models} --outdir {params.outdir} \
            --tier {params.tier} --top {params.top} \
            --min-genes {params.min_genes} > {log} 2>&1
        """


# Self-contained HTML summary comparing every normalization side by side.
# Runs in the snakemake env (python + pandas); charts are inline SVG so the page
# needs no external assets and opens anywhere.
rule diffopen_report:
    input:
        # declared so edits to the script invalidate its outputs
        script  = "workflow/scripts/build_diffopen_report.py",
        summaries = expand(f"{DIFFOPEN_DIR}/{{mode}}/run_summary.txt", mode=DIFFOPEN_MODES),
        hybrid    = f"{DIFFOPEN_DIR}/anchor_shape/run_summary.txt",
    output:
        html = f"{DIFFOPEN_DIR}/diffopen_report.html",
    params:
        indir    = DIFFOPEN_DIR,
        contrast = lambda w: f"{config.get('diffopen_ref_label', 'Control')} (reference)",
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/diffopen/report.log"
    shell:
        """
        mkdir -p logs/diffopen
        python workflow/scripts/build_diffopen_report.py \
            --diffopen-dir {params.indir} \
            --out {output.html} > {log} 2>&1
        """


# Aggregate target: run every normalization -- the three selectable modes AND the
# anchor+shape hybrid -- so all four can be compared side by side (see each
# mode's run_summary.txt). Note the hybrid's table carries an extra
# `excess_over_global` column and no `stat`, so the schemas are not identical.
rule diffopen_all:
    input:
        expand(f"{DIFFOPEN_DIR}/{{mode}}/differential_openness.tsv", mode=DIFFOPEN_MODES),
        expand(f"{DIFFOPEN_DIR}/{{mode}}/run_summary.txt", mode=DIFFOPEN_MODES),
        rules.diffopen_anchor_shape.output.table,
        rules.diffopen_anchor_shape.output.summary,
        rules.diffopen_report.output.html,
        # Downstream runs for the hybrid too: it now emits the same
        # diffopen_{promoter,enhancer}.tsv layout, so the wildcard rules apply
        # unchanged (rule `diffopen` is constrained to none|spikein|ctcf, so
        # there is no ambiguity over who produces the anchor_shape tables).
        expand(f"{DIFFOPEN_DIR}/{{mode}}/genes/annotation_summary.tsv",
               mode=DIFFOPEN_MODES + ["anchor_shape"]),
        expand(f"{DIFFOPEN_DIR}/{{mode}}/enrichment/enrichment_summary.tsv",
               mode=DIFFOPEN_MODES + ["anchor_shape"]),
        expand(f"{DIFFOPEN_DIR}/{{mode}}/tracks/.tracks_done",
               mode=DIFFOPEN_MODES + ["anchor_shape"]),
