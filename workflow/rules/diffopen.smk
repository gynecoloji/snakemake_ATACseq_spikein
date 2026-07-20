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
        unpack(_diffopen_extra_input),
        counts  = f"{CONSENSUS_DIR}/consensus_counts.txt",
        samples = config["samples_table"],
    output:
        table   = f"{DIFFOPEN_DIR}/{{mode}}/differential_openness.tsv",
        factors = f"{DIFFOPEN_DIR}/{{mode}}/size_factors.tsv",
        summary = f"{DIFFOPEN_DIR}/{{mode}}/run_summary.txt",
        ma      = f"{DIFFOPEN_DIR}/{{mode}}/MA_plot.png",
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
            {params.extra} > {log} 2>&1
        """


# Hybrid "anchor + shape" (Method 6): LEVEL from the spike-in, intensity-dependent
# SHAPE from constitutive CTCF anchors (anchors that move between conditions are
# trimmed), combined into per-region DESeq2 normalizationFactors.
rule diffopen_anchor_shape:
    input:
        counts  = f"{CONSENSUS_DIR}/consensus_counts.txt",
        spikein = f"{SPIKEIN_DIR}/normalization_factors.tsv",
        samples = config["samples_table"],
        ctcf    = config.get("ctcf_bed", "ref/constitutive_ctcf_hg38.bed"),
    output:
        table   = f"{DIFFOPEN_DIR}/anchor_shape/differential_openness.tsv",
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
            --span {params.span} \
            --trim-k {params.trim_k} \
            --iter {params.iter} \
            --ref-label '{params.ref_label}' > {log} 2>&1
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
