#!/usr/bin/env python3
"""Numeric TSS enrichment score per sample from a deepTools plotProfile data table.

TSS enrichment = signal at the TSS (center bin) / mean signal in the flanking
background (outer edge bins). ~1 means no enrichment; higher is better (ENCODE
uses this as a key ATAC signal-to-noise metric).
"""
from pathlib import Path


def enrichment(values, edge_frac=0.1):
    """Center-bin value / mean of the outer `edge_frac` bins on each side."""
    n = len(values)
    if n == 0:
        return 0.0
    center = values[n // 2]
    k = max(1, int(n * edge_frac))
    bg = values[:k] + values[-k:]
    m = sum(bg) / len(bg)
    return center / m if m > 0 else 0.0


# deepTools `plotProfile --outFileNameData` writes TWO header rows before the
# per-sample data:
#
#     bin labels<TAB><TAB>-2.0Kb ...        <- tick labels
#     bins<TAB><TAB>1.0<TAB>2.0<TAB>3.0 ... <- bin indices
#     <sample><TAB>genes<TAB>0.38<TAB>...   <- one row per sample
#
# Skipping only the first (`lines[1:]`) left the `bins` row to be parsed as if it
# were a sample, because its bin indices are perfectly good floats. That put a
# phantom row named `bins` into the shipped QC table and MultiQC bargraph, with a
# meaningless "enrichment" computed from bin numbers.
#
# Filtering by label rather than by position, so a future deepTools release that
# adds or reorders a header row cannot silently reintroduce the same class of bug.
_HEADER_LABELS = {"bin labels", "bins"}


def parse_profile(path):
    """Parse `plotProfile --outFileNameData` into {sample: [profile floats]}.

    Each data row is `<sample> <group> <per-bin means...>`; the row's first field
    is the sample and every field that parses as a float is signal.
    """
    out = {}
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    for line in lines:
        f = line.split("\t")
        if f[0].strip() in _HEADER_LABELS:
            continue
        vals = []
        for x in f[1:]:
            try:
                vals.append(float(x))
            except ValueError:
                pass
        if vals:
            out[f[0]] = vals
    return out


def build(profile_path):
    prof = parse_profile(profile_path)
    return sorted((s, round(enrichment(v), 3)) for s, v in prof.items())


def write_tsv(rows, path):
    lines = ["sample\ttss_enrichment"]
    lines += [f"{s}\t{score}" for s, score in rows]
    Path(path).write_text("\n".join(lines) + "\n")


def write_mqc(rows, path):
    lines = [
        "# id: tss_enrichment",
        "# section_name: 'TSS enrichment'",
        "# description: 'TSS enrichment score (center/background). Higher is better; <5 is poor for ATAC.'",
        "# plot_type: 'bargraph'",
        "# pconfig:",
        "#    id: 'tss_enrichment_plot'",
        "#    title: 'TSS enrichment'",
        "sample\tTSS enrichment",
    ]
    lines += [f"{s}\t{score}" for s, score in rows]
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    _rows = build(str(snakemake.input.profile))  # noqa: F821
    write_tsv(_rows, str(snakemake.output.tsv))  # noqa: F821
    write_mqc(_rows, str(snakemake.output.mqc))  # noqa: F821
