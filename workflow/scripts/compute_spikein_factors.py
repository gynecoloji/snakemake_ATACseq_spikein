#!/usr/bin/env python3
"""Compute spike-in normalization factors from per-sample spike-in read counts.

NF_i = min(counts) / count_i  (Active Motif ATAC-seq spike-in protocol).
The sample with the fewest spike-in reads gets NF = 1.0; all others get NF < 1.0.
"""
from pathlib import Path


def read_count(path):
    """Read a single integer read count (samtools view -c output)."""
    return int(Path(path).read_text().strip())


def sample_from_path(path):
    """Recover the sample id from a '<sample>.spikein_count.txt' filename."""
    return Path(path).name[: -len(".spikein_count.txt")]


def compute_factors(counts):
    """Map {sample: count} -> {sample: NF}. Raise ValueError on any count <= 0."""
    if not counts:
        raise ValueError("no spike-in counts provided")
    zero = sorted(s for s, c in counts.items() if c <= 0)
    if zero:
        raise ValueError(f"zero spike-in reads for sample(s): {', '.join(zero)}")
    d_min = min(counts.values())
    return {s: d_min / c for s, c in counts.items()}


def write_table(counts, factors, out_path):
    """Write TSV: sample <TAB> spikein_reads <TAB> norm_factor (sorted by sample)."""
    lines = ["sample\tspikein_reads\tnorm_factor"]
    for s in sorted(counts):
        lines.append(f"{s}\t{counts[s]}\t{factors[s]:.6f}")
    Path(out_path).write_text("\n".join(lines) + "\n")


def main(count_files, out_path):
    counts = {sample_from_path(f): read_count(f) for f in count_files}
    factors = compute_factors(counts)
    write_table(counts, factors, out_path)


if "snakemake" in globals():  # pragma: no cover
    main(count_files=list(snakemake.input), out_path=str(snakemake.output))  # noqa: F821
