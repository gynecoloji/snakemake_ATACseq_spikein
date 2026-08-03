#!/usr/bin/env python3
"""Compute spike-in normalization factors from per-sample spike-in read counts.

NF_i = min(counts) / count_i  (Active Motif ATAC-seq spike-in protocol).
The sample with the fewest spike-in reads gets NF = 1.0; all others get NF < 1.0.

REFUSES to emit factors when the spike-in is too shallow to support them. Before
this guard the only check was `count <= 0`, so a failed or absent spike-in still
produced a full set of factors with no warning. On a public dataset whose deposited
FASTQs contained no Drosophila reads at all, the workflow computed factors from
counts of 20, 20, 36 and 32 -- a spike-in fraction near 0.0001% against the Active
Motif target of 2-10% -- and the resulting 1.8x "size-factor spread", which is pure
Poisson noise at those counts, propagated into 7,488 significant peaks where the
spike-in-free baseline found 115. A 65-fold false-positive inflation, silently.
"""
DEFAULT_MIN_READS = 100_000

from pathlib import Path


def read_count(path):
    """Read a single integer read count (samtools view -c output)."""
    return int(Path(path).read_text().strip())


def sample_from_path(path):
    """Recover the sample id from a '<sample>.spikein_count.txt' filename."""
    return Path(path).name[: -len(".spikein_count.txt")]


def compute_factors(counts, min_reads=DEFAULT_MIN_READS):
    """Map {sample: count} -> {sample: NF}.

    Raises ValueError if any sample has fewer than `min_reads` spike-in reads --
    below that the factors are dominated by counting noise (at n reads the Poisson
    SD is sqrt(n), so at 20 reads the per-sample error is ~22%).
    """
    if not counts:
        raise ValueError("no spike-in counts provided")
    zero = sorted(s for s, c in counts.items() if c <= 0)
    if zero:
        raise ValueError(f"zero spike-in reads for sample(s): {', '.join(zero)}")

    if min_reads and min_reads > 0:
        thin = sorted((c, s) for s, c in counts.items() if c < min_reads)
        if thin:
            detail = ", ".join(f"{s}={c:,}" for c, s in thin)
            raise ValueError(
                f"spike-in too shallow to normalize on: {detail} "
                f"(minimum is {min_reads:,} reads; set `spikein_min_reads` in "
                f"config/config.yaml to override).\n"
                f"Normalization factors derived from counts this low are dominated "
                f"by Poisson noise and will produce large numbers of false "
                f"differential peaks. Check that the spike-in was added, that the "
                f"reads were not filtered out before deposition, and that "
                f"`spikein_prefix` matches the combined index."
            )

    d_min = min(counts.values())
    return {s: d_min / c for s, c in counts.items()}


def write_table(counts, factors, out_path):
    """Write TSV: sample <TAB> spikein_reads <TAB> norm_factor (sorted by sample)."""
    lines = ["sample\tspikein_reads\tnorm_factor"]
    for s in sorted(counts):
        lines.append(f"{s}\t{counts[s]}\t{factors[s]:.6f}")
    Path(out_path).write_text("\n".join(lines) + "\n")


def main(count_files, out_path, min_reads=DEFAULT_MIN_READS):
    counts = {sample_from_path(f): read_count(f) for f in count_files}
    factors = compute_factors(counts, min_reads=min_reads)
    write_table(counts, factors, out_path)


if "snakemake" in globals():  # pragma: no cover
    main(
        count_files=list(snakemake.input),  # noqa: F821
        out_path=str(snakemake.output),  # noqa: F821
        min_reads=int(snakemake.params.get("min_reads", DEFAULT_MIN_READS)),  # noqa: F821
    )
