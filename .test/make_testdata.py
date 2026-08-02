#!/usr/bin/env python3
"""Generate the miniature reference and reads for the executable test case.

Everything here is SYNTHETIC and DETERMINISTIC (fixed seed), so the expected
outputs are known by construction rather than by having run the pipeline once and
recorded whatever came out. That is what makes `.test/assert_outputs.py` able to
assert numbers instead of just file existence.

Design
------
Genome      chr21, 400 kb of random sequence, plus a spike-in contig `2L` of
            120 kb. Small enough that bowtie2-build finishes in seconds; large
            enough that MACS2 has a real background to model.

Peaks       40 accessible regions of 400 bp on chr21. Reads pile up there against
            a low uniform background, so peak calling has something to find.

Truth       Two things are planted and later asserted:
              1. SPIKE-IN FRACTION per sample, set exactly (see SAMPLES below).
                 The workflow must recover it to within a tolerance.
              2. A DIFFERENTIAL BLOCK: peaks 0-9 get 4x more reads in the
                 Treatment samples. Peaks 10-39 are unchanged (true nulls).

Everything is written where the workflow expects it:
    .test/ref/{hg38.fa,dm6.fa,...}   .test/data/{sample}_R{1,2}_001.fastq.gz

Usage:  python .test/make_testdata.py [--outdir .test]
"""
from __future__ import annotations

import argparse
import gzip
import random
import zlib
from pathlib import Path

SEED = 20260801
CHROM = "chr21"
CHROM_LEN = 400_000
SPIKE_CHROM = "2L"
SPIKE_LEN = 120_000

N_PEAKS = 40
PEAK_WIDTH = 400
DIFF_PEAKS = 10          # peaks 0..9 change; 10..39 are true nulls
DIFF_FOLD = 4            # treatment enrichment at those peaks

READ_LEN = 50
FRAG_MIN, FRAG_MAX = 120, 300
N_HUMAN_PAIRS = 60_000   # per sample
BACKGROUND_FRAC = 0.30   # fraction of human pairs placed uniformly, not in peaks

# sample -> (condition, spike-in fraction of total pairs)
# Spread deliberately kept tight WITHIN each condition: the workflow now refuses
# to normalize when within-condition spike-in spread exceeds 2x.
SAMPLES = {
    "test_ctrl_1": ("Control", 0.10),
    "test_ctrl_2": ("Control", 0.11),
    "test_trt_1": ("Treatment", 0.10),
    "test_trt_2": ("Treatment", 0.11),
}

BASES = "ACGT"


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_fasta(path: Path, name: str, seq: str, width: int = 60) -> None:
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i : i + width] + "\n")


def peak_starts() -> list[int]:
    """Evenly spaced peaks, clear of the contig ends."""
    margin = 20_000
    span = CHROM_LEN - 2 * margin
    step = span // N_PEAKS
    return [margin + i * step for i in range(N_PEAKS)]


def draw_fragment(rng, seq, lo, hi):
    """Pick a fragment wholly inside [lo, hi); return (start, length)."""
    flen = rng.randint(FRAG_MIN, FRAG_MAX)
    if hi - lo <= flen:
        return None
    return rng.randrange(lo, hi - flen), flen


def read_name(rid: int) -> str:
    """Illumina-style name: instrument:run:flowcell:lane:tile:x:y.

    Picard's optical-duplicate finder parses the last three colon-separated
    fields as integers and warns on every run if they are not numeric, so the
    synthetic names follow the real layout rather than a bare counter.
    """
    tile = 1101 + (rid // 10_000)
    x = rid % 10_000
    y = 1000 + (rid % 7_000)
    return f"SIM:1:TESTFLOWCELL:1:{tile}:{x}:{y}"


def emit_pair(fh1, fh2, rid, seq, start, flen, chrom):
    frag = seq[start : start + flen]
    if len(frag) < READ_LEN:
        return
    r1 = frag[:READ_LEN]
    r2 = revcomp(frag[-READ_LEN:])
    q = "I" * READ_LEN
    name = read_name(rid)
    fh1.write(f"@{name} {chrom}:{start}/1\n{r1}\n+\n{q}\n")
    fh2.write(f"@{name} {chrom}:{start}/2\n{r2}\n+\n{q}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=str(Path(__file__).parent))
    args = ap.parse_args()
    out = Path(args.outdir)
    (out / "ref").mkdir(parents=True, exist_ok=True)
    (out / "data").mkdir(parents=True, exist_ok=True)

    rng = random.Random(SEED)
    human = random_seq(rng, CHROM_LEN)
    spike = random_seq(rng, SPIKE_LEN)

    write_fasta(out / "ref" / "hg38.fa", CHROM, human)
    write_fasta(out / "ref" / "dm6.fa", SPIKE_CHROM, spike)

    starts = peak_starts()

    # --- annotation: promoters over the first half of the peaks, enhancers the rest,
    #     so the promoter/enhancer split has both classes populated.
    with open(out / "ref" / "promoter_chr1-22X.bed", "w") as fh:
        for i, s in enumerate(starts[: N_PEAKS // 2]):
            fh.write(f"{CHROM}\t{s - 200}\t{s + PEAK_WIDTH + 200}\tPROM{i}\n")
    with open(out / "ref" / "enhancer_chr1-22X.bed", "w") as fh:
        for i, s in enumerate(starts[N_PEAKS // 2 :]):
            fh.write(f"{CHROM}\t{s - 200}\t{s + PEAK_WIDTH + 200}\tENH{i}\n")
    # CTCF anchors: a subset of the NON-differential peaks, so the ctcf mode has
    # genuinely invariant anchors to normalize on.
    with open(out / "ref" / "constitutive_ctcf_hg38.bed", "w") as fh:
        for i, s in enumerate(starts[DIFF_PEAKS:], start=DIFF_PEAKS):
            fh.write(f"{CHROM}\t{s}\t{s + PEAK_WIDTH}\tCTCF{i}\n")
    # Blacklist: one region deliberately away from every peak.
    with open(out / "ref" / "hg38_blacklist_regions.bed", "w") as fh:
        fh.write(f"{CHROM}\t1000\t3000\tblacklist_test\n")

    # --- minimal GTF: one transcript per peak, TSS at the peak start.
    with open(out / "ref" / "gencode.v36.annotation.gtf", "w") as fh:
        fh.write("##description: synthetic test annotation\n")
        for i, s in enumerate(starts):
            gid, tid = f"ENSGTEST{i:05d}", f"ENSTTEST{i:05d}"
            attrs = (
                f'gene_id "{gid}"; transcript_id "{tid}"; gene_type "protein_coding"; '
                f'gene_name "TESTGENE{i}"; transcript_type "protein_coding";'
            )
            fh.write(f"{CHROM}\tTEST\tgene\t{s+1}\t{s+PEAK_WIDTH}\t.\t+\t.\t{attrs}\n")
            fh.write(f"{CHROM}\tTEST\ttranscript\t{s+1}\t{s+PEAK_WIDTH}\t.\t+\t.\t{attrs}\n")
            fh.write(f"{CHROM}\tTEST\texon\t{s+1}\t{s+PEAK_WIDTH}\t.\t+\t.\t{attrs}\n")

    # --- reads
    # The per-sample seed offset must be STABLE ACROSS PROCESSES. Python salts
    # str.__hash__ per interpreter unless PYTHONHASHSEED is set, so hash(sample)
    # would silently make this generator non-deterministic -- defeating the whole
    # point of constructing the expected values rather than recording them.
    # zlib.crc32 is stable by definition.
    for sample, (cond, spike_frac) in SAMPLES.items():
        rs = random.Random(SEED + zlib.crc32(sample.encode()) % 10_000)
        n_spike = int(N_HUMAN_PAIRS * spike_frac / (1 - spike_frac))
        n_bg = int(N_HUMAN_PAIRS * BACKGROUND_FRAC)
        n_peak = N_HUMAN_PAIRS - n_bg

        # weight peaks: the differential block is enriched in Treatment
        weights = [
            (DIFF_FOLD if (i < DIFF_PEAKS and cond == "Treatment") else 1)
            for i in range(N_PEAKS)
        ]
        total_w = sum(weights)

        f1 = out / "data" / f"{sample}_R1_001.fastq.gz"
        f2 = out / "data" / f"{sample}_R2_001.fastq.gz"
        rid = 0
        with gzip.open(f1, "wt") as fh1, gzip.open(f2, "wt") as fh2:
            # peak-associated human pairs
            for i, s in enumerate(starts):
                k = int(n_peak * weights[i] / total_w)
                lo, hi = max(0, s - 100), min(CHROM_LEN, s + PEAK_WIDTH + 100)
                for _ in range(k):
                    d = draw_fragment(rs, human, lo, hi)
                    if d:
                        emit_pair(fh1, fh2, rid, human, d[0], d[1], CHROM)
                        rid += 1
            # uniform human background
            for _ in range(n_bg):
                d = draw_fragment(rs, human, 0, CHROM_LEN)
                if d:
                    emit_pair(fh1, fh2, rid, human, d[0], d[1], CHROM)
                    rid += 1
            # spike-in, uniform over the spike-in contig
            for _ in range(n_spike):
                d = draw_fragment(rs, spike, 0, SPIKE_LEN)
                if d:
                    emit_pair(fh1, fh2, rid, spike, d[0], d[1], SPIKE_CHROM)
                    rid += 1
        print(f"{sample}: {rid} pairs  (spike-in target {spike_frac:.0%})")

    print(f"\nwrote reference + reads under {out}")


if __name__ == "__main__":
    main()
