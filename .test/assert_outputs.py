#!/usr/bin/env python3
"""Assert the NUMERIC outputs of the executable test case.

Run after `snakemake -d .test atacseq_all` (and optionally `diffopen_all`).
Exits non-zero on any failed assertion.

Why this exists
---------------
Every ATAC-seq pipeline surveyed asserts only that the pipeline finished or that
an output file exists. nf-core/atacseq's CI header says verbatim: "runs the
pipeline with the minimal test dataset to check that it completes without any
syntax errors". nf-core/cutandrun -- a JOSS-accepted pipeline with swappable
spike-in normalization -- tests each of its five normalization modes by asserting
that one bedGraph exists.

File existence would not have caught any of the defects found in this workflow:
the TSS parser emitting a phantom `bins` sample, spike-in factors computed from
20 reads, or a differential stage silently unreachable. So this checks values.

The expected values are known BY CONSTRUCTION from .test/make_testdata.py, not
recorded from a previous run.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from make_testdata import (  # noqa: E402
    CHROM, N_HUMAN_PAIRS, N_PEAKS, PEAK_WIDTH, SAMPLES, peak_starts,
)


def planted_spikein_reads(frac: float) -> int:
    """Spike-in READS planted for a sample, mirroring make_testdata exactly.

    make_testdata draws `int(N_HUMAN_PAIRS * f / (1 - f))` spike-in PAIRS; each
    pair is two reads, which is what `samtools view -c` counts.
    """
    return 2 * int(N_HUMAN_PAIRS * frac / (1 - frac))

ROOT = Path(__file__).parent
RES = ROOT / "results"

failures: list[str] = []
checks = 0


def check(cond: bool, label: str, detail: str = "") -> None:
    global checks
    checks += 1
    if cond:
        print(f"  PASS  {label}")
    else:
        print(f"  FAIL  {label}" + (f"  [{detail}]" if detail else ""))
        failures.append(label)


def read_tsv(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().strip().splitlines()]


# ── 1. spike-in fractions must match what was planted ───────────────────────
print("\n[1] spike-in recovery")
nf = RES / "spikein" / "normalization_factors.tsv"
check(nf.exists(), "normalization_factors.tsv written")
if nf.exists():
    rows = read_tsv(nf)[1:]
    counts = {r[0]: int(r[1]) for r in rows}
    check(set(counts) == set(SAMPLES), "one row per sample",
          f"got {sorted(counts)}")
    # Recovery is near-exact on this fixture (measured 99.97-99.98%), so this can
    # be a tight quantitative check rather than a floor. It is the assertion that
    # would have caught the GSE174272 failure, where the workflow computed size
    # factors from ~20 spike-in reads and reported nothing amiss.
    for s, (cond, frac) in sorted(SAMPLES.items()):
        planted = planted_spikein_reads(frac)
        got = counts.get(s, 0)
        pct = 100.0 * got / planted
        check(
            0.95 <= got / planted <= 1.0,
            f"{s}: recovered {got:,}/{planted:,} planted spike-in reads ({pct:.2f}%)",
            f"{pct:.2f}% outside 95-100%",
        )
    # The guard added in this workflow requires <=2x within-condition spread;
    # the fixture is built to satisfy it, so assert that too.
    by_cond: dict[str, list[int]] = {}
    for s, (cond, _) in SAMPLES.items():
        by_cond.setdefault(cond, []).append(counts[s])
    for cond, vals in by_cond.items():
        spread = max(vals) / min(vals)
        check(spread <= 2.0, f"{cond}: within-condition spike-in spread {spread:.2f}x <= 2x")

# ── 2. peaks must be found WHERE they were planted ──────────────────────────
# A count alone is weak: 40 peaks in the wrong places would pass. This checks
# that each planted region is actually covered by a called peak.
print("\n[2] peak calling")
planted = peak_starts()
for s in sorted(SAMPLES):
    p = RES / "peaks" / f"{s}_peaks.narrowPeak"
    if not p.exists():
        check(False, f"{s}: narrowPeak written", "missing")
        continue
    rows = [ln.split("\t") for ln in p.read_text().strip().splitlines() if ln]
    n = len(rows)
    # 40 planted; MACS2 splits or merges a few, so allow a modest window.
    check(0.8 * N_PEAKS <= n <= 1.5 * N_PEAKS,
          f"{s}: called {n} peaks (planted {N_PEAKS})", f"n={n}")

    chroms = {r[0] for r in rows}
    check(chroms <= {CHROM}, f"{s}: peaks only on {CHROM}", f"got {chroms}")

    called = [(int(r[1]), int(r[2])) for r in rows if r[0] == CHROM]
    hit = sum(
        any(a < st + PEAK_WIDTH and b > st for a, b in called) for st in planted
    )
    check(hit >= 0.9 * N_PEAKS,
          f"{s}: {hit}/{N_PEAKS} planted regions recovered",
          f"only {hit}")

# ── 3. consensus set and count matrix ───────────────────────────────────────
print("\n[3] consensus + counts")
cons = RES / "consensus" / "consensus_peaks.bed"
cnt = RES / "consensus" / "consensus_counts.txt"
check(cons.exists() and cons.stat().st_size > 0, "consensus_peaks.bed non-empty")
check(cnt.exists() and cnt.stat().st_size > 0, "consensus_counts.txt non-empty")
if cnt.exists() and cnt.stat().st_size:
    lines = [ln for ln in cnt.read_text().splitlines() if not ln.startswith("#")]
    header = lines[0].split("\t")
    n_rows = len(lines) - 1
    meta = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    sample_cols = [c for c in header if c not in meta]
    check(len(sample_cols) == len(SAMPLES),
          f"count matrix has {len(sample_cols)} sample columns (expect {len(SAMPLES)})")
    check(n_rows >= N_PEAKS * 0.5, f"count matrix has {n_rows} regions", f"n={n_rows}")

# ── 4. TSS enrichment table must contain ONLY samples ───────────────────────
# Regression: the parser skipped one deepTools header line where the format has
# two, so a phantom row named `bins` was emitted as if it were a sample.
print("\n[4] TSS enrichment table (regression: phantom 'bins' row)")
tss = RES / "qc" / "tss_enrichment_scores.tsv"
if tss.exists():
    rows = read_tsv(tss)[1:]
    names = {r[0] for r in rows}
    check(names == set(SAMPLES), f"exactly {len(SAMPLES)} sample rows, no extras",
          f"got {sorted(names)}")
    check("bins" not in names, "no phantom 'bins' row")
else:
    print("  SKIP  TSS table absent (qc_all not run)")

# ── 5. differential stage, if it was run ────────────────────────────────────
print("\n[5] differential openness (if diffopen_all was run)")
dd = RES / "diffopen"
if dd.is_dir() and any(dd.iterdir()):
    for mode in sorted(p.name for p in dd.iterdir() if p.is_dir()):
        t = dd / mode / "differential_openness.tsv"
        if not t.exists():
            check(False, f"{mode}: differential_openness.tsv written", "missing")
            continue
        rows = read_tsv(t)
        hdr = rows[0]
        check("padj" in hdr and "log2FoldChange" in hdr,
              f"{mode}: table has padj + log2FoldChange columns")
        check(len(rows) - 1 >= N_PEAKS * 0.5,
              f"{mode}: {len(rows)-1} regions tested")
else:
    print("  SKIP  diffopen not run")

# ── verdict ─────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
if failures:
    print(f"FAILED {len(failures)}/{checks} checks:")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print(f"All {checks} checks passed.")
