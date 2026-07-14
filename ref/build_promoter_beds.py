#!/usr/bin/env python3
"""Generate promoter BEDs (TSS +/- N bp) under several TSS definitions.

Usage:
    python ref/build_promoter_beds.py <mode> [--windows 3000,5000]
    modes: transcript | unique | mane | cage | all

TSS definitions
---------------
transcript  Every transcript's 5' end (strand-aware) from gencode v36.
            + strand -> transcript start; - strand -> transcript end.
            ~231k windows on chr1-22,X (one per transcript).  <- legacy set
unique      transcript set collapsed to distinct (chrom, TSS, strand) tuples
            (~206k) -- removes isoform double-counting, keeps every alternative
            promoter. A representative transcript supplies the annotation cols.
mane        One canonical TSS per gene: the MANE Select transcript where the
            gene has one, else the 5'-most transcript. ~one row per gene.
cage        Empirical TSS from FANTOM5 robust CAGE peaks (hg38); the peak's
            representative position (BED thickStart) is the TSS. BED6 output.

Window: BED half-open [max(0, tss0-N), min(chromlen, tss0+N)); width up to 2N,
score column = 2N. Kept chroms: chr1-22, chrX. LF line endings.

GTF-derived modes emit the legacy 10-column layout
  tx_id  chrom  start  end  score  strand  tx_id.ver  transcript_type  gene_name  gene_id
sorted by transcript_id. CAGE emits BED6 (chrom start end peak_id score strand)
sorted by genomic position.
"""
import re, sys, gzip, argparse

GTF  = "ref/gencode.v36.annotation.gtf"
FAI  = "ref/hg38.fa.fai"
CAGE = "ref/FANTOM5_CAGE_peaks_hg38.bed.gz"
KEEP = {f"chr{c}" for c in list(range(1, 23)) + ["X"]}

OUT = {
    "transcript": "ref/Promoter_UCSC_hg38_{N}bp_chr1-22X.bed",
    "unique":     "ref/Promoter_uniqueTSS_hg38_{N}bp_chr1-22X.bed",
    "mane":       "ref/Promoter_MANEcanonical_hg38_{N}bp_chr1-22X.bed",
    "cage":       "ref/Promoter_FANTOM5CAGE_hg38_{N}bp_chr1-22X.bed",
}


def chrom_sizes():
    size = {}
    with open(FAI) as f:
        for ln in f:
            p = ln.split("\t")
            if p[0] in KEEP:
                size[p[0]] = int(p[1])
    return size


def parse_gtf_transcripts():
    rx_tid = re.compile(r'transcript_id "([^"]+)"')
    rx_tt  = re.compile(r'transcript_type "([^"]+)"')
    rx_gn  = re.compile(r'gene_name "([^"]+)"')
    rx_gid = re.compile(r'gene_id "([^"]+)"')
    recs = []
    with open(GTF) as f:
        for line in f:
            if "\ttranscript\t" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "transcript" or c[0] not in KEEP:
                continue
            a = c[8]
            tid_v = rx_tid.search(a).group(1); tid = tid_v.split(".")[0]
            m = rx_tt.search(a);  tt = m.group(1) if m else ""
            m = rx_gn.search(a);  gn = m.group(1) if m else ""
            m = rx_gid.search(a); gid = m.group(1).split(".")[0] if m else ""
            strand = c[6]
            tss0 = (int(c[3]) - 1) if strand == "+" else (int(c[4]) - 1)
            recs.append({"chrom": c[0], "tss0": tss0, "strand": strand,
                         "tid_v": tid_v, "tid": tid, "tt": tt, "gn": gn, "gid": gid,
                         "mane": 'tag "MANE_Select"' in a})
    return recs


def select(mode, recs):
    if mode == "transcript":
        return recs
    if mode == "unique":
        seen, out = set(), []
        for r in recs:
            k = (r["chrom"], r["tss0"], r["strand"])
            if k not in seen:
                seen.add(k); out.append(r)
        return out
    if mode == "mane":
        by_gene = {}
        for r in recs:
            by_gene.setdefault(r["gid"], []).append(r)
        out = []
        for gid, rs in by_gene.items():
            mane = [r for r in rs if r["mane"]]
            if mane:
                out.append(mane[0])
            else:                                   # 5'-most transcript of the gene
                plus = rs[0]["strand"] == "+"
                out.append(min(rs, key=lambda r: r["tss0"]) if plus
                           else max(rs, key=lambda r: r["tss0"]))
        return out
    raise ValueError(mode)


def write_gtf_mode(mode, recs, windows, size):
    rows = select(mode, recs)
    rows.sort(key=lambda r: r["tid"])
    for N in windows:
        path = OUT[mode].format(N=N)
        with open(path, "w", newline="\n") as o:
            for r in rows:
                s = max(0, r["tss0"] - N); e = min(size[r["chrom"]], r["tss0"] + N)
                o.write("\t".join(map(str, [r["tid"], r["chrom"], s, e, 2 * N,
                        r["strand"], r["tid_v"], r["tt"], r["gn"], r["gid"]])) + "\n")
        print(f"  wrote {path}: {len(rows)} rows")


def write_cage(windows, size):
    peaks = []
    with gzip.open(CAGE, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 8 or c[0] not in KEEP:
                continue
            peaks.append((c[0], int(c[6]), c[3], c[4], c[5]))   # chrom, thickStart(TSS0), name, score, strand
    peaks.sort(key=lambda p: (p[0], p[1]))
    for N in windows:
        path = OUT["cage"].format(N=N)
        with open(path, "w", newline="\n") as o:
            for chrom, tss0, name, score, strand in peaks:
                s = max(0, tss0 - N); e = min(size[chrom], tss0 + N)
                o.write("\t".join([chrom, str(s), str(e), name, str(score), strand]) + "\n")
        print(f"  wrote {path}: {len(peaks)} rows")


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["transcript", "unique", "mane", "cage", "all"])
    ap.add_argument("--windows", default="3000,5000")
    a = ap.parse_args(argv)
    windows = [int(x) for x in a.windows.split(",")]
    size = chrom_sizes()
    modes = ["transcript", "unique", "mane", "cage"] if a.mode == "all" else [a.mode]
    recs = parse_gtf_transcripts() if any(m != "cage" for m in modes) else None
    for m in modes:
        print(f"[{m}] TSS definition:")
        if m == "cage":
            write_cage(windows, size)
        else:
            write_gtf_mode(m, recs, windows, size)


if __name__ == "__main__":
    main(sys.argv[1:])
