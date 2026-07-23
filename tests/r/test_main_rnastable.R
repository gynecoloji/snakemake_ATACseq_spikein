#!/usr/bin/env Rscript
# End-to-end CLI smoke test (error path). Run in the r-diffopen conda env:
#   Rscript tests/r/test_main_rnastable.R
suppressPackageStartupMessages({ library(GenomicRanges); library(IRanges) })

d <- tempfile("rnastable_fixture_"); dir.create(d)
on.exit(unlink(d, recursive = TRUE), add = TRUE)

# --- featureCounts matrix: 6 peaks x 4 samples ---
fc <- data.frame(
  Geneid = paste0("pk", 1:6), Chr = "chr1",
  Start = c(1000, 2000, 3000, 4000, 5000, 6000),
  End   = c(1200, 2200, 3200, 4200, 5200, 6200),
  Strand = "+", Length = 201,
  `c1_1_S1.nobl.bam` = 100:105, `c2_2_S2.nobl.bam` = 110:115,
  `t1_3_S3.nobl.bam` = 90:95,   `t2_4_S4.nobl.bam` = 95:100,
  check.names = FALSE, stringsAsFactors = FALSE)
write.table(fc, file.path(d, "counts.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- samples sheet: 2 Control, 2 Treat ---
samples <- data.frame(
  sample_id = c("c1_1_S1","c2_2_S2","t1_3_S3","t2_4_S4"),
  type = c("Control","Control","Treat","Treat"),
  group = c("Control","Control","Treat","Treat"), stringsAsFactors = FALSE)
write.csv(samples, file.path(d, "samples.csv"), row.names = FALSE)

# --- DE table + gene models + BEDs ---
write.table(data.frame(gene = "GENEA", baseMean = 100, log2FoldChange = 0, padj = 0.9),
            file.path(d, "de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tx <- GRanges("chr1", IRanges(start = 1050, width = 1), strand = "+")
tx$gene_name <- "GENEA"; tx$gene_type <- "protein_coding"
saveRDS(list(tx = tx, ex = GRanges()), file.path(d, "models.rds"))
writeLines("chr1\t900\t1300", file.path(d, "prom.bed"))
writeLines("chr1\t9000\t9300", file.path(d, "enh.bed"))

# min-anchors is impossibly high -> must stop() at the floor after full dispatch.
rscript <- shQuote(file.path(R.home("bin"), "Rscript"))
cmd <- sprintf(paste(
  "%s workflow/scripts/diffopen.R --mode rnastable --counts %s --samples %s",
  "--outdir %s --ref-label Control --promoter-bed %s --enhancer-bed %s",
  "--rna-table %s --models %s --tss-window 2000 --min-anchors 1000000",
  "--promoter-class-required true --trim-k 2.5 --trim-iter 2 2>&1"),
  rscript,
  file.path(d, "counts.txt"), file.path(d, "samples.csv"), file.path(d, "out"),
  file.path(d, "prom.bed"), file.path(d, "enh.bed"),
  file.path(d, "de.tsv"), file.path(d, "models.rds"))
out    <- suppressWarnings(system(cmd, intern = TRUE))
joined <- paste(out, collapse = "\n")

# Must fail (non-zero exit) AND fail specifically at the min-anchor FLOOR.
# Asserting on the floor message (not just exit code) is what makes this a real
# red/green: before implementation the failure is "--mode must be one of ...",
# which does NOT contain "rnastable anchors", so the test correctly fails.
stopifnot(!is.null(attr(out, "status")) && attr(out, "status") != 0)
stopifnot(grepl("rnastable anchors", joined))

cat("test_main_rnastable.R OK\n")
