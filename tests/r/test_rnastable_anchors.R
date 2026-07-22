#!/usr/bin/env Rscript
# Run in the r-diffopen conda env: Rscript tests/r/test_rnastable_anchors.R
suppressPackageStartupMessages({ library(GenomicRanges); library(IRanges) })
source("workflow/scripts/diffopen.R")

# 3 transcripts: STABLE1 (chr1 +, TSS 1000), UP1 (chr1 -, TSS 5499),
#                STABLE2 (chr2 +, TSS 8000)
tx <- GRanges(c("chr1","chr1","chr2"),
              IRanges(start = c(1000, 5000, 8000), width = 500),
              strand = c("+","-","+"))
tx$gene_name <- c("STABLE1","UP1","STABLE2")

win <- stable_tss_windows(tx, c("STABLE1","STABLE2"), window = 100)
stopifnot(length(win) == 2)                 # UP1 excluded
stopifnot(all(width(win) == 201))           # +/- 100 around a 1bp TSS

# strand-awareness: UP1 is on the '-' strand, so its TSS is the transcript END
# (5000 + 500 - 1 = 5499), NOT its start. A naive start()-based calc would be wrong.
win_minus <- stable_tss_windows(tx, "UP1", window = 100)
stopifnot(length(win_minus) == 1)
stopifnot(start(win_minus) == 5399, end(win_minus) == 5599)   # 5499 +/- 100

# 4 consensus peaks
coords <- data.frame(
  Geneid = paste0("p", 1:4),
  Chr    = c("chr1","chr1","chr2","chr1"),
  Start  = c(1000, 5450, 8050, 2000),
  End    = c(1050, 5500, 8090, 2100),
  stringsAsFactors = FALSE)
# p1 over STABLE1 window; p3 over STABLE2 window; p2 (over UP1) & p4 not.
promoter_is <- c(TRUE, TRUE, FALSE, TRUE)    # p3 is NOT promoter-class

strict  <- rnastable_anchor_idx(coords, promoter_is, win, promoter_class_required = TRUE)
relaxed <- rnastable_anchor_idx(coords, promoter_is, win, promoter_class_required = FALSE)
stopifnot(identical(strict,  1L))            # only p1: promoter AND over-stable
stopifnot(identical(relaxed, c(1L, 3L)))     # p1 and p3: over-stable regardless of class

cat("test_rnastable_anchors.R OK\n")
