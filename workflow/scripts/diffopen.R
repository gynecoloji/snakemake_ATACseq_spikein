#!/usr/bin/env Rscript
# ============================================================================
# Differential chromatin openness with a SELECTABLE normalization mode
# ============================================================================
# Runs DESeq2 on the consensus fragment-count matrix, differing only in how the
# per-sample size factors are established:
#
#   --mode none     Standard DESeq2 median-of-ratios over ALL consensus peaks.
#                   Baseline. Assumes most regions do not change; a true
#                   genome-wide shift is absorbed and reported as zero.
#
#   --mode spikein  Size factors from the Drosophila spike-in read counts
#                   (sf_i proportional to spike-in depth, geometric mean 1).
#                   The only mode that can detect a genuine GLOBAL shift, but
#                   it is only as trustworthy as the spike-in itself: check
#                   the factor spread and whether it separates by condition.
#
#   --mode ctcf     Median-of-ratios restricted to constitutive CTCF anchors
#                   (peaks overlapping the CTCF cCRE BED). Spike-in free.
#                   Robust to a noisy spike-in, but like `none` it defines the
#                   global level as invariant, so a true global shift is
#                   undetectable by construction.
#
# Inputs (produced by this workflow):
#   results/consensus/consensus_counts.txt     featureCounts matrix + coords
#   config/samples.csv                         design (sample_id, type, group)
#   results/spikein/normalization_factors.tsv  mode=spikein only
#   ref/GRCh38-cCREs.CTCF-only.bed             mode=ctcf only
#
# Usage:
#   Rscript workflow/scripts/diffopen.R --mode ctcf \
#       --counts  results/consensus/consensus_counts.txt \
#       --samples config/samples.csv \
#       --ctcf    ref/GRCh38-cCREs.CTCF-only.bed \
#       --outdir  results/diffopen/ctcf [--ref-label Control]
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- tiny --key value arg parser -------------------------------------------
parse_args <- function(args) {
  out <- list(mode = "none", `ref-label` = "Control", `min-anchors` = "200")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); out[[key]] <- args[i + 1]; i <- i + 2
  }
  out
}

# ---- IO (kept identical in shape to spikein_anchor_shape.R) ----------------

#' featureCounts table -> list(coords, counts[G x n], samples). Cols named by sample.
read_featurecounts_matrix <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  meta <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
  count_cols <- setdiff(colnames(df), meta)
  samples <- sub("\\.nobl\\.bam$", "", basename(count_cols))
  counts <- as.matrix(df[, count_cols, drop = FALSE]); storage.mode(counts) <- "integer"
  colnames(counts) <- samples; rownames(counts) <- df$Geneid
  list(coords = df[, c("Geneid", "Chr", "Start", "End")], counts = counts, samples = samples)
}

#' spike-in table -> named vector of spike-in read counts (per sample).
read_spikein_reads <- function(path) {
  df <- read.delim(path, stringsAsFactors = FALSE)
  setNames(as.numeric(df$spikein_reads), df$sample)
}

#' samples.csv -> data.frame(sample_id, condition[factor, ref first], pair[factor]).
read_design <- function(path, ref_label) {
  s <- read.csv(path, stringsAsFactors = FALSE)
  cond_raw <- s$type
  lv <- c(ref_label, setdiff(unique(cond_raw), ref_label))     # reference level first
  condition <- factor(cond_raw, levels = lv)
  pair <- factor(sub(".*_([0-9]+)_S[0-9]+$", "\\1", s$sample_id))
  data.frame(sample_id = s$sample_id, condition = condition, pair = pair,
             stringsAsFactors = FALSE)
}

#' Logical over consensus peaks: TRUE where the peak overlaps a CTCF cCRE.
ctcf_overlap <- function(coords, ctcf_path) {
  suppressPackageStartupMessages({
    library(GenomicRanges); library(IRanges)
  })
  cc <- read.delim(ctcf_path, header = FALSE, stringsAsFactors = FALSE)
  ctcf  <- GRanges(cc[[1]], IRanges(cc[[2]] + 1L, cc[[3]]))        # BED 0-based -> 1-based
  peaks <- GRanges(coords$Chr, IRanges(coords$Start, coords$End))
  overlapsAny(peaks, ctcf)
}

# ---- size factors, one function per mode ------------------------------------

#' Spike-in depth -> size factors, centered to geometric mean 1.
#' A sample with more spike-in reads received proportionally more material, so
#' its counts are divided by a proportionally larger factor.
size_factors_spikein <- function(spike) {
  stopifnot(all(spike > 0))
  spike / exp(mean(log(spike)))
}

#' Median-of-ratios restricted to a subset of rows (the CTCF anchors).
#' sf_i = exp( median_g [ log(K_gi + 0.5) - mean_j log(K_gj + 0.5) ] ), geomean 1.
size_factors_ctcf <- function(counts, anchor_idx) {
  lc  <- log(counts[anchor_idx, , drop = FALSE] + 0.5)
  ref <- rowMeans(lc)                                  # per-anchor reference
  sf  <- exp(apply(lc - ref, 2, stats::median))
  sf / exp(mean(log(sf)))
}

# ---- DESeq2 -----------------------------------------------------------------

run_deseq2 <- function(counts, coords, condition, pair, size_factors) {
  coldata <- data.frame(condition = condition, pair = pair, row.names = colnames(counts))
  # Use the paired design only when every pair is seen exactly once per condition.
  paired <- nlevels(droplevels(pair)) > 1 &&
            all(table(coldata$pair, coldata$condition) == 1)
  design <- if (paired) ~pair + condition else ~condition
  message("design: ", deparse(design), if (paired) "  (paired)" else "  (unpaired)")

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design)
  if (is.null(size_factors)) {
    dds <- DESeq2::estimateSizeFactors(dds)            # mode=none: DESeq2 default
  } else {
    DESeq2::sizeFactors(dds) <- size_factors           # mode=spikein / ctcf
  }
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)

  cf  <- grep("^condition_", DESeq2::resultsNames(dds), value = TRUE)[1]
  res <- DESeq2::results(dds, name = cf)
  shr <- tryCatch(DESeq2::lfcShrink(dds, coef = cf, type = "apeglm"),
                  error = function(e) { message("lfcShrink failed; using unshrunk LFC"); res })
  list(
    table = data.frame(coords,
                       baseMean       = res$baseMean,
                       log2FoldChange = shr$log2FoldChange,
                       lfcSE          = shr$lfcSE,
                       stat           = res$stat,
                       pvalue         = res$pvalue,
                       padj           = res$padj,
                       row.names = NULL, check.names = FALSE),
    contrast = cf,
    sf = DESeq2::sizeFactors(dds)
  )
}

# ---- main -------------------------------------------------------------------

main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  mode <- a$mode
  if (!mode %in% c("none", "spikein", "ctcf"))
    stop("--mode must be one of: none, spikein, ctcf (got '", mode, "')")
  for (k in c("counts", "samples", "outdir"))
    if (is.null(a[[k]])) stop("missing required --", k)
  if (mode == "spikein" && is.null(a$spikein)) stop("--mode spikein requires --spikein")
  if (mode == "ctcf"    && is.null(a$ctcf))    stop("--mode ctcf requires --ctcf")
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)

  fc  <- read_featurecounts_matrix(a$counts)
  des <- read_design(a$samples, a$`ref-label`)

  samp <- des$sample_id
  stopifnot(all(samp %in% colnames(fc$counts)))
  counts <- fc$counts[, samp, drop = FALSE]

  if (nlevels(droplevels(des$condition)) < 2)
    stop("need >= 2 conditions in samples.csv 'type' column for a differential test")

  # ---- size factors for the requested mode ----
  n_anchors <- NA_integer_
  sf <- switch(mode,
    none = NULL,
    spikein = {
      spike <- read_spikein_reads(a$spikein)
      stopifnot(all(samp %in% names(spike)))
      size_factors_spikein(spike[samp])
    },
    ctcf = {
      idx <- which(ctcf_overlap(fc$coords, a$ctcf))
      n_anchors <<- length(idx)
      message(sprintf("CTCF-overlapping consensus peaks: %d / %d", length(idx), nrow(counts)))
      if (length(idx) < as.integer(a$`min-anchors`))
        stop(sprintf("only %d CTCF anchors (< %s); refusing to normalize on so few",
                     length(idx), a$`min-anchors`))
      size_factors_ctcf(counts, idx)
    })

  out <- run_deseq2(counts, fc$coords, des$condition, des$pair, sf)
  res <- out$table

  # ---- outputs ----
  utils::write.table(res, file.path(a$outdir, "differential_openness.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  sft <- data.frame(sample = names(out$sf), size_factor = as.numeric(out$sf),
                    row.names = NULL, check.names = FALSE)
  utils::write.table(sft, file.path(a$outdir, "size_factors.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)

  # MA plot
  grDevices::png(file.path(a$outdir, "MA_plot.png"), 1100, 850, res = 130)
  sig <- !is.na(res$padj) & res$padj < 0.05
  plot(log2(res$baseMean + 1), res$log2FoldChange, pch = 16, cex = 0.3, col = "#9aa8a5",
       xlab = "log2 mean normalized count", ylab = "log2 fold change",
       main = sprintf("MA — normalization: %s", mode))
  graphics::points(log2(res$baseMean[sig] + 1), res$log2FoldChange[sig],
                   pch = 16, cex = 0.4, col = "#b26a1b")
  graphics::abline(h = 0, col = "#c0c0c0", lty = 3)
  grDevices::dev.off()

  # ---- summary (incl. the spike-in trustworthiness diagnostic) ----
  spread <- max(out$sf) / min(out$sf)
  summ <- c(
    sprintf("normalization mode        : %s", mode),
    sprintf("contrast                  : %s", out$contrast),
    sprintf("consensus peaks           : %d", nrow(counts)),
    if (mode == "ctcf") sprintf("CTCF anchors used         : %d", n_anchors) else NULL,
    sprintf("size-factor spread (max/min): %.2fx", spread),
    sprintf("median |log2FC|           : %.4f", stats::median(abs(res$log2FoldChange), na.rm = TRUE)),
    sprintf("median log2FC (global tilt): %.4f", stats::median(res$log2FoldChange, na.rm = TRUE)),
    sprintf("differential (padj<0.05)  : %d", sum(res$padj < 0.05, na.rm = TRUE)),
    sprintf("nominal (pvalue<0.05)     : %d", sum(res$pvalue < 0.05, na.rm = TRUE)),
    "",
    "Interpretation: a large size-factor spread that tracks condition means the",
    "normalization is confounded -- compare the modes before trusting any one of them.")
  writeLines(summ, file.path(a$outdir, "run_summary.txt"))
  cat(paste(summ, collapse = "\n"), "\n")
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
