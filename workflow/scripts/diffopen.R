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
  out <- list(mode = "none", `ref-label` = "Control", `min-anchors` = "200",
              `trim-k` = "2.5", `trim-iter` = "2", `min-class-peaks` = "100")
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

#' Logical over consensus peaks: TRUE where the peak overlaps any feature in a BED.
bed_overlap <- function(coords, bed_path) {
  suppressPackageStartupMessages({
    library(GenomicRanges); library(IRanges)
  })
  b <- read.delim(bed_path, header = FALSE, stringsAsFactors = FALSE)
  feat  <- GRanges(b[[1]], IRanges(b[[2]] + 1L, b[[3]]))           # BED 0-based -> 1-based
  peaks <- GRanges(coords$Chr, IRanges(coords$Start, coords$End))
  overlapsAny(peaks, feat)
}

#' Backwards-compatible alias used by the CTCF size-factor path.
ctcf_overlap <- function(coords, ctcf_path) bed_overlap(coords, ctcf_path)

#' Classify each consensus peak as promoter / enhancer / other.
#'
#' PROMOTER PRECEDENCE: a peak overlapping both a promoter and an enhancer
#' feature is called a promoter. Without this, peaks straddling an adjacent
#' promoter+enhancer pair would be double-counted across the two classes.
classify_peaks <- function(coords, promoter_bed, enhancer_bed) {
  prom <- bed_overlap(coords, promoter_bed)
  enh  <- bed_overlap(coords, enhancer_bed) & !prom
  factor(ifelse(prom, "promoter", ifelse(enh, "enhancer", "other")),
         levels = c("promoter", "enhancer", "other"))
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
median_of_ratios <- function(counts, idx) {
  lc  <- log(counts[idx, , drop = FALSE] + 0.5)
  ref <- rowMeans(lc)                                  # per-anchor reference
  sf  <- exp(apply(lc - ref, 2, stats::median))
  sf / exp(mean(log(sf)))
}

#' CTCF size factors with an iterative *invariance* trim.
#'
#' The anchors are constitutive in ENCODE, but a given experiment can still move
#' a few of them (a CTCF site inside a responsive enhancer, a copy-number
#' difference, a blacklist-adjacent artifact) and a handful of large movers can
#' drag a median-of-ratios estimate. So: estimate -> measure how far each anchor
#' actually shifted between the two conditions -> drop anchors more than
#' `trim_k` MADs from the median shift -> re-estimate on the survivors.
#'
#' Caveat: this is self-referential, so it cannot detect a genuine *uniform*
#' genome-wide shift (neither can plain median-of-ratios). Only the `spikein`
#' and `anchor_shape` modes can.
#'
#' @return list(sf, idx = surviving anchors, n_start)
size_factors_ctcf <- function(counts, anchor_idx, condition,
                              trim_k = 2.5, iter = 2, min_anchors = 200) {
  idx <- anchor_idx
  sf  <- median_of_ratios(counts, idx)
  lv  <- levels(droplevels(condition))
  if (length(lv) < 2) return(list(sf = sf, idx = idx, n_start = length(anchor_idx)))

  isA <- condition == lv[1]          # reference level (contrast denominator)
  isB <- condition == lv[2]          # contrast numerator
  if (!any(isA) || !any(isB))
    return(list(sf = sf, idx = idx, n_start = length(anchor_idx)))

  for (it in seq_len(max(0L, iter))) {
    # depth-corrected signal at the current anchors
    ln    <- sweep(log2(counts[idx, , drop = FALSE] + 0.5), 2, log2(sf), "-")
    delta <- rowMeans(ln[, isB, drop = FALSE]) - rowMeans(ln[, isA, drop = FALSE])
    med   <- stats::median(delta)
    s_mad <- stats::mad(delta)
    if (!is.finite(s_mad) || s_mad == 0) break
    keep <- abs(delta - med) <= trim_k * s_mad
    if (sum(keep) < min_anchors || all(keep)) break     # don't over-trim / nothing to do
    idx <- idx[keep]
    sf  <- median_of_ratios(counts, idx)
  }
  list(sf = sf, idx = idx, n_start = length(anchor_idx))
}

# ---- DESeq2 -----------------------------------------------------------------

#' DESeq2 size factors as DESeq2 itself would estimate them (mode=none), computed
#' ONCE on the full matrix so every peak class shares them.
size_factors_deseq2 <- function(counts, condition) {
  cd  <- data.frame(condition = condition, row.names = colnames(counts))
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, cd, ~condition)
  DESeq2::sizeFactors(DESeq2::estimateSizeFactors(dds))
}

#' Fit one peak class with the GLOBAL size factors injected (never re-estimated:
#' size factors are a library-level property, so promoter/enhancer/all must share
#' them). Each class still gets its own dispersion trend and its own within-class
#' FDR, which is the point of splitting.
#'
#' @param idx row indices of the class (NULL = all peaks)
fit_class <- function(counts, coords, condition, pair, size_factors, idx = NULL,
                      label = "all") {
  if (!is.null(idx)) {
    counts <- counts[idx, , drop = FALSE]
    coords <- coords[idx, , drop = FALSE]
  }
  coldata <- data.frame(condition = condition, pair = pair, row.names = colnames(counts))
  # Use the paired design only when every pair is seen exactly once per condition.
  paired <- nlevels(droplevels(pair)) > 1 &&
            all(table(coldata$pair, coldata$condition) == 1)
  design <- if (paired) ~pair + condition else ~condition
  message(sprintf("[%s] %d peaks | design %s%s", label, nrow(counts),
                  deparse(design), if (paired) " (paired)" else " (unpaired)"))

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design)
  DESeq2::sizeFactors(dds) <- size_factors
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
    n        = nrow(counts),
    n_sig    = sum(res$padj < 0.05, na.rm = TRUE),
    n_nom    = sum(res$pvalue < 0.05, na.rm = TRUE),
    med_lfc  = stats::median(shr$log2FoldChange, na.rm = TRUE),
    up_frac  = mean(shr$log2FoldChange[which(res$pvalue < 0.05)] > 0, na.rm = TRUE)
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
  # if/else rather than switch(): `<<-` inside a switch branch would assign to the
  # global env, not this frame, silently leaving the counters at NA.
  n_anchors <- NA_integer_
  n_kept    <- NA_integer_
  if (mode == "none") {
    # Estimated once on the FULL matrix, then shared by every peak class.
    sf <- size_factors_deseq2(counts, des$condition)
  } else if (mode == "spikein") {
    spike <- read_spikein_reads(a$spikein)
    stopifnot(all(samp %in% names(spike)))
    sf <- size_factors_spikein(spike[samp])
  } else {                                            # ctcf
    idx <- which(ctcf_overlap(fc$coords, a$ctcf))
    n_anchors <- length(idx)
    message(sprintf("CTCF-overlapping consensus peaks: %d / %d", n_anchors, nrow(counts)))
    if (n_anchors < as.integer(a$`min-anchors`))
      stop(sprintf("only %d CTCF anchors (< %s); refusing to normalize on so few",
                   n_anchors, a$`min-anchors`))
    fit <- size_factors_ctcf(counts, idx, des$condition,
                             trim_k      = as.numeric(a$`trim-k`),
                             iter        = as.integer(a$`trim-iter`),
                             min_anchors = as.integer(a$`min-anchors`))
    sf     <- fit$sf
    n_kept <- length(fit$idx)
    message(sprintf("anchors kept after invariance trim: %d / %d", n_kept, n_anchors))
  }

  # ---- fit: pooled first, then promoter / enhancer separately ----
  # Each class is fit on its own so it gets its own dispersion trend and its own
  # within-class FDR; the size factors above are shared, never re-estimated.
  fits <- list(all = fit_class(counts, fc$coords, des$condition, des$pair, sf,
                               NULL, "all"))
  res  <- fits$all$table

  cls <- NULL
  if (!is.null(a$`promoter-bed`) && !is.null(a$`enhancer-bed`)) {
    cls <- classify_peaks(fc$coords, a$`promoter-bed`, a$`enhancer-bed`)
    message("peak classes: ",
            paste(sprintf("%s=%d", levels(cls), as.integer(table(cls))), collapse = "  "))
    for (k in c("promoter", "enhancer")) {
      idx <- which(cls == k)
      if (length(idx) < as.integer(a$`min-class-peaks`)) {
        message(sprintf("skipping %s: only %d peaks", k, length(idx))); next
      }
      fits[[k]] <- fit_class(counts, fc$coords, des$condition, des$pair, sf, idx, k)
    }
  }

  # ---- outputs ----
  utils::write.table(res, file.path(a$outdir, "differential_openness.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  for (k in c("promoter", "enhancer")) {
    if (!is.null(fits[[k]]))
      utils::write.table(fits[[k]]$table,
                         file.path(a$outdir, sprintf("diffopen_%s.tsv", k)),
                         sep = "\t", quote = FALSE, row.names = FALSE)
  }
  sft <- data.frame(sample = names(sf), size_factor = as.numeric(sf),
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
  spread <- max(sf) / min(sf)
  class_lines <- character(0)
  if (!is.null(cls)) {
    tb <- table(cls)
    class_lines <- c(
      "",
      sprintf("peak classes (promoter precedence): promoter=%d  enhancer=%d  other=%d",
              tb[["promoter"]], tb[["enhancer"]], tb[["other"]]),
      "  class                n      padj<0.05   p<0.05   median log2FC   %% up (of p<0.05)")
    for (k in c("all", "promoter", "enhancer")) {
      f <- fits[[k]]
      if (is.null(f)) next
      class_lines <- c(class_lines,
        sprintf("  %-10s %8d   %8d %8d %14.4f %14.1f%%",
                k, f$n, f$n_sig, f$n_nom, f$med_lfc, 100 * f$up_frac))
    }
    class_lines <- c(class_lines,
      "  (each class fit separately -> own dispersion trend and within-class FDR;",
      "   size factors are shared, computed once on all peaks)")
  }
  summ <- c(
    sprintf("normalization mode        : %s", mode),
    sprintf("contrast                  : %s", fits$all$contrast),
    sprintf("consensus peaks           : %d", nrow(counts)),
    if (mode == "ctcf")
      sprintf("CTCF anchors              : %d overlapping -> %d kept after invariance trim (%.1f%% dropped)",
              n_anchors, n_kept, 100 * (1 - n_kept / n_anchors)) else NULL,
    sprintf("size-factor spread (max/min): %.2fx", spread),
    sprintf("median |log2FC|           : %.4f", stats::median(abs(res$log2FoldChange), na.rm = TRUE)),
    sprintf("median log2FC (global tilt): %.4f", stats::median(res$log2FoldChange, na.rm = TRUE)),
    sprintf("differential (padj<0.05)  : %d", sum(res$padj < 0.05, na.rm = TRUE)),
    sprintf("nominal (pvalue<0.05)     : %d", sum(res$pvalue < 0.05, na.rm = TRUE)),
    class_lines,
    "",
    "Interpretation: a large size-factor spread that tracks condition means the",
    "normalization is confounded -- compare the modes before trusting any one of them.")
  writeLines(summ, file.path(a$outdir, "run_summary.txt"))
  cat(paste(summ, collapse = "\n"), "\n")
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
