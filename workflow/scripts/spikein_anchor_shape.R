#!/usr/bin/env Rscript
# ============================================================================
# Method 6 — spike-in-anchored, CTCF-shape-corrected differential openness
# ============================================================================
# Hybrid "anchor + shape" normalization for ATAC-seq (see
# atacseq_spikein_Dx/method6_spikein_anchor_shape.md):
#   * LEVEL  o_i     : from the Drosophila spike-in (per-sample scalar).
#   * SHAPE  s_i(A)  : mean-zero, intensity-dependent curve fit from constitutive
#                      CTCF anchors (ENCODE SCREEN CTCF-only cCREs), then applied
#                      genome-wide. CTCF anchors that MOVE between conditions are
#                      trimmed (self-consistency) so only invariant sites shape it.
# The combined per-region offset feeds DESeq2 as normalizationFactors.
#
# Consumes snakemake_ATACseq_spikein outputs:
#   results/consensus/consensus_counts.txt        (featureCounts: K + coords)
#   results/spikein/normalization_factors.tsv     (spike-in reads -> level o_i)
#   config/samples.csv                            (design: condition + pair)
#   ref/GRCh38-cCREs.CTCF-only.bed                (constitutive CTCF anchors)
#
# Usage:
#   Rscript workflow/scripts/spikein_anchor_shape.R \
#       --counts   results/consensus/consensus_counts.txt \
#       --spikein  results/spikein/normalization_factors.tsv \
#       --samples  config/samples.csv \
#       --ctcf     ref/GRCh38-cCREs.CTCF-only.bed \
#       --outdir   results/diffopen/anchor_shape \
#       [--span 0.6] [--trim-k 2.5] [--iter 2] [--ref-label Control]
#
# Or via Snakemake (recommended):  snakemake --use-conda diffopen_anchor_shape
#
# Requires (workflow/envs/r-diffopen.yaml): DESeq2, apeglm, GenomicRanges,
# IRanges.  NOTE: these are NOT in ATACSeq_Dx.yaml, which contains no R at all.
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(GenomicRanges)
  library(IRanges)
})

# ---- tiny --key value arg parser -------------------------------------------
parse_args <- function(args) {
  out <- list(span = "0.6", `trim-k` = "2.5", iter = "2", `ref-label` = "Control",
              `min-refs` = "200")
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i]); out[[key]] <- args[i + 1]; i <- i + 2
  }
  out
}

# ---- IO --------------------------------------------------------------------

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

# ---- CTCF anchor selection --------------------------------------------------

#' Logical over consensus peaks: TRUE where the peak overlaps a CTCF cCRE.
ctcf_overlap <- function(coords, ctcf_path) {
  cc <- read.delim(ctcf_path, header = FALSE, stringsAsFactors = FALSE)
  ctcf  <- GRanges(cc[[1]], IRanges(cc[[2]] + 1L, cc[[3]]))        # BED 0-based -> 1-based
  peaks <- GRanges(coords$Chr, IRanges(coords$Start, coords$End))
  overlapsAny(peaks, ctcf)
}

# ---- core: level + (trimmed) CTCF shape -> normalizationFactors -------------

#' @return list(NF, o, shape, A, ybar, ref_idx_final, g0, ctcf_level_divergence)
anchor_shape_offsets <- function(counts, spike_reads, ref_idx0, condition,
                                 span = 0.6, trim_k = 2.5, iter = 2, min_refs = 200) {
  n <- ncol(counts); G <- nrow(counts)
  # LEVEL from spike-in: o_i estimates (const - c_i); reference const is arbitrary
  # (DESeq2 row-centering removes it). median() picks a symmetric reference sample.
  o <- log2(stats::median(spike_reads) / spike_reads)

  # absolute-scale signal and intensity coordinate
  y    <- sweep(log2(counts + 0.5), 2, o, "+")     # y_gi = log2(K)+o_i
  ybar <- rowMeans(y)
  A    <- ybar
  dd   <- stats::density(A)
  wall <- stats::approx(dd$x, dd$y, xout = A, rule = 2)$y   # genome A-density weights
  mean_counts <- rowMeans(counts)
  isT  <- as.integer(condition) == 2               # treatment (non-reference) samples

  fit_shape <- function(ref_idx) {
    shape <- matrix(0, G, n)
    for (i in seq_len(n)) {
      ri  <- y[ref_idx, i] - ybar[ref_idx]
      Ai  <- A[ref_idx]
      wts <- 1 / (1 / (counts[ref_idx, i] + 1) + 1 / (mean_counts[ref_idx] + 1))
      fit <- stats::loess(ri ~ Ai, span = span, family = "symmetric", weights = wts,
                          control = stats::loess.control(surface = "direct"))
      s <- stats::predict(fit, newdata = data.frame(Ai = A)); s[is.na(s)] <- 0
      shape[, i] <- s - sum(s * wall) / sum(wall)   # mean-zero over genome A-density
    }
    shape
  }

  ref_idx <- ref_idx0
  g0 <- NA_real_
  for (it in seq_len(max(1L, iter))) {
    shape <- fit_shape(ref_idx)
    ycorr <- y - shape                              # corrected absolute signal
    delta <- rowMeans(ycorr[ref_idx, isT,  drop = FALSE]) -
             rowMeans(ycorr[ref_idx, !isT, drop = FALSE])
    med <- stats::median(delta); s_mad <- stats::mad(delta)
    g0  <- med                                      # recovered global shift @ anchors
    if (it == iter || s_mad == 0) break
    keep <- abs(delta - med) <= trim_k * s_mad      # trim CTCF anchors that MOVED
    if (sum(keep) < min_refs) break                 # don't over-trim
    ref_idx <- ref_idx[keep]
  }

  # combined normalizationFactors:  log2(s_gi) = shape_i(A_g) - o_i, row geomean 1
  NF <- 2^(sweep(shape, 2, o, "-"))
  NF <- NF / exp(rowMeans(log(NF)))
  list(NF = NF, o = o, shape = shape, A = A, ybar = ybar,
       ref_idx_final = ref_idx, g0 = g0, ycorr = y - shape)
}

# ---- DESeq2 with the region-specific offset --------------------------------

run_deseq2 <- function(counts, coords, condition, pair, NF, g0) {
  coldata <- data.frame(condition = condition, pair = pair, row.names = colnames(counts))
  paired  <- nlevels(droplevels(pair)) > 1 &&
             all(table(coldata$pair, coldata$condition) == 1)
  design  <- if (paired) ~pair + condition else ~condition
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design)
  DESeq2::normalizationFactors(dds) <- NF          # NOT estimateSizeFactors: use ours
  dds <- DESeq2::estimateDispersions(dds)
  dds <- DESeq2::nbinomWaldTest(dds)
  cf  <- grep("^condition_", DESeq2::resultsNames(dds), value = TRUE)[1]
  res <- DESeq2::results(dds, name = cf)
  shr <- tryCatch(DESeq2::lfcShrink(dds, coef = cf, type = "apeglm"),
                  error = function(e) res)
  data.frame(coords,
             baseMean          = res$baseMean,
             log2FoldChange    = shr$log2FoldChange,   # ABSOLUTE (spike-anchored)
             excess_over_global= shr$log2FoldChange - g0,  # local rho effect
             lfcSE             = shr$lfcSE,
             pvalue            = res$pvalue,
             padj              = res$padj,
             row.names = NULL, check.names = FALSE)
}

# ---- diagnostics ------------------------------------------------------------

plot_diagnostics <- function(outdir, fit, res, ref_final_names, coords, condition, ycorr) {
  is_ref <- coords$Geneid %in% ref_final_names
  # 1) Anchored MA: reference line at the recovered global shift g0
  grDevices::png(file.path(outdir, "anchored_MA.png"), 1100, 850, res = 130)
  plot(fit$A, res$log2FoldChange, pch = 16, cex = 0.3, col = "#9aa8a5",
       xlab = "mean absolute intensity  A (log2)", ylab = "log2FC  (absolute, spike-anchored)",
       main = "Anchored MA — line at recovered global shift g0")
  graphics::points(fit$A[is_ref], res$log2FoldChange[is_ref], pch = 16, cex = 0.4, col = "#0b6e7c")
  graphics::abline(h = fit$g0, col = "#b26a1b", lwd = 2, lty = 2)
  graphics::abline(h = 0, col = "#c0c0c0", lwd = 1, lty = 3)
  graphics::legend("topright", bty = "n",
                   legend = c("all peaks", "CTCF anchors (kept)", sprintf("g0 = %.3f", fit$g0)),
                   col = c("#9aa8a5", "#0b6e7c", "#b26a1b"), pch = c(16, 16, NA),
                   lty = c(NA, NA, 2), lwd = c(NA, NA, 2))
  grDevices::dev.off()

  # 2) Absolute-scale ECDF: per-condition mean corrected signal per region
  isT <- as.integer(condition) == 2
  mC  <- rowMeans(ycorr[, !isT, drop = FALSE]); mT <- rowMeans(ycorr[, isT, drop = FALSE])
  grDevices::png(file.path(outdir, "absolute_ecdf.png"), 1000, 800, res = 130)
  plot(stats::ecdf(mC), col = "#0b6e7c", lwd = 2, main = "Absolute-scale accessibility ECDF",
       xlab = "corrected absolute signal (log2)", ylab = "cumulative fraction", do.points = FALSE)
  graphics::lines(stats::ecdf(mT), col = "#b26a1b", lwd = 2)
  graphics::legend("bottomright", bty = "n", legend = levels(condition),
                   col = c("#0b6e7c", "#b26a1b"), lwd = 2)
  grDevices::dev.off()

  # 3) Shape curves per sample
  grDevices::png(file.path(outdir, "shape_curves.png"), 1000, 800, res = 130)
  o <- order(fit$A)
  matplot(fit$A[o], fit$shape[o, ], type = "l", lty = 1,
          xlab = "A (log2 intensity)", ylab = "shape_i(A)  (mean-zero)",
          main = "Intensity-dependent technical shape, per sample")
  graphics::abline(h = 0, col = "#c0c0c0", lty = 3)
  grDevices::dev.off()
}

# ---- main ------------------------------------------------------------------

main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  for (k in c("counts", "spikein", "samples", "ctcf", "outdir"))
    if (is.null(a[[k]])) stop("missing required --", k)
  dir.create(a$outdir, showWarnings = FALSE, recursive = TRUE)

  fc    <- read_featurecounts_matrix(a$counts)
  spike <- read_spikein_reads(a$spikein)
  des   <- read_design(a$samples, a$`ref-label`)

  # align all three on the same sample order
  samp <- des$sample_id
  stopifnot(all(samp %in% colnames(fc$counts)), all(samp %in% names(spike)))
  counts <- fc$counts[, samp, drop = FALSE]
  spike  <- spike[samp]
  condition <- des$condition; pair <- des$pair

  ref0 <- which(ctcf_overlap(fc$coords, a$ctcf))
  message(sprintf("CTCF-overlapping consensus peaks: %d / %d", length(ref0), nrow(counts)))
  if (length(ref0) < as.integer(a$`min-refs`))
    warning("few CTCF anchors — shape estimate will be noisy")

  fit <- anchor_shape_offsets(counts, spike, ref0, condition,
                              span = as.numeric(a$span), trim_k = as.numeric(a$`trim-k`),
                              iter = as.integer(a$iter), min_refs = as.integer(a$`min-refs`))
  res <- run_deseq2(counts, fc$coords, condition, pair, fit$NF, fit$g0)

  # outputs
  utils::write.table(res, file.path(a$outdir, "differential_openness.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  ref_names <- fc$coords$Geneid[fit$ref_idx_final]
  writeLines(ref_names, file.path(a$outdir, "invariant_ctcf_anchors.txt"))
  lvl <- data.frame(sample = colnames(counts), spikein_reads = spike,
                    o_level = fit$o, check.names = FALSE)
  utils::write.table(lvl, file.path(a$outdir, "spikein_level.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  plot_diagnostics(a$outdir, fit, res, ref_names, fc$coords, condition, fit$ycorr)

  # run summary + the key diagnostic (spike-in vs CTCF-invariance agreement)
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  summ <- c(
    sprintf("consensus peaks           : %d", nrow(counts)),
    sprintf("CTCF anchors (initial)    : %d", length(ref0)),
    sprintf("CTCF anchors (after trim) : %d", length(fit$ref_idx_final)),
    sprintf("recovered global shift g0 : %.4f  log2  (accessibility %s under treatment)",
            fit$g0, ifelse(fit$g0 < 0, "globally DOWN", "globally UP")),
    sprintf("  |g0| interpretation     : ~0 => spike-in agrees with CTCF-invariance (robust);"),
    sprintf("                            large => spike-in is rescuing a global shift TMM/CTCF-only would delete"),
    sprintf("differential (padj<0.05)  : %d", n_sig))
  writeLines(summ, file.path(a$outdir, "run_summary.txt"))
  cat(paste(summ, collapse = "\n"), "\n")
  message("done -> ", a$outdir)
}

if (sys.nframe() == 0) main()
