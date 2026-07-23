#!/usr/bin/env Rscript
# Run in the r-diffopen conda env: Rscript tests/r/test_size_factors_rnastable.R
suppressPackageStartupMessages({ library(GenomicRanges); library(IRanges) })
source("workflow/scripts/diffopen.R")

set.seed(1)
G <- 300; n <- 4
counts <- matrix(rpois(G * n, 200), nrow = G,
                 dimnames = list(paste0("p", 1:G), paste0("s", 1:n)))
coords <- data.frame(Geneid = paste0("p", 1:G), Chr = "chr1",
                     Start = seq(1, by = 100, length.out = G),
                     End   = seq(60, by = 100, length.out = G),
                     stringsAsFactors = FALSE)
# All peaks promoter-class; every peak sits under one stable gene's TSS window.
promoter_is <- rep(TRUE, G)
tx <- GRanges("chr1", IRanges(start = coords$Start + 5, width = 1), strand = "+")
tx$gene_name <- paste0("GENE", 1:G)
de <- data.frame(gene = paste0("GENE", 1:G), baseMean = 100,
                 log2FoldChange = 0, padj = 0.9, stringsAsFactors = FALSE)
condition <- factor(c("Control","Control","Treat","Treat"),
                    levels = c("Control","Treat"))

res <- size_factors_rnastable(counts, coords, promoter_is, tx, de,
        gene_col = "gene", lfc_col = "log2FoldChange", padj_col = "padj",
        basemean_col = "baseMean", basemean_min = 10, padj_min = 0.5, lfc_max = 0.5,
        window = 50, min_anchors = 100, promoter_class_required = TRUE,
        condition = condition, trim_k = 2.5, iter = 2)

stopifnot(length(res$sf) == n)
stopifnot(abs(exp(mean(log(res$sf))) - 1) < 1e-8)   # geometric mean 1
stopifnot(res$n_stable == G)
stopifnot(abs(res$match_rate - 1) < 1e-9)
stopifnot(res$n_anchor == G)

# too few anchors -> stop()
too_few <- tryCatch({
  size_factors_rnastable(counts, coords, promoter_is, tx, de,
    gene_col = "gene", lfc_col = "log2FoldChange", padj_col = "padj",
    basemean_col = "baseMean", basemean_min = 10, padj_min = 0.5, lfc_max = 0.5,
    window = 50, min_anchors = 100000, promoter_class_required = TRUE,
    condition = condition, trim_k = 2.5, iter = 2); FALSE
}, error = function(e) grepl("anchors", conditionMessage(e)))
stopifnot(too_few)

cat("test_size_factors_rnastable.R OK\n")
