#!/usr/bin/env Rscript
# Unit test for stable_genes_from_de(). Run in the r-diffopen conda env:
#   Rscript tests/r/test_stable_genes.R
source("workflow/scripts/diffopen.R")

de <- data.frame(
  gene           = c("STABLE1","UP1","LOWEXPR","SIG1","NAPADJ","BORDER","NALFC","NABM","NEGBORDER"),
  baseMean       = c(100,      100,  5,        100,   100,     10,       100,    NA,   100),
  log2FoldChange = c(0.1,      2.0,  0.1,      0.2,   0.1,     0.5,      NA,     0.1,  -0.5),
  padj           = c(0.9,      0.001,0.9,      0.01,  NA,      0.5,      0.9,    0.9,  0.9),
  stringsAsFactors = FALSE)

got <- sort(stable_genes_from_de(de))
want <- sort(c("STABLE1","NAPADJ","BORDER","NEGBORDER"))
stopifnot(identical(got, want))

# column-name overrides
de2 <- de; names(de2)[names(de2) == "gene"] <- "symbol"
stopifnot(identical(sort(stable_genes_from_de(de2, gene_col = "symbol")), want))

# missing column errors
ok <- tryCatch({ stable_genes_from_de(de, gene_col = "nope"); FALSE },
               error = function(e) TRUE)
stopifnot(ok)

cat("test_stable_genes.R OK\n")
