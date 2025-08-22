#!/usr/bin/env Rscript
suppressPackageStartupMessages({
library(optparse); library(data.table)
})


option_list <- list(
make_option(c('--sumstats'), type='character'),
make_option(c('--ldref'), type='character'),
make_option(c('--out'), type='character')
)
opt <- parse_args(OptionParser(option_list=option_list))


ss <- fread(opt$sumstats)
# TODO: implement full LDpred2 using bigsnpr/bigstatsr with external LD ref.
# For toy flow, write passthrough weights to keep pipeline running.
out <- ss[, .(SNP, A1, BETA)]
fwrite(out, opt$out, sep='\t')
cat('LDpred2 stub: copied GWAS betas (replace with real LDpred2 implementation).\n')

