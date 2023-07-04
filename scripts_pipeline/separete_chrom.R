#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "separates each chrom")
parser$add_argument('--freqs', '-freqs',
                    help = 'outputs from separate_time_pops_script.R')
parser$add_argument('--output', '-o',
                    nargs = '+',
                    help = 'table with exh chrom, .tsv')

xargs<- parser$parse_args()

freqs <- fread(xargs$freqs)

chr <- c("2L", "2R", "3L", "3R", "X")

files <- vector("list", length(chr))
for (i in seq_along(chr)) {
  files[[i]] <- freqs[CHROM == chr[[i]],]
}


for (i in seq_along(files)) {
  fwrite(files[[i]], xargs$output, sep = "\t")
}

