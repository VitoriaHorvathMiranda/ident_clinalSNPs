#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Computes effective number of chrom per position")
parser$add_argument('--freqs', '-f', help = 'output from n_chrom.R')
parser$add_argument('--output', '-o', help= 'tsv with effective number of chromossomes for each site')
xargs<- parser$parse_args()


#gets data
freqs <- fread(file = xargs$freqs)

#calcula o NE 
freqs[, NE := ((1/depth) + (1/n_chrom))^-1]

#save
output <- xargs$output
write.table(merged, file = output, sep = "\t", row.names = FALSE)

