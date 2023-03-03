#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "gets latitude and Nº of chrom")
parser$add_argument('--metadata', '-meta', help = 'metadata table')
parser$add_argument('--DepthFreq', '-df', help= 'output from freq_extraction.R')
parser$add_argument('--output', '-o', help= 'tsv with depth, latitude and n of chrom for each snp')
xargs<- parser$parse_args()

#get data
metadata <- fread(file = xargs$metadata)
depths <- fread(file = xargs$DepthFreq)

depths <- depths[depth != ".", ]
depths[, depth := as.double(depth)]

#get relevant columns 
I_info <- metadata[, .(population, n_females, latitude)]

#compute n of chrom
I_info <- I_info[, n_chrom := n_females*2][,!c("n_females")]

#join 
full <- merge.data.table(I_info, depths, by = "population")

#calcula o NE 
full[, NE := ((1/depth) + (1/n_chrom))^-1]

#save
output <- xargs$output
write.table(full, file = output, sep = "\t", row.names = FALSE)



