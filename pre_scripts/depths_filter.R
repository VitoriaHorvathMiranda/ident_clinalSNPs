#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Filter samtools depth output, only positions called will be reported")
parser$add_argument('--SNPsposition', '-sp', help = 'positions of all snps called - output from freq_extraction.R')
parser$add_argument('--DepthFile', '-df', help= 'samtools depth file')
parser$add_argument('--output', '-o', help = 'tidy table of all depths of all called snps')
xargs<- parser$parse_args()

#read data  
positions <- fread(xargs$SNPsposition)

depths <- fread(xargs$DepthFile)

#fix header
oldnames <- paste0("V", 1:17)
newnames <- c("CHROM", "POS", "ESC97", "CMD97B", "WVT97", "MFL97", "RVA97", "SVT09",
              "CMD10", "SNC10", "JFL10", "CMA97", "HMA09", "HMA17",
              "MCT09", "MCT17", "MCT97")

setnames(depths, oldnames, newnames)

#filter positions
setkey(depths, CHROM, POS)
setkey(positions, CHROM, POS)
filtered_depths <- data.table::merge.data.table(positions, depths)

#
filtered_depths <- filtered_depths[, data.table::melt(.SD, measure.vars = 3:17,
                                   variable.name = "population",
                                   value.name = "depth")]

output <- xargs$output
write.table(filtered_depths, file = output, sep = "\t", row.names = FALSE)




