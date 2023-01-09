#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Computes effective number of chrom per position")
parser$add_argument('--SnapeFile', '-sf', help = 'snape filtered file')
parser$add_argument('--DepthFile', '-df', help= 'samtools depth file')
parser$add_argument('--ChromNumber', '-chrom', type = "numeric",
                    help= 'number of chromossomes in that pool')
parser$add_argument('--PopName', '-pop', help= 'population name')
parser$add_argument('--Latitude', '-lat', help= 'population latitude', 
                    type = "numeric")
parser$add_argument('--output', '-o', help= 'tsv with effective number of chromossomes for each site')
xargs<- parser$parse_args()


#gets snps freqs
snapePool <- fread(file = xargs$SnapeFile,
                   sep = "\t", header = FALSE)

snape_col_names <-  c("chrom", "position",
                      "ref_base", "ref_count",
                      "alt_count", "ref_qual",
                      "alt_qual", "bases",
                      "prob", "p(1)", "freq")

setnames(snapePool, paste("V", 1:11, sep = ""), snape_col_names)

snapePool <- snapePool[, .(chrom, position, freq)]

#gets depth
samtoolsDepth <- fread(file = xargs$DepthFile,
                       sep = "\t", header = FALSE)
setnames(samtoolsDepth, paste("V", 1:3, sep = ""), c("chrom", "position", "depth"))

#gets depth for snapePool table
merged <- merge(snapePool, samtoolsDepth, by = c("chrom", "position"), all.x = TRUE)

#gets chrom number
chrom_sampled <- xargs$ChromNumber

#calcula o NE 
merged[, NE := ((1/depth) + (1/chrom_sampled))^-1]
merged[, latitude := xargs$Latitude]
merged[, population := xargs$PopName]

#save
output <- xargs$output
write.table(merged, file = output, sep = "\t", row.names = FALSE)

