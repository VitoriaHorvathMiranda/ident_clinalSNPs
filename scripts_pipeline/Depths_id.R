#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "gets pop name and latitude")
parser$add_argument('--DepthFile', '-df', help= 'samtools depth file')
parser$add_argument('--PopName', '-pop', help= 'population name')
parser$add_argument('--Latitude', '-lat', help= 'population latitude', 
                    type = "numeric")
parser$add_argument('--output', '-o', help= 'tsv with depth name and latitude')
xargs<- parser$parse_args()

#gets depth
samtoolsDepth <- fread(file = xargs$DepthFile,
                       sep = "\t", header = FALSE)
setnames(samtoolsDepth, paste("V", 1:3, sep = ""), c("chrom", "position", "depth"))

#adds name and latitude
samtoolsDepth[, latitude := xargs$Latitude]
samtoolsDepth[, population := xargs$PopName]

#save
output <- xargs$output
write.table(samtoolsDepth, file = output, sep = "\t", row.names = FALSE)




