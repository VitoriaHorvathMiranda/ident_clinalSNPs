#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(splitstackshape)

#parse arguments
parser <- ArgumentParser(description= "Gets SNPs freq from PoolSNP vcf")
parser$add_argument('--PoolSNPvcf', '-vcf', help = 'output from PoolSNP')
parser$add_argument('--calledSNPs', '-snps', help= 'output of all snps called')
parser$add_argument('--output', '-o', help= 'a tidy table with freqs and depths for each SNP')
xargs<- parser$parse_args()

#reads vcf
vcf <- fread(xargs$PoolSNPvcf, skip = "##")
setnames(vcf, "#CHROM", "CHROM")

#filter chrom
vcf <- vcf[CHROM %in% c("2L", "2R", "3L", "3R", "X"), ]

#split pop columns 
freqs <- vcf[, cSplit(.SD ,splitCols = c("ESC97","CMD97B","WVT97","MFL97","RVA97","SVT09","CMD10","SNC10","JFL10",
                                         "CMA97","HMA09","HMA17","MCT09","MCT17","MCT97"),
                      sep = ":")]
#drop all non freq columns 
drop_cols <- grep(c("_[123]"), colnames(freqs))
freqs[, (drop_cols) := NULL]
freqs <- freqs[, !c("ID", "QUAL", "FILTER", "INFO", "FORMAT")]

#keep only biallelic snps
freqs <- freqs[!(ALT %like% ","),]

#gets all called snps
#it will be used to filter the depths (from samtools depth)
called_snps <- freqs[, .(CHROM, POS)]

#fix names
newnames <- colnames(freqs) %>% str_replace_all(c("_5" = "_freq", "_4" = "_depth")) 
oldnames <- colnames(freqs)
setnames(freqs, oldnames, newnames)

#subset freqs and depths
depths <- freqs[, !(names(freqs) %like% "freq"), with = FALSE]
freqs <- freqs[, !(names(freqs) %like% "depth"), with = FALSE]

#fix names again
newnames2 <- colnames(freqs) %>% str_replace("_freq", "") 
oldnames_freq <- colnames(freqs)
oldnames_depth <- colnames(depths)

setnames(freqs, oldnames_freq, newnames2)
setnames(depths, oldnames_depth, newnames2)

#make freqs tidy
freqs <- freqs[, data.table::melt(.SD, measure.vars = 5:19,
                                  variable.name = "population",
                                  value.name = "freq")]

#make depths tidy
depths <- depths[, data.table::melt(.SD, measure.vars = 5:19,
                                  variable.name = "population",
                                  value.name = "depth")]

#merge
freq_depth <- merge.data.table(depths, freqs)

#saves
output <- xargs$output
snps <- xargs$calledSNPs
write.table(called_snps, file = snps, sep = "\t", row.names = FALSE)
write.table(freq_depth, file = output, sep = "\t", row.names = FALSE)



