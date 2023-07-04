#!/usr/bin/env Rscript
library(data.table)
library(argparse)

#parse arguments
parser <- ArgumentParser(description= "Filters clinal snps from vcf")
parser$add_argument('clinalSNPs', '-cSNPs', help = 'output from FDR_script.R')
parser$add_argument('--vcf', '-vcf', help= 'original vcf')
parser$add_argument('--output', '-o', help= 'a pre_vcf (vcf without header) with only clinal SNPs')
xargs<- parser$parse_args()

# get data
clinal_snps <- fread(xargs$clinalSNPs)
vcf <- fread(xargs$vcf)
setnames(vcf, "#CHROM", "CHROM")

# fix clinal_snps
clinal_snps[, c("CHROM", "POS") := tstrsplit(position2, ":")]
clinal_snps <- clinal_snps[, .(CHROM, POS)]
clinal_snps[, POS := as.integer(POS)]

#filter
pre_clinal_vcf <- vcf[clinal_snps, on = c("CHROM", "POS")]
setnames(pre_clinal_vcf, "CHROM", "#CHROM")

#save
fwrite(pre_clinal_vcf, file = xargs$output,
       sep = "\t")

