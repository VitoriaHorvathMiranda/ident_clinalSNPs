#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(qvalue)

#parse arguments
parser <- ArgumentParser(description= "compute significative SNPs for a given FDR")
parser$add_argument('--SNPs', '-snps',
                    help = 'outputs from glm_script.R')
parser$add_argument('--freqs', '-f',
                    help = 'outputs from separate_time_pops.R')
parser$add_argument('--FDR', '-q', default=0.1,
                    help = 'false discovery rate, ex:0.1')
parser$add_argument('--output', '-o',
                    help = 'table with significant snps for a given FDR, .tsv')
parser$add_argument('--histPlot', '-hp',
                    help = 'p-value histogram output, .png')
parser$add_argument('--qPlot', '-qp',
                    help = 'p-plot output, .png')
xargs<- parser$parse_args()

#reads data
file <- xargs$SNPs
SNPs <- fread(file = file)
freqs <- fread(xargs$freqs)

#finds SNPs that aren't present in that period or that are present in only
#one population
freqs[, position2 := paste(CHROM, POS, sep = ":")]
freqs[, freq_status := 1*(freq!= 0)]
no_0freq <- freqs[, .(freq_status_sum = sum(freq_status)), 
                  by = "position2"][freq_status_sum <= 1,]

#filter those SNPs out
SNPs <- SNPs[!no_0freq, on = "position2"]


#calculates the p-value cutoff for a 10% FDR
#based on Storey 2003 paper
FDR <- xargs$FDR
p_values <- SNPs$p_value
qobj <- qvalue(p = p_values)
p_value_cutoff <- max(qobj$pvalues[qobj$qvalues <= FDR])

#save important FDR plots
hist <- hist(qobj)
hist_plot <- xargs$histPlot
ggsave(hist_plot)

q_plots <- plot(qobj)
q_plots_output <- xargs$qPlot
ggsave(q_plots_output)

#gets snps below that p-value
significant_snps <- SNPs[p_value <= p_value_cutoff,]

#save table with significant snps
output_path <- xargs$output
fwrite(significant_snps, output_path, sep = "\t")





