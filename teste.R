library(tidyverse)
library(data.table)


#gets snps freqs
snapePool <- fread(file = "/dados/time_clines/data/seqs/calls/HMA09_L002_filtered.pool",
      sep = "\t", header = FALSE)

snape_col_names <-  c("chrom", "position",
              "ref_base", "ref_count",
              "alt_count", "ref_qual",
              "alt_qual", "bases",
              "prob", "p(1)", "freq")

setnames(snapePool, paste("V", 1:11, sep = ""), snape_col_names)

snapePool <- snapePool[, .(chrom, position, freq)]

#gets depth
samtoolsDepth <- fread(file = "/dados/vitoria/HMA09_samtoolsDepth.tsv",
                        sep = "\t", header = FALSE)
setnames(samtoolsDepth, paste("V", 1:3, sep = ""), c("chrom", "position", "depth"))

#gets depth for snapePool table
merged <- merge(snapePool, samtoolsDepth, by = c("chrom", "position"), all.x = TRUE)

#
chrom_sampled <- 184
latitude <- 1274
pop <- "HMA_09"

#calcula o NE 
merged[, NE := ((1/depth) + (1/chrom_sampled))^-1]
merged[, latitude := latitude]
merged[, pop := "HMA09"]

tab1_test <- fread(file = "tab1_test.tsv",
                   sep = "\t")

