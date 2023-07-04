library(tidyverse)
library(data.table)
library(splitstackshape)

vcf <- fread("/dados/time_clines/data/seqs/calls/PoolSNP_no_ESC97_mincount5_minfreq0.005_cov15.vcf", skip = "##")

setnames(vcf, "#CHROM", "CHROM")
chrom <- c("2L", "2R", "3L", "3R", "X")

vcf <- vcf[CHROM %in% c("2L", "2R", "3L", "3R", "X"), ]

vcf[, unique(CHROM)]

freqs <- vcf[, cSplit(.SD ,splitCols = c("CMD97B","WVT97","MFL97","RVA97","SVT09","CMD10","SNC10","JFL10",
                        "CMA97","HMA09","HMA17","MCT09","MCT17","MCT97"),
                      sep = ":")]

drop_cols <- grep(c("_[123]"), colnames(freqs))
freqs[, (drop_cols) := NULL]

called_snps <- freqs[, .(CHROM, POS)]

freqs <- freqs[, !c("ID", "QUAL", "FILTER", "INFO", "FORMAT")]

newnames <- colnames(freqs) %>% str_replace_all(c("_5" = "_freq", "_4" = "_depth")) 
oldnames <- colnames(freqs)
setnames(freqs, oldnames, newnames)

tiny_freqs <- freqs[1:10000]
tiny_depths <- tiny_freqs[, !(names(tiny_freqs) %like% "freq"), with = FALSE]
tiny_freqs <- tiny_freqs[, !(names(tiny_freqs) %like% "depth"), with = FALSE]

newnames2 <- colnames(tiny_freqs) %>% str_replace("_freq", "") 
oldnames_freq <- colnames(tiny_freqs)
oldnames_depth <- colnames(tiny_depths)
  
setnames(tiny_freqs, oldnames_freq, newnames2)
setnames(tiny_depths, oldnames_depth, newnames2)


tiny_freqs <- tiny_freqs[, data.table::melt(.SD, measure.vars = 5:19,
                        variable.name = "population",
                        value.name = "freq")]

tiny_depths <- tiny_depths[, data.table::melt(.SD, measure.vars = 5:19,
                                            variable.name = "population",
                                            value.name = "depth")]

merge.data.table(tiny_depths, tiny_freqs)

write.table(freq_depth, file = "/dados/time_clines/data/seqs/calls/freqs_and_depths_no_ESC97_mincount5_minfreq0.005_cov15.tsv", sep = "\t", row.names = FALSE)









