library(data.table)
library(tidyverse)

freqs_cov10 <- fread("/dados/time_clines/data/seqs/align/NE_mincount5_minfreq0.005_cov10.tsv")

ggplot(freqs_cov10) +
  geom_histogram(aes(x= freq)) +
  facet_wrap(~ factor(population,
                      levels = c("ESC97", "CMD97B", "CMA97", "MCT97",
                                 "MFL97", "RVA97", "WVT97", "CMD10",
                                 "HMA09", "JFL10", "MCT09", "SNC10",
                                 "SVT09", "HMA17", "MCT17")))

freqs_cov15 <- fread("/dados/time_clines/data/seqs/align/NE_mincount10_minfreq0.01_cov15.tsv")

ggplot(freqs_cov15) +
  geom_histogram(aes(x= freq)) +
  facet_wrap(~ factor(population,
                      levels = c("ESC97", "CMD97B", "CMA97", "MCT97",
                                 "MFL97", "RVA97", "WVT97", "CMD10",
                                 "HMA09", "JFL10", "MCT09", "SNC10",
                                 "SVT09", "HMA17", "MCT17")))

