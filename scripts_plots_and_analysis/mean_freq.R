library(data.table)
library(tidyverse)

## 97 cov10
positions_97_cov10 <- fread("/dados/time_clines/output/SNPs/joined97_mincount10_minfreq0.01_cov15.tsv")

mean_freq_97_cov10 <- positions_97_cov10[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_freq_97_cov10) +
  geom_histogram(aes(x = V1), binwidth = 0.02)

ggplot(mean_freq_97_cov10) +
  geom_histogram(aes(x = V1), binwidth = 0.005) +
  coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1e+05))

## 0910 cov10
positions_0910_cov10 <- fread("/dados/time_clines/output/SNPs/joined0910_mincount10_minfreq0.01_cov15.tsv")

mean_freq_0910_cov10 <- positions_0910_cov10[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_freq_0910_cov10) +
  geom_histogram(aes(x = V1), binwidth = 0.02)

ggplot(mean_freq_97_cov10) +
  geom_histogram(aes(x = V1), binwidth = 0.005) +
  coord_cartesian(xlim = c(0, 0.1))

## all cov10
all <- fread("/dados/time_clines/data/seqs/align/NE_mincount5_minfreq0.005_cov10.tsv")

mean_all <- all[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_all) +
  geom_histogram(aes(x = V1), binwidth = 0.02) +
  labs(x = "mean freq")
  

ggplot(mean_all) +
  geom_histogram(aes(x = V1), binwidth = 0.003) +
  coord_cartesian(xlim = c(0, 0.1))

ggplot(mean_all) +
  geom_boxplot(aes(y = V1))

## cov15
all_cov15 <- fread("/dados/time_clines/data/seqs/align/NE_mincount10_minfreq0.01_cov15.tsv")

mean_all_cov15 <- all_cov15[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_all_cov15) +
  geom_histogram(aes(x = V1), binwidth = 0.01) +
  labs(x = "mean freq", title = "All Pops",
       caption = "min-cov = 15, min-freq = 0.01")

#1997
pops_97_cov15 <- all[population %like% "97",]
mean_pops_97_cov15 <- pops_97_cov15[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_pops_97_cov15) +
  geom_histogram(aes(x = V1), binwidth = 0.01) +
  labs(x = "mean freq", title = "1997",
       caption = "min-cov = 15, min-freq = 0.01")

#2009/2010
pops_0910cov15 <- all_cov15[population %like% "09" | population %like% "10",]
mean_pops_0910_cov15 <- pops_0910cov15[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_pops_0910_cov15) +
  geom_histogram(aes(x = V1), binwidth = 0.01) +
  labs(x = "mean freq", title = "2009/2010",
       caption = "min-cov = 10, min-freq = 0.01")






