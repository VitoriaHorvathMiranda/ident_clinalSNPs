library(data.table)
library(tidyverse)
library(qvalue)


SNPs_97_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_97_mincount5_minfreq0.005_cov10.tsv")
SNPs_97_cov15 <- fread("/dados/time_clines/output/SNPs/p-values_97_mincount10_minfreq0.01_cov15.tsv")

SNPs_0910_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_0910_mincount5_minfreq0.005_cov10.tsv")
SNPs_0910_cov15 <- fread("/dados/time_clines/output/SNPs/p-values_0910_mincount10_minfreq0.01_cov15.tsv")

SNPs_noESC97 <- fread("/dados/time_clines/output/SNPs/p-values_97_no_ESC97_mincount5_minfreq0.005_cov15.tsv")
SNPs_noESC97_0910 <- fread("/dados/time_clines/output/SNPs/p-values_0910_no_ESC97_mincount5_minfreq0.005_cov15.tsv")


# histogramas -------------------------------------------------------------
SNPs_0910_cov10 %>%
  ggplot() +
  geom_histogram(aes(p_value)) +
  labs(x = "p-value", title = "2009/2010",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

SNPs_0910_cov15 %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "2009/2010",
       caption = "min-freq = 0.01, min-cov = 15, min-count = 10") +
  theme_minimal()


SNPs_97_cov10 %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "1997",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

SNPs_97_cov15 %>%
  ggplot() +
  geom_histogram(aes(p_value)) +
  labs(x = "p-value", title = "1997",
       caption = "min-freq = 0.01, min-cov = 15, min-count = 10") +
  theme_minimal()


# q-value -----------------------------------------------------------------
##q-value
Pvalue_97_cov10 <- SNPs_97_cov10$p_value
Pvalue_97_cov15 <- SNPs_97_cov15$p_value

qobj_97_cov10 <- qvalue(p = Pvalue_97_cov10)
qobj_97_cov15 <- qvalue(p = Pvalue_97_cov15)


p_value_cutoff_0.1 <- max(qobj_97_cov10$pvalues[qobj_97_cov10$qvalues <= 0.1])

hist <- hist(qobj_97_cov15)
hist_plot <- xargs$histPlot
ggsave(hist_plot)


q_plots <- plot(qobj)
q_plots_output <- xargs$qPlot
ggsave(q_plots_output)

significant_snps <- SNPs_97_cov15[p_value <= p_value_cutoff_0.01,]

# test --------------------------------------------------------------------
big_p_values_97_cov10 <- SNPs_97_cov10[p_value >= 0.99, ]
big_p_values_97_cov15 <- SNPs_97_cov15[p_value >= 0.99, ]
big_p_values_0910_cov10 <- SNPs_0910_cov10[p_value >= 0.99, ]
big_p_values_0910_cov15 <- SNPs_0910_cov15[p_value >= 0.99, ]

#09/10 vs 09/10
big_p_values_0910_cov10 %>%
  anti_join(big_p_values_0910_cov15, by = "position2")#tem 74751 snps com pvalor alto em 0910 10 que não existem em 0910 15

big_p_values_0910_cov15 %>%
  anti_join(big_p_values_0910_cov10, by = "position2")#tem 84 snps que tem p-valor alto em 0910 15 e não 0910 10

#97 vs 09/10 15
big_p_values_97_cov15 %>%
  anti_join(big_p_values_0910_cov10, by = "position2")

big_p_values_0910_cov15 %>%
  anti_join(big_p_values_97_cov10, by = "position2")

#97 vs 09/10 10
big_p_values_97_cov10 %>%
  anti_join(big_p_values_0910_cov10, by = "position2") # 279649

big_p_values_0910_cov10 %>%
  anti_join(big_p_values_97_cov10, by = "position2") #66874



