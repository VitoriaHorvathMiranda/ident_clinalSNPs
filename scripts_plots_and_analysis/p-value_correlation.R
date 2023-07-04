library(data.table)
library(tidyverse)

### cov10
SNPs_97_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_97_mincount5_minfreq0.005_cov10.tsv")
SNPs_0910_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_0910_mincount5_minfreq0.005_cov10.tsv")

new_names_97 <- c("position2", "coefficients_97", "p_value_97")
new_names_0910 <- c("position2", "coefficients_0910", "p_value_0910")

setnames(SNPs_97_cov10, names(SNPs_97_cov10), new_names_97)
setnames(SNPs_0910_cov10, names(SNPs_0910_cov10), new_names_0910)

merged_cov10 <- merge.data.table(SNPs_noESC97, SNPs_noESC97_0910)

ggplot(merged_cov10) +
  geom_bin2d(aes(x = p_value_97, y = p_value_0910), binwidth = c(0.01, 0.01))

ggplot(merged_cov10) +
  geom_bin2d(aes(x = coefficients_97, y = coefficients_0910), binwidth = c(2, 2))

cor.test(merged_cov10$p_value_97, merged_cov10$p_value_0910, method = "pearson")
cor.test(merged_cov10$coefficients_97, merged_cov10$coefficients_0910, method = "pearson")

p_value_lm <- lm(p_value_97 ~ p_value_0910, data = merged_cov10)
summary(p_value_lm)

coefficient_lm <- lm(coefficients_97 ~ coefficients_0910, data = merged_cov10)
summary(coefficient_lm)

ggplot(merged_cov10) +
  geom_bin2d(aes(x = p_value_97, y = p_value_0910),
             binwidth = c(0.001, 0.001)) +
  coord_cartesian(xlim = c(0.9925, 1), ylim = c(0.9925, 1))

merged_cov10 %>%
  mutate(coefficients_0910 = abs(coefficients_0910),
         coefficients_97 = abs(coefficients_97)) %>%
  ggplot() +
  geom_bin2d(aes(x = coefficients_97, y = coefficients_0910))
  

### cov15
SNPs_97_cov15 <- fread("/dados/time_clines/output/SNPs/p-values_97_mincount10_minfreq0.01_cov15.tsv")
SNPs_0910_cov15 <- fread("/dados/time_clines/output/SNPs/p-values_0910_mincount10_minfreq0.01_cov15.tsv")

new_names_97 <- c("position2", "coefficients_97", "p_value_97")
new_names_0910 <- c("position2", "coefficients_0910", "p_value_0910")

setnames(SNPs_97_cov15, names(SNPs_97_cov15), new_names_97)
setnames(SNPs_0910_cov15, names(SNPs_0910_cov15), new_names_0910)

merged_cov15 <- merge.data.table(SNPs_97_cov15, SNPs_0910_cov15)

ggplot(merged_cov15) +
  geom_bin2d(aes(x = p_value_97, y = p_value_0910), binwidth = c(0.01, 0.01))

ggplot(merged_cov15) +
  geom_bin2d(aes(x = coefficients_97, y = coefficients_0910), binwidth = c(2, 2))

cor.test(merged_cov15$p_value_97, merged_cov15$p_value_0910, method = "pearson")
cor.test(merged_cov15$coefficients_97, merged_cov15$coefficients_0910, method = "pearson")

p_value_lm <- lm(p_value_97 ~ p_value_0910,
                 data = merged_cov15)
summary(p_value_lm)

coefficient_lm <- lm(coefficients_97 ~ coefficients_0910,
                     data = merged_cov15)
summary(coefficient_lm)

ggplot(merged_cov15) +
  geom_bin2d(aes(x = p_value_97, y = p_value_0910),
             binwidth =  c(0.001, 0.001)) +
  coord_cartesian(xlim = c(0.9925, 1), ylim = c(0.9925, 1))


