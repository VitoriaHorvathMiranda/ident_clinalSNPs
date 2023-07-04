library(data.table)
library(qvalue)

snps <- fread("/dados/time_clines/output/SNPs/p-values_0910_no_ESC97_mincount5_minfreq0.001_cov15.tsv")

p_values <- snps$p_value




cutoffs <- c(0.1,0.05,0.01)

output <- vector("double", length(cutoffs))
for (i in seq_along(cutoffs)) {
  output[[i]] <- max(qobj_97_cov10$pvalues[qobj_97_cov10$qvalues <= cutoffs[[i]]])
}