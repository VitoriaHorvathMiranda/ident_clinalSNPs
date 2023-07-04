library(data.table)
library(tidyverse)
library(qvalue)


SNPs_97_cov15 <- fread("/dados/time_clines/output/SNPs/p-values_97_no_ESC97_mincount5_minfreq0.001_cov15.tsv")

SNPs_97_cov15[, c("CHROM", "POS") := tstrsplit(position2, ":", fixed = TRUE)]

SNPs_97_cov15 <- SNPs_97_cov15[, .(CHROM, POS, p_value)]
SNPs_97_cov15[, POS := as.double(POS)]
setkey(SNPs_97_cov15, CHROM)
typeof(SNPs_97_cov15$POS)
SNPs_97_cov15[, ID := .I]

#### manhattan plot with qqman-----------------
setkey(SNPs_97_cov15, CHROM, POS)


SNPs_97_cov15[, CHROM_N := case_when(CHROM == "2L" ~ 1,
                                     CHROM == "2R" ~ 2,
                                     CHROM == "3L" ~ 3,
                                     CHROM == "3R" ~ 4,
                                     CHROM == "X" ~ 5)]

manhattan(SNPs_97_cov15, chr="CHROM_N", bp="POS", snp="ID", p="p_value" )

######get FDR cutof ---------------
Pvalue_97_cov15 <- SNPs_97_cov15$p_value
qobj_97_cov15 <- qvalue(p = Pvalue_97_cov15)

#cutoff 1%
p_value_cutoff_0.1 <- max(qobj_97_cov15$pvalues[qobj_97_cov15$qvalues <= 0.1])

#cutoff 1%
p_value_cutoff_0.01 <- max(qobj_97_cov15$pvalues[qobj_97_cov15$qvalues <= 0.01])

#cutoff 0.05%
p_value_cutoff_0.005 <- max(qobj_97_cov15$pvalues[qobj_97_cov15$qvalues <= 0.005])


#### manhattan plot manual ---------------
tiny_SNPs <- SNPs_97_cov15[-log10(p_value) > 1, ] 

axisdf <-  SNPs_97_cov15 %>% group_by(CHROM) %>% summarize(center=( max(ID) + min(ID) ) / 2 )

SNPs_97_cov15 %>%
  group_by(CHROM) %>%
  summarise(max = max(ID),
            min = min(ID),
            center = (max + min)/2)

tiny_SNPs %>%
  ggplot(aes(x = ID, y =-log10(p_value))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=0.9) +
  geom_line(aes(y = -log10(p_value_cutoff_0.1)), color = "green") +
  geom_line(aes(y = -log10(p_value_cutoff_0.01)), color = "blue") +
  geom_line(aes(y = -log10(p_value_cutoff_0.005)), color = "red") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.1) + 0.25,
                label = "FDR 10%"),
            size = 2.5, color = "green") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.01) + 0.25,
                label = "FDR 1%"),
            size = 2.5, color = "blue") +
  geom_text(aes(x = max(ID)+10,
                y = -log10(p_value_cutoff_0.005) + 0.25,
                label = "FDR 0.5%"),
            size = 2.5, color = "red") +
  scale_color_manual(values = rep(c("lavender", "lavenderblush1"), 5)) +
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) + 
  labs(title = "1997", y = "-log10(p-value)")

-log10(p_value_cutoff_0.1)


?labs
