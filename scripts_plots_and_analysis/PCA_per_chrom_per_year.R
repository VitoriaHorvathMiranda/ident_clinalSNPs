library(data.table)
library(tidyverse)
seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)

freqs <- fread("/dados/time_clines/data/seqs/align/NE_no_ESC97_mincount5_minfreq0.001_cov15.tsv")
all_samples <- freqs[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
#all_samples <- all_samples[, !c("CHROM", "POS")]

freqs_97 <- all_samples[population %like% "97"]
unique(freqs_97$population)

freqs_0910 <- all_samples[population %like% "09" | population %like% "10"]
unique(freqs_0910$population)


##new chuk
pre_pca_table <- dcast.data.table(freqs_0910, CHROM + POS + position2 ~ population,
                                  value.var = "freq")
no_NA_table <- na.omit(pre_pca_table)

chroms <- c("2L", "2R", "3L", "3R", "X")

#separate each chrom
output <- vector("list", length = length(chroms))
for (i in seq_along(chroms)) {
  output[[i]] <- no_NA_table[CHROM == chroms[[i]],][, !c("CHROM", "POS")]
}

#pre pca and pca
for (i in seq_along(output)) {
  output[[i]] <- column_to_rownames(output[[i]], var = "position2")
}

for (i in seq_along(output)) {
  output[[i]] <- t(output[[i]])
}

pca_per_chrom <- vector("list", length = length(chroms))
for (i in seq_along(output)) {
  pca_per_chrom[[i]] <- prcomp(output[[i]])
}

#gets pc1 and 2 for each chrom
pcs1e2_per_chrom <- vector("list", length = length(chroms))
for (i in seq_along(pca_per_chrom)) {
  pcs1e2_per_chrom[[i]] <- pca_per_chrom[[i]]$x[,1:2]
}

for (i in seq_along(pcs1e2_per_chrom)) {
  pcs1e2_per_chrom[[i]] <- as.data.frame(pcs1e2_per_chrom[[i]])
}

#change cols names
for (i in seq_along(pcs1e2_per_chrom)) {
  pcs1e2_per_chrom[[i]] <- setnames(pcs1e2_per_chrom[[i]],
                                    colnames(pcs1e2_per_chrom[[i]]),
                                    paste(colnames(pcs1e2_per_chrom[[i]]), chroms[[i]], sep = "_"))
}

all_chrom_pcs1e2 <- reduce(pcs1e2_per_chrom, cbind) %>%
  rownames_to_column("population")

#gets sdv per chrom
sdv_per_chrom <- vector("list", length = length(pca_per_chrom))
for (i in seq_along(pca_per_chrom)) {
  sdv_per_chrom[[i]] <- (pca_per_chrom[[i]]$sdev)^2
}

for (i in seq_along(sdv_per_chrom)) {
  sdv_per_chrom[[i]] <- round((sdv_per_chrom[[i]]/sum(sdv_per_chrom[[i]]))*100, 1)
}

sd_pcs_per_chrom <- reduce(sdv_per_chrom, cbind)
sd_pcs_per_chrom <- as_tibble(sd_pcs_per_chrom)
setnames(sd_pcs_per_chrom, colnames(sd_pcs_per_chrom), chroms)


tidy_sd_per_chrom <- sd_pcs_per_chrom %>%
  rownames_to_column("PCs") %>%
  pivot_longer(2:6, names_to = "chrom", values_to = "PCs_sd")

graph_table <- all_chrom_pcs1e2 %>%
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  )) %>%
  pivot_longer(cols = 2:11, names_to = "PCs",
               values_to = "var") %>%
  separate(col = PCs, into = c("PC", "chrom"), sep = "_") %>%
  pivot_wider(names_from = PC, values_from = var)


graph_table %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color = latitude), size = 2, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  facet_wrap(~ chrom, nrow = 2) +
  #xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  #ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_bw()

tidy_sd_per_chrom %>%
  mutate(PCs = as.integer(PCs)) %>%
  ggplot() +
  geom_bar(aes(x = PCs, y = PCs_sd), stat = "identity") +
  facet_wrap(~ chrom) +
  labs(x = "PCs", y = "var %")
