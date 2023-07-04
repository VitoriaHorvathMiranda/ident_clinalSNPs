library(data.table)
library(tidyverse)
library(umap)

p_value_97 <- fread("/dados/time_clines/output/SNPs/p-values_97_joined_no_ESC97_mincount5_minfreq0.001_cov15.tsv")
p_value_0910 <- fread("/dados/time_clines/output/SNPs/p-values_0910_joined_no_ESC97_mincount5_minfreq0.001_cov15.tsv")
setnames(p_value_97, names(p_value_97), c("position2", "coefficient", "p_value"))
setnames(p_value_0910, names(p_value_0910), c("position2", "coefficient", "p_value"))

freqs <- fread("/dados/time_clines/data/seqs/align/NE_no_ESC97_mincount5_minfreq0.001_cov15.tsv")
seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)


# p-value hist ------------------------------------------------------------

p_value_0910 %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "2009/2010",
       caption = "min-freq = 0.005, min-cov = 20, min-count = 5") +
  theme_minimal()

p_value_97 %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "1997",
       caption = "min-freq = 0.005, min-cov = 20, min-count = 5") +
  theme_minimal()



# freq dist ---------------------------------------------------------------

mean_all <- freqs[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_all) +
  geom_histogram(aes(x = V1), binwidth = 0.02) +
  labs(x = "mean freq", title = "All Pops",
       caption = "min-cov = 10, min-freq = 0.005")

#1997
pops_97 <- freqs[population %like% "97",]
mean_pops_97 <- pops_97[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_pops_97) +
  geom_histogram(aes(x = V1), binwidth = 0.02) +
  labs(x = "mean freq", title = "1997",
       caption = "min-cov = 10, min-freq = 0.005")

#2009/2010
pops_0910 <- freqs[population %like% "09" | population %like% "10",]
mean_pops_0910 <- pops_0910[, mean(freq), by = c("CHROM", "POS")]

ggplot(mean_pops_0910) +
  geom_histogram(aes(x = V1), binwidth = 0.02) +
  labs(x = "mean freq", title = "2009/2010",
       caption = "min-cov = 10, min-freq = 0.005")

ggplot(freqs) +
  geom_histogram(aes(x= freq)) +
  facet_wrap(~ factor(population,
                      levels = c("ESC97", "CMD97B", "CMA97", "MCT97",
                                 "MFL97", "RVA97", "WVT97", "CMD10",
                                 "HMA09", "JFL10", "MCT09", "SNC10",
                                 "SVT09", "HMA17", "MCT17"))) +
  labs(caption = "min-cov = 10; min-count = 5, min-freq = 0.5%")


# heatmap -----------------------------------------------------------------

new_names_97 <- c("position2", "coefficients_97", "p_value_97")
new_names_0910 <- c("position2", "coefficients_0910", "p_value_0910")

setnames(p_value_97, names(p_value_97), new_names_97)
setnames(p_value_0910, names(p_value_0910), new_names_0910)

merged_cov20 <- merge.data.table(p_value_97, p_value_0910)

ggplot(merged_cov20) +
  geom_bin2d(aes(x = p_value_97, y = p_value_0910), binwidth = c(0.01, 0.01))


# PCA ---------------------------------------------------------------------

all_samples <- freqs[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
all_samples <- all_samples[, !c("CHROM", "POS")]
pre_pca_table <- dcast.data.table(all_samples, position2 ~ population,
                                  value.var = "freq")

no_NA_table <- na.omit(pre_pca_table)

pca_table <- no_NA_table %>%
  column_to_rownames(var = "position2")

pca_table <- t(pca_table)

pca <- prcomp(pca_table)

#var per pc

pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

tibble(var = pca_var_per) %>%
  rownames_to_column() %>%
  mutate(PCs = as.integer(rowname)) %>%
  select(-rowname) %>%
  ggplot(aes(x = PCs, y = var)) +
  geom_bar(aes(), stat = "identity" ) +
  labs(y = "var %")

#pca total plot
graph_table <- tibble(population = rownames(pca$x), 
                      pc1 = pca$x[,1], #gets the first four PCs
                      pc2 = pca$x[,2],
                      pc3 = pca$x[,3],
                      pc4 = pca$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  ))

graph_table %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_bw()

#97
pre_pca_table97 <- pre_pca_table[, !(names(pre_pca_table) %like% "09" | names(pre_pca_table) %like% "10" | names(pre_pca_table) %like% "17"),
                                 with = FALSE]

no_NA_table97 <- na.omit(pre_pca_table97)

pca_table97 <- no_NA_table97 %>%
  column_to_rownames(var = "position2")

pca_table97 <- t(pca_table97)

pca97 <- prcomp(pca_table97)

pca_var97 <- pca97$sdev^2
pca_var_per97 <- round(pca_var97/sum(pca_var97)*100, 1)

graph_table97 <- tibble(population = rownames(pca97$x), 
                        pc1 = pca97$x[,1], #gets the first four PCs
                        pc2 = pca97$x[,2],
                        pc3 = pca97$x[,3],
                        pc4 = pca97$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = as.character(collection_year))

graph_table97 %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per97[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per97[2], "%", sep = "")) +
  theme_bw()

#2009/2010
pre_pca_table0910 <- pre_pca_table[, !(names(pre_pca_table) %like% "97" | names(pre_pca_table) %like% "17"),
                                   with = FALSE]

no_NA_table0910 <- na.omit(pre_pca_table0910)

pca_table0910 <- no_NA_table0910 %>%
  column_to_rownames(var = "position2")

pca_table0910 <- t(pca_table0910)

pca0910 <- prcomp(pca_table0910)

pca_var0910 <- pca0910$sdev^2
pca_var_per0910 <- round(pca_var0910/sum(pca_var0910)*100, 1)

graph_table0910 <- tibble(population = rownames(pca0910$x), 
                          pc1 = pca0910$x[,1], #gets the first four PCs
                          pc2 = pca0910$x[,2],
                          pc3 = pca0910$x[,3],
                          pc4 = pca0910$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = as.character(collection_year))

graph_table0910 %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per0910[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per0910[2], "%", sep = "")) +
  theme_bw()


# UMAP --------------------------------------------------------------------
##n_neighbors = 2

SNPs_umap <- umap(pca_table, n_neighbors = 2)

UMAP_tibble <- SNPs_umap$layout %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename("population" = "rowname",
         "map1" = "V1",
         "map2" = "V2") %>%
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  ))


UMAP_tibble %>%
  ggplot(aes(x = map1, y = map2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.1) +
  labs(caption = "min-cov = 20; min-count = 5, min-freq = 0.5%,
       n_neighbors = 2")

