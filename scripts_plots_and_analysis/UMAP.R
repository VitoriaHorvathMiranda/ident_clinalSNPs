library(tidyverse)
library(data.table)
library(umap)


# wrangle -----------------------------------------------------------------
all <- fread("/dados/time_clines/data/seqs/align/NE_mincount10_minfreq0.01_cov15.tsv")
all_samples <- all[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
all_samples <- all_samples[, !c("CHROM", "POS")]

pre_pca_table <- dcast.data.table(all_samples, position2 ~ population,
                                  value.var = "freq")

pre_pca_table <- pre_pca_table[,!c("ESC97")]

no_NA_table <- na.omit(pre_pca_table)

pca_table <- no_NA_table %>%
  column_to_rownames(var = "position2")

pca_table <- t(pca_table)


# n_neighbors = 2 ---------------------------------------------------------

SNPs_umap <- umap(pca_table, n_neighbors = 2)

seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)


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
  labs(caption = "min-cov = 15; min-count = 10, min-freq = 1%,
       n_neighbors = 2")
  

# n_neighbors= 4 ----------------------------------------------------------

SNPs_umap4 <- umap(pca_table, n_neighbors = 4)

UMAP_tibble4 <- SNPs_umap4$layout %>%
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


UMAP_tibble4 %>%
  ggplot(aes(x = map1, y = map2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.1) +
  labs(caption = "min-cov = 15; min-count = 10, min-freq = 1%,
       n_neighbors = 4")


# cov 10 ------------------------------------------------------------------

all_cov10 <- fread("/dados/time_clines/data/seqs/align/NE_mincount5_minfreq0.005_cov10.tsv")
all_samples <- all_cov10[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
all_samples <- all_samples[, !c("CHROM", "POS")]

pre_pca_table <- dcast.data.table(all_samples, position2 ~ population,
                                  value.var = "freq")

pre_pca_table <- pre_pca_table[,!c("ESC97")]

no_NA_table <- na.omit(pre_pca_table)

pca_table <- no_NA_table %>%
  column_to_rownames(var = "position2")

pca_table <- t(pca_table)

# n_neighbors = 2 ---------------------------------------------------------

SNPs_umap <- umap(pca_table, n_neighbors = 2)

seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)


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
  labs(caption = "min-cov = 10; min-count = 5, min-freq = 0.5%,
       n_neighbors = 2")


# n_neighbors= 4 ----------------------------------------------------------

SNPs_umap4 <- umap(pca_table, n_neighbors = 4)

UMAP_tibble4 <- SNPs_umap4$layout %>%
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


UMAP_tibble4 %>%
  ggplot(aes(x = map1, y = map2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.1) +
  labs(caption = "min-cov = 10; min-count = 5, min-freq = 0.5%,
       n_neighbors = 4")




