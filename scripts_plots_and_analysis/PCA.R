library(tidyverse)
library(data.table)
library(viridis)

#gets other info that we want
seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)

# wrangle -----------------------------------------------------------------
all <- fread("/dados/time_clines/data/seqs/align/NE_mincount5_minfreq0.005_cov10.tsv")
all_samples <- all[, .(population, CHROM, POS, freq)]
all_samples[, position2 := paste(CHROM, POS, sep = ":")]
#all_samples <- all_samples[, !c("CHROM", "POS")]

pre_pca_table <- dcast.data.table(all_samples, CHROM + POS + position2 ~ population,
                 value.var = "freq")

pre_pca_table <- pre_pca_table[,!c("ESC97")]

no_NA_table <- na.omit(pre_pca_table)

pca_table <- no_NA_table %>%
  column_to_rownames(var = "position2")

pca_table_t <- t(pca_table)

# pre_pca and pca -----------------------------------------------------------------

pca <- prcomp(pca_table)
pca$x

#gets the percentage of variation that each PC accounts
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)

#plots that variation
tibble(var = pca_var_per) %>%
  rownames_to_column() %>%
  mutate(PCs = as.integer(rowname)) %>%
  select(-rowname) %>%
  ggplot(aes(x = PCs, y = var)) +
  geom_bar(aes(), stat = "identity" ) +
  labs(x = "PCs", y = "var %")


#joins all info in a single table
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


#plot colored by latitude
graph_table %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  theme_bw()


#plot colored by collection_year
graph_table %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = collection_year), size = 1.3, alpha = 0.5) +
  #geom_text(aes(label = population), size = 2, nudge_y = 12) +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = ""))


# pc2_pc3 -----------------------------------------------------------------

#plot colored by latitude
graph_table %>%
  ggplot(aes(x = pc2, y = pc3)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_var_per[3], "%", sep = ""))


#plot colored by collection_year
graph_table %>%
  ggplot(aes(x = pc2, y = pc3)) +
  geom_point(aes(color = collection_year), size = 1.3, alpha = 0.5) +
  geom_text(aes(label = population), size = 2, nudge_y = 12) +
  xlab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  ylab(paste("PC3 - ", pca_var_per[3], "%", sep = ""))


# pc3_pc4 -----------------------------------------------------------------

#plot colored by latitude
graph_table %>%
  ggplot(aes(x = pc3, y = pc4)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC3 - ", pca_var_per[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_var_per[4], "%", sep = ""))


#plot colored by collection_year
graph_table %>%
  ggplot(aes(x = pc3, y = pc4)) +
  geom_point(aes(color = collection_year), size = 1.3, alpha = 0.5) +
  geom_text(aes(label = population), size = 2, nudge_y = 12) +
  xlab(paste("PC3 - ", pca_var_per[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_var_per[4], "%", sep = ""))


# pca without CMD97B -------------------------------------------------------

no_CMD97B <- pre_pca_table[,!c("CMD97B")]

no_NA_CMD97B <- na.omit(no_CMD97B)

pca_table_CMD97B <- no_NA_CMD97B %>%
  column_to_rownames(var = "position2")

pca_table_CMD97B <- t(pca_table_CMD97B)

pca_no_CMD97B <- prcomp(pca_table_CMD97B)

#gets the percentage of variation that each PC accounts
pca_no_CMD97B_var <- pca_no_CMD97B$sdev^2
pca_no_CMD97B_var_per <- round(pca_no_CMD97B_var/sum(pca_no_CMD97B_var)*100, 1)

#plots that variation
tibble(var = pca_no_CMD97B_var_per) %>%
  rownames_to_column() %>%
  mutate(rowname = as.integer(rowname)) %>%
  ggplot(aes(x = rowname, y = var)) +
  geom_bar(aes(), stat = "identity" ) +
  labs(y = "var %", x = "PCs")
  


#joins all info in a single table
graph_table_noCMD97B <- tibble(population = rownames(pca_no_CMD97B$x), 
                      pc1 = pca_no_CMD97B$x[,1], #gets the first two PCs
                      pc2 = pca_no_CMD97B$x[,2],
                      pc3 = pca_no_CMD97B$x[,3],
                      pc4 = pca_no_CMD97B$x[,4]) %>%
  mutate(population = str_replace(population, "_freq", "")) %>%
  inner_join(pop_info, by = "population") %>% #joins tables
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  ))


#plot colored by latitude
graph_table_noCMD97B %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_no_CMD97B_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_no_CMD97B_var_per[2], "%", sep = ""))





# pca 97 ------------------------------------------------------------------
pre_pca_table97 <- pre_pca_table[, !(names(pre_pca_table) %like% "09" | names(pre_pca_table) %like% "10" | names(pre_pca_table) %like% "17"),
                                 with = FALSE]

no_NA_table97 <- na.omit(pre_pca_table97)

pca_table97 <- no_NA_table97 %>%
  column_to_rownames(var = "position2")

pca_table97 <- t(pca_table97)

pca97 <- prcomp(pca_table97)

#gets the percentage of variation that each PC accounts
pca_var97 <- pca97$sdev^2
pca_var_per97 <- round(pca_var97/sum(pca_var97)*100, 1)

#plots that variation
tibble(var = pca_var_per97) %>%
  rownames_to_column() %>%
  mutate(PCs = as.integer(rowname)) %>%
  select(-rowname) %>%
  ggplot(aes(x = PCs, y = var)) +
  geom_bar(aes(), stat = "identity") +
  labs(y = "var %", x = "PCs", title = "1997")

#joins all info in a single table
graph_table97 <- tibble(population = rownames(pca97$x), 
                      pc1 = pca97$x[,1], #gets the first four PCs
                      pc2 = pca97$x[,2],
                      pc3 = pca97$x[,3],
                      pc4 = pca97$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = as.character(collection_year))


#plot colored by latitude
graph_table97 %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per97[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per97[2], "%", sep = "")) +
  theme_bw()

graph_table97 %>%
  ggplot(aes(x = pc3, y = pc4)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC3 - ", pca_var_per97[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_var_per97[4], "%", sep = "")) +
  theme_bw()


# pca 0910 ----------------------------------------------------------------
pre_pca_table0910 <- pre_pca_table[, !(names(pre_pca_table) %like% "97" | names(pre_pca_table) %like% "17"),
                                 with = FALSE]

no_NA_table0910 <- na.omit(pre_pca_table0910)

pca_table0910 <- no_NA_table0910 %>%
  column_to_rownames(var = "position2")

pca_table0910 <- t(pca_table0910)

pca0910 <- prcomp(pca_table0910)

#gets the percentage of variation that each PC accounts
pca_var0910 <- pca0910$sdev^2
pca_var_per0910 <- round(pca_var0910/sum(pca_var0910)*100, 1)

#plots that variation
tibble(var = pca_var_per0910) %>%
  rownames_to_column() %>%
  mutate(PCs = as.integer(rowname)) %>%
  select(-rowname) %>%
  ggplot(aes(x = PCs, y = var)) +
  geom_bar(aes(), stat = "identity") +
  labs(y = "var %", x = "PCs", title = "2009/2010")

#joins all info in a single table
graph_table0910 <- tibble(population = rownames(pca0910$x), 
                        pc1 = pca0910$x[,1], #gets the first four PCs
                        pc2 = pca0910$x[,2],
                        pc3 = pca0910$x[,3],
                        pc4 = pca0910$x[,4]) %>% 
  inner_join(pop_info, by = "population") %>%  #joins tables
  mutate(collection_year = as.character(collection_year))


#plot colored by latitude
graph_table0910 %>%
  ggplot(aes(x = pc1, y = pc2)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC1 - ", pca_var_per0910[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per0910[2], "%", sep = "")) +
  theme_bw()

graph_table0910 %>%
  ggplot(aes(x = pc3, y = pc4)) +
  geom_point(aes(color = latitude), size = 4, alpha = 0.5) +
  scale_color_viridis_c() +
  geom_text(aes(label = population), size = 1.5, nudge_y = 10) +
  xlab(paste("PC3 - ", pca_var_per0910[3], "%", sep = "")) +
  ylab(paste("PC4 - ", pca_var_per0910[4], "%", sep = "")) +
  theme_bw()




# PCA per Chrom -----------------------------------------------------------

colnames(pca_table)
setnames(pca_table, colnames(pca_table), colnames(pca_table))

pca_table <- as.data.table(pca_table)

