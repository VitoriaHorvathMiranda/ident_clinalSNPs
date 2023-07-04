library(data.table)
library(tidyverse)

p_value_97_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_97_mincount5_minfreq0.005_cov10.tsv")
NE_table <- fread("/dados/time_clines/output/SNPs/joined97_mincount5_minfreq0.005_cov10.tsv")

p_value_0910_cov10 <- fread("/dados/time_clines/output/SNPs/p-values_0910_mincount5_minfreq0.005_cov10.tsv")
NE_table_0910 <- fread("/dados/time_clines/output/SNPs/joined0910_mincount5_minfreq0.005_cov10.tsv")

seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)



# 2L test -----------------------------------------------------------------

p_value_97_lite <- p_value_97_cov10[position2 %like% "2L",]
NE_table_lite <- NE_table[CHROM == "2L",]

NE_table_lite[, position2 := paste(CHROM, POS, sep = ":")]

no_freq_lite <- NE_table_lite[, .(mean = mean(freq)), by = c("position2")][
  mean == 0,]

p_value_97_no0freq <- p_value_97_lite[!no_freq_lite, on = "position2"]

# all chrom 97 ---------------------------------------------------------------

NE_table[, position2 := paste(CHROM, POS, sep = ":")]

#no_freq mean==0 -> 187259
#no_freq mean <=0.005 -> 260257
no_freq<- NE_table[, .(mean = mean(freq)), by = c("position2")][
  mean <= 0.005,]

p_value_97_no0freq <- p_value_97_cov10[!no_freq, on = "position2"]


p_value_97_no0freq %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "1997",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

# all chorm 0910 ----------------------------------------------------------

NE_table_0910[, position2 := paste(CHROM, POS, sep = ":")]

no_freq_0910 <- NE_table_0910[, .(mean = mean(freq)), by = c("position2")][
  mean == 0,]

p_value_0910_no0freq <- p_value_0910_cov10[!no_freq_0910, on = "position2"]


p_value_0910_no0freq %>%
  ggplot() +
  geom_histogram(aes(p_value)) +
  labs(x = "p-value", title = "2009/2010",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()


# No unique snps 97 ----------------------------------------------------------

NE_table[, freq_status := 1*(freq!= 0)]

unique_snps_97 <- NE_table[, .(freq_status_sum = sum(freq_status)),
                           by = "position2"][freq_status_sum ==0,]

p_value_97_no_unique <- p_value_97_cov10[!unique_snps_97, on = "position2"]

p_value_97_no_unique %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "1997",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

p_value_97_only_unique <- p_value_97_cov10[unique_snps_97, on = "position2"]

p_value_97_only_unique %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "1997 _unique",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

unique_NE <- NE_table[unique_snps_97, on = "position2"]

unique_NE[, .(mean_freq = mean(freq)), by = "position2"] %>%
  ggplot() +
  geom_histogram(aes(mean_freq))

unique_NE %>%
  ggplot() +
  geom_histogram(aes(freq, fill = population)) +
  facet_wrap(~ population)


# No unique 0910 ----------------------------------------------------------

NE_table_0910[, freq_status := 1*(freq!= 0)]

unique_snps_0910 <- NE_table_0910[, .(freq_status_sum = sum(freq_status)),
                           by = "position2"][freq_status_sum ==1,]

p_value_0910_no_unique <- p_value_0910_cov10[!unique_snps_0910, on = "position2"]

p_value_0910_no_unique %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "2009/2010",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()

p_value_0910_only_unique <- p_value_0910_cov10[unique_snps_0910, on = "position2"]

p_value_0910_only_unique %>%
  ggplot() +
  geom_histogram(aes(p_value))+
  labs(x = "p-value", title = "2009/2010 _unique",
       caption = "min-freq = 0.005, min-cov = 10, min-count = 5") +
  theme_minimal()


unique_NE_0910 <- NE_table_0910[unique_snps_0910, on = "position2"]

unique_NE_0910[, .(mean_freq = mean(freq)), by = "position2"] %>%
  ggplot() +
  geom_histogram(aes(mean_freq))

unique_NE_0910 %>%
  ggplot() +
  geom_histogram(aes(freq))

unique_NE_0910[, .(mean_depth = mean(depth)), by = "position2"] %>%
  ggplot() +
  geom_histogram(aes(mean_depth))

# change in unique snps? --------------------------------------------------

merged_unique <- merge.data.table(unique_NE, unique_NE_0910, all=TRUE)

merged_unique_nofre0 <- merged_unique[freq != 0,]

merged_unique_nofre0 %>%
  ggplot() +
  geom_point(aes(x=latitude, y=freq), size = 0.01)


##filtra snps presententes nas duas épocas
merged_unique_nofre0 <- merged_unique_nofre0[, !c("freq_status_sum"), with=FALSE]
both_periods_unique <- merged_unique_nofre0[, .(freq_status_sum = sum(freq_status)), by = "position2"][
  freq_status_sum >= 2,
  ]

filtered_merged <- merged_unique_nofre0[both_periods_unique, on = "position2"]

ggplot(filtered_merged) +
  geom_point(aes(x = latitude, y = freq), size = 0.01)


filtered_merged[, year := ifelse(population %like% "97", "Y_1997", "Y_2009_2010")]
filtered_merged[, collection_area := case_when(latitude > 40 ~ "high",
                                               latitude > 30 & latitude < 40 ~ "medium",
                                               latitude < 30 ~ "low")]

dcasted_table <- dcast.data.table(filtered_merged, collection_area + position2 ~ year, value.var = "freq")

changed_snps <- dcasted_table[!(is.na(Y_1997)),][!(is.na(Y_2009_2010))]

dcasted_table[is.na(Y_1997) | is.na(Y_2009_2010),]

changed_snps[, difference := abs(Y_1997 - Y_2009_2010)]

changed_snps %>%
  ggplot() +
  geom_histogram(aes(difference))

changed_snps[difference >= 0.2,] %>%
  ggplot() +
  geom_histogram(aes(difference, fill = collection_area))


changed_snps %>%
  ggplot() +
  geom_point(aes(x = collection_area, y = difference), size = 0.05)


# dist_unique -------------------------------------------------------------

dist_unique <- tibble(snps = c("total", "zero_pops", "one_pop"),
       Y_1997 = c(187259+297701, 187259, 297701),
       Y_2009_2010 = c(28625+173978, 28625, 173978))

tidy_dist_unique <- dist_unique %>%
  pivot_longer(cols = c(Y_1997, Y_2009_2010), 
               names_to = "year",
               values_to = "n_snps")


tidy_dist_unique %>%
  filter(snps != "total") %>%
  ggplot(aes(x = year, y = n_snps)) +
  geom_bar(aes(fill = snps), stat = "identity", alpha = 0.6) +
  scale_fill_discrete(labels = c("one_pop" = "one population only",
                                 "zero_pops" = "not present in that period")) +
  labs(y = "nº of snps") +
  scale_x_discrete(labels = c("1997", "2009/2010")) +
  theme_bw()
  
  
  
  
  


many97_one0910 <- NE_table[p_value_0910_only_unique, on = "position2"]

many97_one0910[, .(freq_status_sum = sum(freq_status)),
               by = "position2"][freq_status_sum > 1,]

many0910_one97 <- NE_table_0910[p_value_97_only_unique, on = "position2"]
many0910_one97[, .(freq_status_sum = sum(freq_status)),
               by = "position2"][freq_status_sum >1]



