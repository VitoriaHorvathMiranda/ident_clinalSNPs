library(tidyverse)
library(data.table)

inv_freq <- fread("/dados/vitoria/invertions_freqs/inversion.freq")

inv_freq <- melt.data.table(inv_freq, measure.vars =  2:15,
                variable.name = "population",
                value.name = "freq")

seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)

complete <- inv_freq %>%
  inner_join(pop_info, by = "population") %>%
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  ))
  
complete %>%
  filter(collection_year != "2017") %>%
  ggplot(aes(x = latitude, y = freq)) +
  geom_point(aes(color = collection_year), alpha = 0.5) +
  geom_smooth(aes(color = collection_year), method = 'lm', alpha = 0.2) +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.07) +
  facet_wrap(~ Inv, nrow = 3) +
  theme_light() +
  labs(title = "Inversions Frequencies") +
  scale_fill_discrete(name = "Collection Year")


# 1º Call -----------------------------------------------------------------

inv_freq <- fread("/dados/vitoria/invertions_freqs/PoolSNP_mincount5_minfreq0.005_cov10_inversion.freq")

inv_freq <- melt.data.table(inv_freq, measure.vars =  2:15,
                            variable.name = "population",
                            value.name = "freq")

seq_metadata <- read_tsv("/home/vitoria/clinas/seq_metadata.tsv")
pop_info <- seq_metadata %>%
  select(population, collection_year, latitude)

complete <- inv_freq %>%
  inner_join(pop_info, by = "population") %>%
  mutate(collection_year = case_when( #treats collection_year as character
    collection_year == 1997 ~ "1997", 
    collection_year == 2009 ~ "2009/2010", #and join 2009 and 2010
    collection_year == 2010 ~ "2009/2010",
    collection_year == 2017 ~ "2017"
  ))

typeof(complete$collection_year)

complete %>%
  filter(collection_year != "2017") %>%
  ggplot(aes(x = latitude, y = freq)) +
  geom_point(aes(color = collection_year), alpha = 0.5) +
  geom_smooth(aes(color = collection_year), method = 'lm', alpha = 0.2) +
  geom_text(aes(label = population), size = 1.5, nudge_y = 0.07) +
  facet_wrap(~ Inv, nrow = 3) +
  theme_light() +
  labs(title = "Inversions Frequencies") +
  scale_fill_discrete(name = "Collection Year")



