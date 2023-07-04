library(tidyverse)
n_snps <- read_csv("n_snps.tsv")

n_snps %>%
  mutate(min_pop = as.integer(min_pop)) %>%
  ggplot(aes()) +
  geom_line(aes(x = min_pop, y = total_snps, linetype = year), color = "darkcyan") +
  geom_line(aes(x = min_pop, y = FDR_01, linetype = year), color = "springgreen4") +
  geom_line(aes(x = min_pop, y = FDR_005, linetype = year), color = "red4") +
  scale_y_log10()
  
n_snps %>%
  mutate(min_pop = as.integer(min_pop)) %>%
  ggplot(aes()) +
  #geom_line(aes(x = min_pop, y = total_snps, linetype = year), color = "darkcyan") +
  geom_line(aes(x = min_pop, y = FDR_01, linetype = year), color = "springgreen4") +
  geom_line(aes(x = min_pop, y = FDR_005, linetype = year), color = "red4")


n_snps %>%
  mutate(min_pop = as.integer(min_pop)) %>%
  ggplot(aes()) +
  geom_point(aes(x = min_pop, y = total_snps, shape = year), color = "darkcyan") +
  geom_point(aes(x = min_pop, y = FDR_01, shape = year), color = "springgreen4") +
  geom_point(aes(x = min_pop, y = FDR_005, shape = year), color = "red4") +
  scale_y_log10()

