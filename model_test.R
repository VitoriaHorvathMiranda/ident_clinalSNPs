library(data.table)
library(tidyverse)
library(broom)


pop_test <- fread("/dados/time_clines/output/SNPs/joined97.tsv")

pop_test[, position2 := paste(chrom, position, sep = ":")]

pop_test <- pop_test[order(position2),] 

pop_test <- pop_test[, !c("chrom", "position")] 

pop_test_lite <- pop_test[0:500,]


#tidyverse try:
#fited_models <- 
magica_table <- pop_test_lite %>%
  group_by(position2) %>% 
  nest() %>% 
  mutate(models = map(data, ~ glm(freq~latitude, weights = NE,
                data = .x,
                family = binomial()))) 


magica_table %>%
  mutate(n_pops = map_dbl(data, nrow)) %>%
  filter(n_pops > 2) %>%
  select(!n_pops) %>%
  mutate(coefficients = map_dbl(models, ~coef(.x) %>% pluck("latitude"))) %>%
  mutate(p_value = map_dbl(models, 
                       ~summary(.x) %>% 
                         pluck("coefficients") %>% pluck(8))) 
  


magica_table %>%
  mutate(n_pops = map_dbl(data, nrow)) %>%
  filter(n_pops > 2) %>%
  select(!n_pops)











