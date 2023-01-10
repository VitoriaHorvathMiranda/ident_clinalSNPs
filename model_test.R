library(data.table)
library(tidyverse)
library(broom)

#11994070 - 2L
pop_test <- fread("/dados/time_clines/output/SNPs/joined97.tsv")

pop_test_lite <- pop_test[chrom == "2L",]

pop_test_lite[, position2 := paste(chrom, position, sep = ":")]

pop_test_lite <- pop_test_lite[, !c("chrom", "position")] 



#tidyverse and data.table try:
#data.table function for nesting:
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

nested_test <- pop_test_lite[depth>0,][, group_nest_dt(.SD, position2)]

#test[10000:10010,][position2 == "2L:622893", ] %>% pluck(2)

nested_test_filtered <- nested_test[,n_pops := purrr::map_dbl(data, nrow)][
  n_pops>2,][,!c("n_pops"), with = FALSE]


nested_test_filtered[, models := purrr::map(data, ~ glm(freq~latitude, 
                                    weights = NE,
                                    data = .x,
                                    family = binomial()))]

nested_test_filtered[, coefficients := purrr::map_dbl(models, ~coef(.x) %>% pluck("latitude")) ]
nested_test_filtered[, p_value := purrr::map_dbl(models, 
                                         ~summary(.x) %>% 
                                           pluck("coefficients") %>% pluck(8))]

#fited_models <- 
magica_table <- pop_test_lite[depth>0,] %>%
  group_by(position2) %>% 
  nest() %>% 
  mutate(n_pops = map_dbl(data, nrow)) %>%
  filter(n_pops > 2) %>%
  select(!n_pops) %>%
  mutate(models = map(data, ~ glm(freq~latitude, weights = NE,
                data = .x,
                family = binomial()))) 


magica_table %>%
  mutate(coefficients = map_dbl(models, ~coef(.x) %>% pluck("latitude"))) %>%
  mutate(p_value = map_dbl(models, 
                       ~summary(.x) %>% 
                         pluck("coefficients") %>% pluck(8))) 
  


magica_table %>%
  mutate(n_pops = map_dbl(data, nrow)) %>%
  filter(n_pops > 2) %>%
  select(!n_pops)




GLM_TABLE <- pop_test %>%
  group_by(position2) %>% 
  nest() %>% 
  mutate(n_pops = map_dbl(data, nrow)) %>%
  filter(n_pops > 2) %>%
  select(!n_pops) %>%
  mutate(models = map(data, ~ glm(freq~latitude, weights = NE,
                                  data = .x,
                                  family = binomial()))) %>%
  mutate(coefficients = map_dbl(models, ~coef(.x) %>% pluck("latitude"))) %>%
  mutate(p_value = map_dbl(models, 
                           ~summary(.x) %>% 
                             pluck("coefficients") %>% pluck(8))) 






