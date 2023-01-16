library(data.table)
library(tidyverse)
library(qvalue)


#11994070 - 2L
pop_test <- fread("/dados/time_clines/output/SNPs/joined97.tsv")

pop_test_lite <- pop_test[chrom == "2L",]
rm(pop_test)
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
#for unnesting
unnest_dt <- function(dt, col, id){
  stopifnot(is.data.table(dt))
  by <- substitute(id)
  col <- substitute(unlist(col, recursive = FALSE))
  dt[, eval(col), by = eval(by)]
}

nested_test <- pop_test_lite[depth>0,][, group_nest_dt(.SD, position2)]

#test[10000:10010,][position2 == "2L:622893", ] %>% pluck(2)

nested_test_filtered <- nested_test[,n_pops := purrr::map_dbl(data, nrow)][
  n_pops>2,][,!c("n_pops"), with = FALSE]

nested_test_filtered_lite <- nested_test_filtered[1:1000,] 
nested_test_filtered[, models := purrr::map(data, ~ glm(freq~latitude, 
                                    weights = NE,
                                    data = .x,
                                    family = binomial()))]

nested_test_filtered[, coefficients := purrr::map_dbl(models, ~coef(.x) %>% pluck("latitude")) ]
nested_test_filtered[, p_value := purrr::map_dbl(models, 
                                         ~summary(.x) %>% 
                                           pluck("coefficients") %>% pluck(8))]
nested_test_filtered <- nested_test_filtered[, !c("models")]

nested_test_filtered <- nested_test_filtered[, !c("data")]

typeof(nested_test_filtered)
class(nested_test_filtered[1])


unnest_dt(nested_test_filtered_lite,
          col = data,
          id = list(position2, coefficients, p_value))


# only using tidyverse ----------------------------------------------------
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









p_values <- nested_test_filtered$p_value
qobj <- qvalue(p = p_values)


hist <- hist(qobj)
ggsave("hist_test.png")
plot(qobj)

#Estimating a p-value cut-off for a given false discovery rate level
p_value_cutoff <- max(qobj$pvalues[qobj$qvalues <= 0.1])

#Estimating the false discovery rate for a given p-value cut-off
max(qobj$qvalues[qobj$pvalues <= 0.01]) 

significant_snps <- nested_test_filtered[p_value <= p_value_cutoff,]
head(nested_test_filtered, n = 10)

class(significant_snps)
typeof(significant_snps)

fwrite(significant_snps, "fwrite_test.tsv", sep = "\t")

class(significant_snps$position2)
class(significant_snps$coefficients)
class(significant_snps$p_value)

typeof(significant_snps$position2)
typeof(significant_snps$coefficients)
typeof(significant_snps$p_value)
