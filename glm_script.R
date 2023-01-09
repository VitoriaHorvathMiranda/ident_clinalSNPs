library(argparse)
library(tidyverse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "computes a glm for each snp")
parser$add_argument('--timePops', '-tPops',
                    help = 'outputs from join_time_pops_script.R')
parser$add_argument('output', '-o',
                    help = 'table with glm, coefficient and p-value')
xargs<- parser$parse_args()

joined_pops <- fread(file = xargs$timePops)
joined_pops[, position2 := paste(chrom, position, sep = ":")]
joined_pops <- joined_pops[, !c("chrom", "position")]

GLM_TABLE <- joined_pops %>%
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

glm_table <- xargs$output

saveRDS(GLM_TABLE, file = glm_table)




