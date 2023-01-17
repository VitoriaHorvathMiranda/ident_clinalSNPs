library(argparse)
library(tidyverse)
library(data.table)
library(qvalue)

#parse arguments
parser <- ArgumentParser(description= "computes a glm for each snp")
parser$add_argument('--timePops', '-tPops',
                    help = 'outputs from join_time_pops_script.R')
parser$add_argument('--minPOP', '-minpop',
                    help = 'Minimum number of populations in which a SNP
                    must be present to be considered')
parser$add_argument('--output', '-o',
                    help = 'table with inclination coefficient and 
                    p-values for each SNP, .tsv')

xargs<- parser$parse_args()

joined_pops <- xargs$timePops

#reads data and joins chrom and position
joined_pops <- fread(file = joined_pops)
joined_pops[, position2 := paste(chrom, position, sep = ":")]
joined_pops <- joined_pops[, !c("chrom", "position")]

#data.table function for nesting:
group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

#filters calls with depth = 0
#nest data by snp
#filter out snps that were only called in two or less pops
minPop <- xargs$minPOP
nested_snps <- joined_pops[depth>=minPop,][,group_nest_dt(.SD, position2)][, 
                                       n_pops := purrr::map_dbl(data, nrow)][
                                       n_pops>2,][,
                                       !c("n_pops"), with = FALSE]
#runs glm for each snp
nested_snps[, models := purrr::map(data, ~ glm(freq~latitude, 
                                               weights = NE,
                                               data = .x,
                                               family = binomial()))]

#drop data column 
nested_snps <- nested_snps[, !c("data")]

#gets lat coefficient
nested_snps[, coefficients := purrr::map_dbl(models, ~coef(.x) %>% pluck("latitude")) ]

#gets p-value
nested_snps[, p_value := purrr::map_dbl(models, 
                                         ~summary(.x) %>% 
                                           pluck("coefficients") %>% pluck(8))]

#drops models column (it is too heavy)
nested_snps <- nested_snps[, !c("models")]


#saves
output_path <- xargs$output
fwrite(nested_snps, output_path, sep = "\t")

