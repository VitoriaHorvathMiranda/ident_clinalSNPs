library(data.table)
library(tidyverse)

POP_09_10 <- fread("/dados/time_clines/output/SNPs/joined0910.tsv")
POP_09_102L <- POP_09_10[CHROM == "2L"]
rm(POP_09_102L)
POP_09_10X <- POP_09_10[CHROM == "X"]
POP_09_10[POS == 6398710]

n_snps <- POP_09_10[, nrow(.SD), by = position2]
n_snps[V1 < 5]

POP_09_10X[, position2 := paste(CHROM, POS, sep = ":")]
POP_09_10X <- POP_09_10X[, !c("CHROM", "POS")]
typeof(POP_09_10$freq)
POP_09_10[depth < 10,]

POP_09_10 <- POP_09_10[freq != "."][depth != 0]
#13312507
twofreqs <- POP_09_10[freq %like% ",",]
twoalts <- POP_09_10[!(ALT %like% ","),]

POP_09_10[ALT %like% ",",]

merge.data.table(twofreqs, twoalts, by = c("POS", "CHROM"))
anti_join(twoalts, twofreqs, by = c("POS", "CHROM"))



POP_09_10[, freq := as.double(freq)]


head(POP_09_10)
typeof(POP_09_10$freq)
typeof(POP_09_10$latitude)
typeof(POP_09_10$NE)

POP_09_10[is.na(freq), ]
POP_09_10[NE == 0]

#114480

POP_09_10[POS == 20242]


group_nest_dt <- function(dt, ..., .key = "data"){
  stopifnot(is.data.table(dt))
  by <- substitute(list(...))
  dt <- dt[, list(list(.SD)), by = eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

nested_snps <- POP_09_10X[,group_nest_dt(.SD, position2)]

nested_snps[524277,]

test <- nested_snps$data %>% pluck(524277) %>% as_tibble() 
test <- test %>%
  mutate(freq = as.double(freq))


glm(freq~latitude, 
    weights = NE,
    data = test2,
    family = binomial())


tiny_nested_snps <- nested_snps[position2 == "2R"]

POP_09_10[position2 == "2R:6398710", ]

test2 <- tiny_nested_snps$data %>% pluck(1103)

nested_snps[, models := purrr::map(data, ~ glm(freq~latitude, 
                                               weights = NE,
                                               data = .x,
                                               family = binomial()))]
?glm
#drop data column 
tiny_nested_snps <- tiny_nested_snps[, !c("data")]

#gets lat coefficient
nested_snps[, coefficients := purrr::map_dbl(models, ~coef(.x) %>% pluck("latitude")) ]

#gets p-value
nested_snps[, p_value := purrr::map_dbl(models, 
                                        ~summary(.x) %>% 
                                          pluck("coefficients") %>% pluck(8))]




