
library(data.table)

args <- c("/dados/time_clines/output/SNPs/CMA97_L002.NE.tsv", 
          "/dados/time_clines/output/SNPs/HMA09_L002.NE.tsv",
          "/dados/time_clines/output/SNPs/MCT97_L002.NE.tsv",
         "/dados/time_clines/output/SNPs/HMA17_L002.NE.tsv",
         "/dados/time_clines/output/SNPs/MCT09_L002.NE.tsv",
          "/dados/time_clines/output/SNPs/MCT17_L002.NE.tsv")



output <- vector("list", length(args))
for (i in seq_along(args)) {
  output[[i]] <- fread(file = args[[i]])
}

full_merge <- function(...) {
  merge.data.table(..., all = TRUE)
}

ALL_POP <- Reduce(full_merge, output)

POP_97 <- ALL_POP[population %like% "97"]

POP_09_10 <- ALL_POP[population %like% "09" | population %like% "10"]

POP_17 <- ALL_POP[population %like% "17"]

