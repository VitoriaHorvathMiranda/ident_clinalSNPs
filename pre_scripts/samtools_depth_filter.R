library(data.table)

positions <- fread("/dados/time_clines/data/seqs/calls/called_snps.tsv")

depths <- fread("/dados/vitoria/all.samtoolsDepth_10000Lines.tsv")

oldnames <- paste0("V", 1:17)
newnames <- c("CHROM", "POS", "ESC97", "CMD97B", "WVT97", "MFL97", "RVA97", "SVT09",
              "CMD10", "SNC10", "JFL10", "CMA97", "HMA09", "HMA17",
              "MCT09", "MCT17", "MCT97")

setnames(depths, oldnames, newnames)       
setkey(depths, CHROM, POS)
setkey(positions, CHROM, POS)

test <- data.table::merge.data.table(positions, depths)

head(test)
head(positions)

test <- test[, data.table::melt(.SD, measure.vars = 3:17,
                        variable.name = "population",
                        value.name = "depth")]
