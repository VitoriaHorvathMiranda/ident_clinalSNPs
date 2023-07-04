library(data.table)

metadata <- fread("/dados/time_clines/data/meta/seq_metadata.tsv")

depths <- fread("/dados/time_clines/data/seqs/calls/freqs_and_depths.tsv")

depths <- depths[depth != ".", ]
depths <- depths[, depth := as.double(depth)]

typeof(depths$depth)

I_info <- metadata[, .(population, n_females, latitude)]
I_info <- I_info[, n_chrom := n_females*2][,!c("n_females")]

merged <- merge.data.table(I_info, depths, by = "population")

freqs <- fread("/dados/time_clines/data/seqs/calls/freqs.tsv")

merge.data.table(freqs, merged, allow.cartesian=TRUE)





