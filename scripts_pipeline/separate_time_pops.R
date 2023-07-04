#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Separate pops from different times")
parser$add_argument('--NEtable', '-NE', help = 'outputs from NE.R')
parser$add_argument('--output97', '-o97', 
                    help= 'table with all pops from 1997')
parser$add_argument('--output0910', '-o0910', 
                    help= 'table with all pops from 2009|2010')
parser$add_argument('--output17', '-o17', 
                    help= 'table with all pops from 2017')
xargs<- parser$parse_args()
pops <- xargs$populations

#reads all pops
ALL_POP <- fread(file = xargs$NEtable)

#separate pops based on year
POP_97 <- ALL_POP[population %like% "97"]

POP_09_10 <- ALL_POP[population %like% "09" | population %like% "10"]

POP_17 <- ALL_POP[population %like% "17"]

otp97 <- xargs$output97
otp0910 <- xargs$output0910
otp17 <- xargs$output17

fwrite(POP_97, otp97, sep = "\t")
fwrite(POP_09_10, otp0910, sep = "\t")
fwrite(POP_17, otp17, sep = "\t")
