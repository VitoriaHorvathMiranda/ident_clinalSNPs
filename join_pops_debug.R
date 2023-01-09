#!/usr/bin/env Rscript

library(argparse)
library(data.table)

#parse arguments
parser <- ArgumentParser(description= "Unites pops from different times")
parser$add_argument('--populations', '-pops',
                    help = 'outputs from script_NE.R',
                    nargs = '+')

parser$add_argument('--output97', '-o97', 
                    help= 'table with all pops from 1997')

parser$add_argument('--output0910', '-o0910', 
                    help= 'table with all pops from 2009|2010')

parser$add_argument('--output17', '-o17', 
                    help= 'table with all pops from 2017')

xargs<- parser$parse_args()

pops <- xargs$populations

pops

length(pops)

xargs$output97





