#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()

parser$add_argument('--SnapeFile', '-sf', help = 'snape filtered file')
parser$add_argument('--DepthFile', '-df', help= 'samtools depth file')

xargs<- parser$parse_args()

xargs$sf