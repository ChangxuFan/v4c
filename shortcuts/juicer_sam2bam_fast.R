library(dplyr)
options(stringsAsFactors = F)
options(scipen = 999)
source("/bar/cfan/R_packages/v4c/R/juicer_sam2bam.R")
args <- commandArgs(trailingOnly=TRUE)
juicer.sam2bam(sam = args[1],
               genome = args[2], threads = args[3], space.to.tab = as.logical(args[4]))
