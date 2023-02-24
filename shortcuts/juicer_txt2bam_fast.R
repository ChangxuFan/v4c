library(dplyr)
options(stringsAsFactors = F)
options(scipen = 999)
source("/bar/cfan/R_packages/v4c/R/juicer_txt2bam.R")
source("/bar/cfan/R_packages/v4c/R/juicer_sam2bam.R")
args <- commandArgs(trailingOnly=TRUE)

juicer.txt2bam.2(juicer.txt = args[1], mapq = args[2], genome = args[3], threads = args[4])
