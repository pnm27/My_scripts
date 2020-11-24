library(data.table)
library(tidyverse)


setwd("../Single_cell/GTEx comparisons/Clean_latest_VAFs/")

a <- fread("N7_naive_10_VAF_filt_avemed_025_distr40-60_025_20c_harmonized.txt")

b <- list(colnames(a)[-1])

fwrite(b, "N7_naive_10_bc.txt", sep = "\t")
