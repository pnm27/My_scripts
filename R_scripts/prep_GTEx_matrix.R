library(data.table)
library(tidyverse)

a <- fread("D:/Virtualbox/Ubuntu/Shared Folder/Entcheva/Full-data.txt")
colnames(a)[1] <- "genes"

d <- a %>% distinct_at(vars(genes), .keep_all = T)
b <- a[duplicated(a$genes), 1]

d <- d[!which(genes %in% b$genes),]

g_list <- fread("D:/Virtualbox/Ubuntu/Shared Folder/Hs_allchr_MT.txt")
rm(b)
colnames(g_list)[2:5] <- c("genes", "start", "end", "chr")
combo <- inner_join(g_list, d, by = "genes")

check <- combo[duplicated(combo$genes),2]
check <- check %>% distinct()
combo <- combo[!which(genes %in% check$genes),]
combo <- combo %>% select(-`Gene stable ID`)
combo <- combo %>% filter(chr != "MT") %>% filter(chr != "X") %>% filter(chr != "Y")
combo <- na.omit(combo)
pos <- combo[, c(4, 2, 3)]



fwrite(combo, "D:/Virtualbox/Ubuntu/Shared Folder/Entcheva/GTEx_GE_no_dup_no_na.txt", sep = "\t")
fwrite(pos, "D:/Virtualbox/Ubuntu/Shared Folder/Entcheva/GTEX_pos_no_na.txt", sep = "\t")
