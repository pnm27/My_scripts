###
# This script serves to retain only the known clusters to be present in those genes that have
# been annotated in the GTEx portal for having cis-correlation
# for all the three sets of samples

###

library(data.table)
library(tidyverse)
library(GenomicRanges)

setwd("D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/")
sample <- "N5"
r_thresh <- "5"

n8_c <- fread(paste0(sample, "_20_perind_numers.counts"))
colnames(n8_c)[1] <- "cluster"
col <- ncol(n8_c)
n8_c <- n8_c %>% separate(cluster, into = c("contig", "coord1", "coord2", "clust"), sep = ":", remove = F)
col <- ncol(n8_c)
#n8_c <- n8_c[, c(2:4, 1, 5:8923)]
colnames(n8_c)[1:4] <- c("gene", "chrom", "start", "end")
n8_gl <- GRanges(n8_c)
n8_c <- n8_c %>% select(-chrom, -start, -end, -clust)



my_list <- fread(paste0("common_ReQTL/sc_ol_GTEx_sQTL_43tissues/ol_", r_thresh, "/res_", sample, "_VAFs_", r_thresh, ".txt"))
my_list$gene <- gsub("-ENSG.*", "", my_list$gene)
my_list <- my_list %>% select(gene)
colnames(my_list) <- "ensembl_gene"
gene_list <- fread("D:/Virtualbox/Ubuntu/Shared Folder/RsQTL-master/RsQTL-master/data/gene_locations_hg38.txt")
new_list <- inner_join(my_list, gene_list)
#in_j <- inner_join(n8_c, gene_list, by = c("contig" = "chrom"))
# n8_c$gene <- ifelse(sapply(n8_c[, c(4, 1, 2)], function (x) ))

gl <- GRanges(new_list)

a <- suppressWarnings(findOverlaps(n8_gl, gl, type = "within", select = "last", ignore.strand = TRUE))
# table(a)
b <- as.data.frame(a)
# check <- as.data.frame(a)
# for (i in 1:nrow(n8_c)) {
  # n8_c$gene[i] <- gene_list[a[i], 1]  
# }

d <- cbind.data.frame(n8_c, b)
e <- d[!is.na(a), ]
e <- e %>% select(-a)
colnames(e)[1] <- ""

fwrite(e, paste0(sample, "_20reads_", r_thresh, ".counts"), sep = "\t")
