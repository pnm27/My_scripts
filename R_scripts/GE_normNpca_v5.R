library(tidyverse)
library(data.table)
library(ggplot2)


# setwd("~/Documents/OtherPIs/EEntcheva/QTL/")
# ge = fread('Heart_LV_female_GE.txt')
setwd("D:/Virtualbox/Ubuntu/Shared Folder/")
ge = fread('Heart_LV_female_GE_2.txt')


###REMOVE GAPDH
###retain only genes that are expressed at TPM >=1 in at least 20% of the samples
###from liam
# get number of samples

ge_long <- ge %>% gather(sample, TPM, -gene_id)
num_samples <- length(levels(factor(ge_long$sample)))
ge_long <- ge_long %>% group_by(gene_id) %>% mutate(count_zero = length(which(TPM < 1)),
                                            perc_zero = count_zero / num_samples)

# remove duplicate entries for each gene if necessary
df <- ge_long %>% filter(perc_zero < 0.8) %>% 
  select(sample, TPM, gene_id) %>% 
  group_by(sample, gene_id) %>% 
  # filter(row_number() == 1) %>%
  spread(sample, TPM)
###quantile transform (see Liam's build GE script)


# quantile normalize gene expression values
# (code adapted from MatrixEQTL sample code)
df_w <- df[-1]

for(gene in 1:nrow(df_w)) {
  mat = df_w[gene,]
  mat = apply(mat, 1, rank, ties.method = "average")
  mat = qnorm(mat / (ncol(df_w) + 1))
  df_w[gene,] = t(mat)
}

df_w <- cbind(df$gene_id, df_w)
colnames(df_w)[1] <- "gene_id"
fwrite(df_w, "Heart_LV_GE_q_transf.txt", sep = "\t") 
# Remove temp created variables
# rm(ge, ge_n, ge_df, ge_clean)

# build up PC
pca <- prcomp(na.omit(df_w[,-1]))
pc <- pca$rotation
png("PCAs_GE.png")
screeplot(pca, npcs = 20)
dev.off()
pc <- t(pc) %>% as.data.frame()
pc$id = row.names(pc)
pc = pc[,c(ncol(pc),1:(ncol(pc)-1))]
#selectNumberPCs,3 below
fwrite(pc[1:15,], "Heart_LV_cov.txt", sep = "\t")


qtls = fread("output/ReQTL_ee__all_ReQTLs.txt") %>% filter(FDR<0.0005)
qtls_norecipr <- qtls[qtls$SNP!=qtls$gene,]
top10000 = head()
write.table(qtls_norecipr, "ee_GE_QTLs005.txt", sep = "\t")
