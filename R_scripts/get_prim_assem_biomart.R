library(biomaRt)
library(data.table)
library(tidyverse)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
f1 <- fread("D:/Virtualbox/Ubuntu/Shared Folder/Archive (1)/SRR3192364_gene_abund.tab")
a <- colnames(f1)[1]
colnames(f1)[1] <- "Gene_ID"

f2 <- fread("D:/Virtualbox/Ubuntu/Shared Folder/Archive (1)/SRR3192365_gene_abund.tab")
colnames(f2)[1] <- "Gene_ID"

# If Biomart works

t<- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "external_gene_source"), filters = "chromosome_name",
          values = c(1:22, "X", "Y"),
          mart= mart)
t <- t %>% filter(external_gene_source == "HGNC Symbol")



# Otherwise
f3 <- fread("D:/Virtualbox/Ubuntu/Shared Folder/Useful R Scripts/Hs_gene_id_name_source.txt")
# 4 possible sources for gene names: Clone-based (Ensembl) gene, HGNC Symbol, miRBase NCBI gene (formerly Entrezgene), RFAM
f3 <- f3 %>% filter(`Source of gene name` == "HGNC Symbol")
colnames(f3)[1] <- "Gene_ID"



combo_1 <- inner_join(f1, f3)
combo_2 <- inner_join(f2, f3)
colnames(combo_1)[1] <- a
colnames(combo_2)[1] <- a

combo_1 <- combo_1 %>% select(-`Gene Name`, -`Source of gene name`)
combo_2 <- combo_2 %>% select(-`Gene Name`, -`Source of gene name`)



fwrite(combo_1, "D:/Virtualbox/Ubuntu/Shared Folder/Archive (1)/SRR3192364_gene_abund_prim.txt", sep = "\t")
fwrite(combo_2, "D:/Virtualbox/Ubuntu/Shared Folder/Archive (1)/SRR3192365_gene_abund_prim.txt", sep = "\t")
