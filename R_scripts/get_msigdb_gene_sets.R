library(msigdbr)
library(data.table)
library(tidyverse)


gene_sets <- c("M26610", "M18473", "M26496", "M18553", "M26833", "M26373", "M26298", "M17888")
genes_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF" ) %>% dplyr::filter(gs_id %in% gene_sets) #%>% dplyr::select(c("entrez_gene", "gene_symbol"))
colnames(genes_df)[1] <- "ENTREZID"
genes_df$ENTREZID <- as.character(genes_df$ENTREZID)

fwrite(genes_df, "D:/Virtualbox/Ubuntu/Shared Folder/Entcheva/all_ion_channels.txt", sep = "\t")
