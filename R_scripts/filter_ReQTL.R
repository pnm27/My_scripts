library(data.table)
library(tidyverse)

a <- fread("Heart_LV_female_combo_ReQTL.txt")

a_no_recip <- a[a$SNP!=a$gene,]
# gene_l <- fread("2Prashant_09232020/circ_and interesting_genes.txt")
# a_gene <- a_no_recip %>% filter(SNP %in% gene_l$gene | gene %in% gene_l$gene)

a_gene <- a_gene %>% filter(grepl("^[^MT-]", SNP, perl = T)) %>%  filter(grepl("^[^MT-]", gene, perl = T))
a_no_recip <- a_no_recip %>% filter(grepl("^[^MT-]", SNP, perl = T)) %>%  filter(grepl("^[^MT-]", gene, perl = T))

fwrite(a_no_recip, "Heart_LV_female_combo_ReQTL_filt.txt", sep = "\t")
