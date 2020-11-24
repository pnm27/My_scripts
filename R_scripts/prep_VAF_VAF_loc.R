library(tidyverse)
library(data.table)


r_thresh <- "5"
path = paste0('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/common_ReQTL/sc_ol_GTEx_sQTL_43tissues/ol_', r_thresh,'/')
exp_files <- list.files(path = path, pattern = '*ol_minR.*.txt')

df <- bind_rows(lapply(exp_files, function(f) {
  # select the necessary columns, including gene position and TPM value
  # fix sample name formatting, and remove genes without an annotation
  fread(paste0(path, f)) %>% 
    select(1, 2, 9, 25, 36)
}))

# df_n8 <- df[which(df$Donor == "N8"), 1:4] %>% distinct(SNV, .keep_all = T)
df_n7 <- df[which(df$Donor == "N7"), 1:4] %>% distinct(SNV, .keep_all = T)
# df_n5 <- df[which(df$Donor == "N5"), 1:4] %>% distinct(SNV, .keep_all = T)
colnames(vaf_n7)[1] <- "SNV"
colnames(vaf_n8)[1] <- "SNV"
vaf_n8 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/N7_min5_adip_VAF_filt.txt')
vaf_n7 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/N7_min5_ery_VAF_filt.txt')
vaf_n5 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/N7_min5_naive_VAF_filt.txt')


# For making new VAF files
vaf_n8_joined <- inner_join(vaf_n8, df_n7, by ="SNV") %>% select(1:ncol(vaf_n8))
vaf_n7_joined <- inner_join(vaf_n7, df_n7, by ="SNV") %>% select(1:ncol(vaf_n7))
vaf_n5_joined <- inner_join(vaf_n5, df_n7, by ="SNV") %>% select(1:ncol(vaf_n5))

# new VAF files
fwrite(vaf_n8_joined, paste0("inner_joined_N7_adip_VAFs_", r_thresh,".txt"), sep = "\t")
fwrite(vaf_n7_joined, paste0("inner_joined_N7_ery_VAFs_", r_thresh,".txt"), sep = "\t")
fwrite(vaf_n5_joined, paste0("inner_joined_N7_naive_VAFs_", r_thresh,".txt"), sep = "\t")

# For making new correlation file to plot in the script called Plot_ReQTLS_VAF_star.r (put in as in_file)
res_n8_joined <- inner_join(df_n7, vaf_n8, by ="SNV") %>% select(1:ncol(df_n7)) %>% mutate(gene = paste(geneList, gene_id, sep = "-")) %>% select(SNV, gene)
res_n7_joined <- inner_join(df_n7, vaf_n7, by ="SNV") %>% select(1:ncol(df_n7)) %>% mutate(gene = paste(geneList, gene_id, sep = "-")) %>% select(SNV, gene)
res_n5_joined <- inner_join(df_n7, vaf_n5, by ="SNV") %>% select(1:ncol(df_n7)) %>% mutate(gene = paste(geneList, gene_id, sep = "-")) %>% select(SNV, gene)

colnames(res_n8_joined)[1] <- "SNP"
colnames(res_n7_joined)[1] <- "SNP"
colnames(res_n5_joined)[1] <- "SNP"

# For making new correlation file
fwrite(res_n8_joined, paste0("res_N7_adip_VAFs_", r_thresh,".txt"), sep = "\t")
fwrite(res_n7_joined, paste0("res_N7_ery_VAFs_", r_thresh,".txt"), sep = "\t")
fwrite(res_n5_joined, paste0("res_N7_naive_VAFs_", r_thresh,".txt"), sep = "\t")
