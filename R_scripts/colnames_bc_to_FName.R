library(tidyverse)
library(data.table)

# r_thresh <- "10"

# vaf_n8 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_latest_VAFs/N7_adip_10_VAF_filt_avemed_025_distr40-60_025_20c_harmonized.txt')
# vaf_n7 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_latest_VAFs/N7_ery_10_VAF_filt_avemed_025_distr40-60_025_20c_harmonized.txt')
# vaf_n5 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_latest_VAFs/N7_naive_10_VAF_filt_avemed_025_distr40-60_025_20c_harmonized.txt')

vaf_n8 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_adip_VAFs_5.txt')
vaf_n7 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_ery_VAFs_5.txt')
vaf_n5 <- fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_naive_VAFs_5.txt')


# bc_n8 <-fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/filename_v_barcode_N8.txt', col.names = c("file", "bc"), header = F)
bc_n7 <-fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/filename_v_barcode_N7.txt', col.names = c("file", "bc"), header = F)
# bc_n5 <-fread('D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/filename_v_barcode_N5.txt', col.names = c("file", "bc"), header = F)
# bc_n8$file <- as.character(gsub("\\.fastq", "", bc_n8$file))
bc_n7$file <- as.character(gsub("\\.fastq", "", bc_n7$file))
# bc_n5$file <- as.character(gsub("\\.fastq", "", bc_n5$file))



colnames(vaf_n8)[2:ncol(vaf_n8)] <- sapply(colnames(vaf_n8)[2:ncol(vaf_n8)], function (x) {
  as.character(bc_n7[which(bc_n7$bc == x), 1])
})

colnames(vaf_n7)[2:ncol(vaf_n7)] <- sapply(colnames(vaf_n7)[2:ncol(vaf_n7)], function (x) {
  as.character(bc_n7[which(bc_n7$bc == x), 1])
})
colnames(vaf_n5)[2:ncol(vaf_n5)] <- sapply(colnames(vaf_n5)[2:ncol(vaf_n5)], function (x) {
  as.character(bc_n7[which(bc_n7$bc == x), 1])
})


fwrite(vaf_n8, 'D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_adip_VAFs_5_latest.txt', sep = "\t")
fwrite(vaf_n7, 'D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_ery_VAFs_5_latest.txt', sep = "\t")
fwrite(vaf_n5, 'D:/Virtualbox/Ubuntu/Shared Folder/Single_cell/GTEx comparisons/Clean_filtered_VAFs/inner_joined_N7_naive_VAFs_5_latest.txt', sep = "\t")
