library(tidyverse)
library(data.table)

exp_path <- "D:/Virtualbox/Ubuntu/Shared Folder/N_Lee_Nebraska/fC_DE/GTF-Stringtie/"
pattern <- ".gtf"
exp_files <- list.files(path = exp_path, pattern = pattern)

s <- data.table(TranscriptID = character(), GeneID = character(), GeneName = character())
# load in the readcounts
for (i in 1:length(exp_files)) {
  # select the necessary columns, including gene position and TPM value
  # fix sample name formatting, and remove genes without an annotation
  a <- fread(paste0(exp_path, exp_files[i]), skip = 2, sep = "\t", select = c(3, 9), col.names = c("type", "info")) %>%
    filter(grepl("*reference_id", info)) %>% filter(grepl("transcript", type)) %>% select(info) %>% 
    separate(info, into = c("GID", "TID", "TranscriptID", "GeneID", "GeneName", "cov", "FPKM", gsub(".gtf", "", paste0("TPM_", exp_files[i]), fixed = TRUE)), sep = ";")
  a <- data.frame(lapply(a, function(x) gsub(".*\\s\"{1}", "", x)))
  
  a<- data.frame(lapply(a, function(x) gsub("\"", "", x))) %>% select(3, 4, 5, 8)
  a$TPM <- as.numeric(a$TPM)
  s <- full_join(s, a, by = c("TranscriptID", "GeneID", "GeneName"))
  s <- s %>% select(-ncol(s))
    #mutate(sample = gsub(x = f, pattern = pattern, replacement = '')) %>% # remove sample name suffix specified in "pattern"
    #filter(`Gene Name` != '-')
  print(colnames(a))
}

fwrite(s, "TPM_counts.txt", sep = "\t")
