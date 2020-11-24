library(tidyverse)
library(data.table)
library(readxl)
library(plyr)

a <- fread("../Bukrinsky/MockNef_TPMs_short.txt")

b <- fread("../Hs_allchr_MT.txt")              # Biomart input of Gene stable ID, Gene name, Gene start (bp), Gene end (bp), Chromosome/scaffold name

d <- inner_join(a, b, by = c("GID" = "Gene stable ID"))

d <- d %>% select(-Gene_name, -`Gene start (bp)`, -`Gene end (bp)`)
d <- d[,c(1, 18, 2:17, 19)]
colnames(d)[2] <- "Gene_name"

e = copy(d)
setorder(e, Gene_name)
test <- as.data.frame(table(e$Gene_name)) # Frequency table of Gene names
table(table(e$Gene_name)) # shows how many genes are unique, how many are present twice, thrice....

e <- e[!duplicated(e$Gene_name), ]  # Removes the duplicated Genes (keeps the one with highest as we sorted it using setorder function)

fwrite(e, "../Bukrinsky/MockNef_TPMs_unique_highest.txt", sep = "\t")

f <- ddply(d, "Gene_name", numcolwise(sum))  # Sum duplicated values column-wise, split at the Gene_name column

fwrite(f, "../Bukrinsky/MockNef_TPMs_unique_sum.txt", sep = "\t")
