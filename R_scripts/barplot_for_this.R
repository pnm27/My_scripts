library(data.table)
library(tidyverse)
library(gridExtra)
library(grid)

############################################################
#Input format for matrix should be: colnames as "SNV, sample1, sample2.....sampleN (doesn't depend on the naming style except at the simplification step
# i.e gather function)"
#with SNV having a format like 1:10000_A>T
############################################################


setwd('../Single_cell/Fig7_2Prashant/Fig7_2Prashant/')


sample_name = "N5"
#Input VAF file, split snv, deselect ref and alt and simplify colnames
vaf_data <- fread(paste0('../../Bar plots/', sample_name, '_10_VAF_matrix_CLEAN_10c_nExtrRmean.txt')) %>% separate(SNV, into = c("chr", "loc", "ref", "alt"))
vaf_data$loc <- as.integer(vaf_data$loc)
colnames(vaf_data) <- sapply(colnames(vaf_data), function (x) gsub("_wasp_AlignedByCoord_vW_filt_dedupped_sorted", "", x))

snv_set <- fread(paste0(sample_name, "_10_mono_2plot.txt"), col.names = "SNV", header = FALSE) %>% separate(SNV, into = c("chr", "loc", "ref", "alt"))
snv_set$loc <- as.integer(snv_set$loc)
snv_set2 <- fread(paste0(sample_name, "_10_skewed_2plot.txt"), col.names = "SNV", header = FALSE) %>% separate(SNV, into = c("chr", "loc", "ref", "alt"))
snv_set2$loc <- as.integer(snv_set2$loc)

vaf_data1 <- inner_join(vaf_data, snv_set, by = c("chr", "loc", "ref", "alt"))
vaf_data2 <- inner_join(vaf_data, snv_set2, by = c("chr", "loc", "ref", "alt"))

# Create gene list (for current use we have a customised set of 57 genes from literature) with colnames as "gene, chr, start, end"
# Input gene reference from biomart with 4 columns: "gene, chr, start, end"
gene_ref <- fread('../../../RsQTL-master/RsQTL-master/data/gene_locations_hg38.txt')
names(gene_ref) <- c("gene", "chr", "start", "end")
gene_ref <- gene_ref %>% filter(grepl("^[1-9XY]|^1[0-9]|^2[0-2]", chr))


# Order by chr and loc, create a new "formatted" location column, and rename it while removing the old columns (chr and loc)
vaf_data1 <- inner_join(vaf_data1, gene_ref, by = "chr") %>% filter(loc >= start & loc <= end) %>% distinct() %>% select(-start, -end)
vaf_data2 <- inner_join(vaf_data2, gene_ref, by = "chr") %>% filter(loc >= start & loc <= end) %>% distinct() %>% select(-start, -end)

# Annotate the VAFs with customised gene list, order by chr and loc, create a new "formatted" location column, and rename it, keep chr column for grouping
setorder(vaf_data1, chr, loc)
vaf_data1 <- vaf_data1 %>% unite("loc1", chr:loc, remove = T)
vaf_data1$loc1 <- sapply(vaf_data1$loc1, function (x) paste0("chr", x))
names(vaf_data1)[1] <- "loc"

setorder(vaf_data2, chr, loc)
vaf_data2 <- vaf_data2 %>% unite("loc1", chr:loc, remove = T)
vaf_data2$loc1 <- sapply(vaf_data2$loc1, function (x) paste0("chr", x))
names(vaf_data2)[1] <- "loc"


#Tidy the data: count for each bins, group by snvs
# If the bins have changed from <0.2, 0.2-0.4....then replace '4' with (# of bins-1)
vaf_data1 <- vaf_data1 %>% group_by(loc)
vaf_data1$`<0.2` <- rowSums(vaf_data1[, 4:(ncol(vaf_data1)-1)] < 0.2, na.rm = TRUE)
vaf_data1$`0.2-0.4` <- rowSums(vaf_data1[, 4:(ncol(vaf_data1)-2)] >= 0.2 & vaf_data1[, 4:(ncol(vaf_data1)-2)] < 0.4 , na.rm = TRUE)
vaf_data1$`0.4-0.6` <- rowSums(vaf_data1[, 4:(ncol(vaf_data1)-3)] >= 0.4 & vaf_data1[, 4:(ncol(vaf_data1)-3)] < 0.6, na.rm = TRUE)
vaf_data1$`0.6-0.8` <- rowSums(vaf_data1[, 4:(ncol(vaf_data1)-4)] >= 0.6 & vaf_data1[, 4:(ncol(vaf_data1)-4)] < 0.8, na.rm = TRUE)
vaf_data1$`>0.8` <- rowSums(vaf_data1[, 4:(ncol(vaf_data1)-5)] >= 0.8, na.rm = TRUE)

vaf_data2 <- vaf_data2 %>% group_by(loc)
vaf_data2$`<0.2` <- rowSums(vaf_data2[, 4:(ncol(vaf_data2)-1)] < 0.2, na.rm = TRUE)
vaf_data2$`0.2-0.4` <- rowSums(vaf_data2[, 4:(ncol(vaf_data2)-2)] >= 0.2 & vaf_data2[, 4:(ncol(vaf_data2)-2)] < 0.4 , na.rm = TRUE)
vaf_data2$`0.4-0.6` <- rowSums(vaf_data2[, 4:(ncol(vaf_data2)-3)] >= 0.4 & vaf_data2[, 4:(ncol(vaf_data2)-3)] < 0.6, na.rm = TRUE)
vaf_data2$`0.6-0.8` <- rowSums(vaf_data2[, 4:(ncol(vaf_data2)-4)] >= 0.6 & vaf_data2[, 4:(ncol(vaf_data2)-4)] < 0.8, na.rm = TRUE)
vaf_data2$`>0.8` <- rowSums(vaf_data2[, 4:(ncol(vaf_data2)-5)] >= 0.8, na.rm = TRUE)

vaf_data1 <- vaf_data1 %>% unite("loc1", loc:alt, sep = ">", remove = T)
colnames(vaf_data1)[1] <- "loc"
vaf_data1$loc <- sub(">", ":", vaf_data1$loc)

vaf_data2 <- vaf_data2 %>% unite("loc1", loc:alt, sep = ">", remove = T)
colnames(vaf_data2)[1] <- "loc"
vaf_data2$loc <- sub(">", ":", vaf_data2$loc)

# get into "tidy form" for ggplot, select only those SNVs that you want to visualize (enter the threshold value obtained from previous histogram)
# Split datasets into 2 parts: 1) where the sum of cells is above threshold and 2) where the sum is below or equal to threshold
# If the bins have changed from <0.2, 0.2-0.4....then replace '6' with (# of bins+1)
vaf_data1 <- vaf_data1 %>% select(1, tail(seq_along(vaf_data1), 6))
vaf_data1 <- vaf_data1 %>% gather("VAF", "counts", -loc, -gene)

vaf_data2 <- vaf_data2 %>% select(1, tail(seq_along(vaf_data2), 6))
vaf_data2 <- vaf_data2 %>% gather("VAF", "counts", -loc, -gene)


# Prepare for labels (if required) and facettting, create stacked ggpot, flip coord, remove unwanted gridlines
# If the bins have changed from <0.2, 0.2-0.4....then replace all three '5's with (# of bins)
vaf_data1$VAF <- factor(vaf_data1$VAF, levels = unique(vaf_data1$VAF))  # if the levels of factoring is weird, reorder using the levels option
vaf_data1$num_label <- c(rep(c(1:(nrow(vaf_data1)/5)), times=5))
#vaf_data$chr <- factor(vaf_data$chr, levels = c(1:22, "X"))

vaf_data2$VAF <- factor(vaf_data2$VAF, levels = unique(vaf_data2$VAF))  # if the levels of factoring is weird, reorder using the levels option
vaf_data2$num_label <- c(rep(c(1:(nrow(vaf_data2)/5)), times=5))

# To reverse the ordering of labels just add reverse=TRUE to position_fill() as an argument
g <- ggplot(vaf_data1, aes(x = num_label, y=counts, fill=VAF, width=1)) + geom_bar(position = position_fill(), stat = "identity") + coord_flip()
g + theme_classic() + ylab("percent cell count") + scale_x_continuous("SNV", breaks = 1:(nrow(vaf_data1)/5), labels = vaf_data1$loc[1:(nrow(vaf_data1)/5)], sec.axis = dup_axis(name = "GENE", labels = vaf_data1$gene[1:(nrow(vaf_data1)/5)])) + 
  theme(axis.title.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold"))


g2 <- ggplot(vaf_data2, aes(x = num_label, y=counts, fill=VAF, width=1)) + geom_bar(position = position_fill(), stat = "identity") + coord_flip()
g2 + theme_classic() + ylab("percent cell count") + scale_x_continuous("SNV", breaks = 1:(nrow(vaf_data2)/5), labels = vaf_data2$loc[1:(nrow(vaf_data2)/5)], sec.axis = dup_axis(name = "GENE", labels = vaf_data2$gene[1:(nrow(vaf_data2)/5)])) + 
  theme(axis.ticks.y = element_blank(), axis.title.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold"))



