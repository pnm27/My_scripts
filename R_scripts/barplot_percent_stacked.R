library(data.table)
library(tidyverse)
library(gridExtra)
library(grid)

############################################################
#Input format for matrix should be: colnames as "SNV, sample1, sample2.....sampleN (doesn't depend on the naming style except at the simplification step
# i.e gather function)"
#with SNV having a format like 1:10000_A>T
############################################################


setwd('../Single_cell/Bar plots/')

#Input VAF file, split snv, deselect ref and alt and simplify colnames
vaf_data <- fread('N8_3_VAF_matrix_CLEAN_10c_nExtrRmean.txt') %>% separate(SNV, into = c("chr", "loc", "ref", "alt")) %>% select(-ref, -alt)
vaf_data$loc <- as.integer(vaf_data$loc)
colnames(vaf_data) <- sapply(colnames(vaf_data), function (x) gsub("_wasp_AlignedByCoord_vW_filt_dedupped_sorted", "", x))


# Annotate the VAFs with customised gene list, order by chr and loc, create a new "formatted" location column, and rename it, keep chr column for grouping
#vaf_data <- left_join(vaf_data, gene_ref, by = "chr") %>% filter(loc >= start & loc <= end) %>% distinct() %>% select(-start, -end)
setorder(vaf_data, chr, loc)
vaf_data <- vaf_data %>% unite("loc1", chr:loc, remove = F) %>% select(-loc)
vaf_data$loc1 <- sapply(vaf_data$loc1, function (x) paste0("chr", x))
names(vaf_data)[1] <- "loc"


#Tidy the data: count for each bins, group by snvs
# If the bins have changed from <0.2, 0.2-0.4....then replace '4' with (# of bins-1)
vaf_data <- vaf_data %>% group_by(loc)
vaf_data$`<0.2` <- rowSums(vaf_data < 0.2, na.rm = TRUE)
vaf_data$`0.2-0.4` <- rowSums(vaf_data >= 0.2 & vaf_data < 0.4 , na.rm = TRUE)
vaf_data$`0.4-0.6` <- rowSums(vaf_data >= 0.4 & vaf_data < 0.6, na.rm = TRUE)
vaf_data$`0.6-0.8` <- rowSums(vaf_data >= 0.6 & vaf_data < 0.8, na.rm = TRUE)
vaf_data$`>0.8` <- rowSums(vaf_data[, 3:(ncol(vaf_data)-4)] >= 0.8, na.rm = TRUE)


# get into "tidy form" for ggplot, select only those SNVs that you want to visualize (enter the threshold value obtained from previous histogram)
# Split datasets into 2 parts: 1) where the sum of cells is above threshold and 2) where the sum is below or equal to threshold
# If the bins have changed from <0.2, 0.2-0.4....then replace '6' with (# of bins+1)
vaf_data <- vaf_data %>% select(1, 2, tail(seq_along(vaf_data), 5))
vaf_data <- vaf_data %>% gather("VAF", "counts", -loc, -chr)


# Prepare for labels (if required) and facettting, create stacked ggpot, flip coord, remove unwanted gridlines
# If the bins have changed from <0.2, 0.2-0.4....then replace all three '5's with (# of bins)
vaf_data$VAF <- factor(vaf_data$VAF, levels = unique(vaf_data$VAF))  # if the levels of factoring is weird, reorder using the levels option
vaf_data$num_label <- c(rep(c(1:(nrow(vaf_data)/5)), times=5))
vaf_data$chr <- factor(vaf_data$chr, levels = c(1:22, "X"))

# To reverse the ordering of labels just add reverse=TRUE to position_fill() as an argument
g <- ggplot(vaf_data, aes(x = num_label, y=counts, fill=VAF, width=1)) + geom_bar(position = position_fill(), stat = "identity") + coord_flip()
g + theme_classic() + ylab("percent cell count") + scale_x_continuous("", breaks = 1:(nrow(vaf_data)/5), labels = vaf_data$loc[1:(nrow(vaf_data)/5)]) + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold")) + facet_wrap(~chr, scales = "free", ncol = 8)




