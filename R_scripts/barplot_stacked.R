library(data.table)
library(tidyverse)
library(gridExtra)
library(grid)

############################################################
#Input format for matrix should be: colnames as "SNV, sample1, sample2.....sampleN (doesn't depend on the naming style except at the simplification step
# i.e gather function)"
#with SNV having a format like 1:10000_A>T
############################################################

# Set up working directory (VAF file is expected to be here)
setwd('../Single_cell/Bar plots/')

#Input VAF file, split snv, deselect ref and alt and simplify colnames, if needed
vaf_data <- fread('N8_3_VAF_matrix_CLEAN_10c_nExtrRmean.txt') %>% separate(SNV, into = c("chr", "loc", "ref", "alt")) %>% select(-ref, -alt)
vaf_data$loc <- as.integer(vaf_data$loc)
colnames(vaf_data) <- sapply(colnames(vaf_data), function (x) gsub("_wasp_AlignedByCoord_vW_filt_dedupped_sorted", "", x))


# Order by chr and loc, create a new "formatted" location column, and rename it while removing the old columns (chr and loc)
#vaf_data <- left_join(vaf_data, gene_ref, by = "chr") %>% filter(loc >= start & loc <= end) %>% distinct() %>% select(-start, -end)
setorder(vaf_data, chr, loc) # If you want to reproduce the graph in CSHL poster don't run this line
vaf_data <- vaf_data %>% unite("loc1", chr:loc, remove = T)
#vaf_data <- vaf_data[, c(1, ncol(vaf_data), 3:(ncol(vaf_data)-1))]
vaf_data$loc1 <- sapply(vaf_data$loc1, function (x) paste0("chr", x))
names(vaf_data)[1] <- "loc"


#Tidy the data: count for each bins, group by snvs
# If the bins have changed from <0.2, 0.2-0.4....then replace '4' with (# of bins-1)
vaf_data <- vaf_data %>% group_by(loc)
vaf_data$`<0.2` <- rowSums(vaf_data < 0.2, na.rm = TRUE)
vaf_data$`0.2-0.4` <- rowSums(vaf_data >= 0.2 & vaf_data < 0.4 , na.rm = TRUE)
vaf_data$`0.4-0.6` <- rowSums(vaf_data >= 0.4 & vaf_data < 0.6, na.rm = TRUE)
vaf_data$`0.6-0.8` <- rowSums(vaf_data >= 0.6 & vaf_data < 0.8, na.rm = TRUE)
vaf_data$`>0.8` <- rowSums(vaf_data[, 2:(ncol(vaf_data)-4)] >= 0.8, na.rm = TRUE)


# Plot a histogram to select a threshold that separates the most informative set of genes from others
# If the bins have changed from <0.2, 0.2-0.4....then replace '4' with (# of bins-1)
vaf_data$sum <- rowSums(vaf_data[,(ncol(vaf_data)-4):ncol(vaf_data)], na.rm = TRUE)
#hist(vaf_data$sum)


# Set the threshold
thresh = 1000

# get into "tidy form" for ggplot, select only those SNVs that you want to visualize (enter the threshold value obtained from previous histogram)
# Split datasets into 2 parts: 1) where the sum of cells is above threshold and 2) where the sum is below or equal to threshold
# If the bins have changed from <0.2, 0.2-0.4....then replace '6' with (# of bins+1)
vaf_data1 <- vaf_data[which(vaf_data$sum > thresh),]
vaf_data1 <- vaf_data1 %>% select(1, tail(seq_along(vaf_data), 6), -ncol(vaf_data))
vaf_data1 <- vaf_data1 %>% gather("VAF", "counts", -loc)

vaf_data2 <- vaf_data[which(vaf_data$sum <= thresh),]
vaf_data2 <- vaf_data2 %>% select(1, tail(seq_along(vaf_data), 6), -ncol(vaf_data))
vaf_data2 <- vaf_data2 %>% gather("VAF", "counts", -loc)

# Prepare for labels (if required), create stacked ggpot, flip coord, remove unwanted gridlines
# If the bins have changed from <0.2, 0.2-0.4....then replace all three '5's with (# of bins)
vaf_data1$VAF <- factor(vaf_data1$VAF, levels = unique(vaf_data1$VAF))  # if the levels of factoring is weird, reorder using the levels option
vaf_data1$num_label <- c(rep(c(1:(nrow(vaf_data1)/5)), times=5))
vaf_data2$VAF <- factor(vaf_data2$VAF, levels = unique(vaf_data2$VAF))  # if the levels of factoring is weird, reorder using the levels option
vaf_data2$num_label <- c(rep(c(1:(nrow(vaf_data2)/5)), times=5))

# To reverse the ordering of labels just add reverse=TRUE to position_stack() as an argument
g1 <- ggplot(vaf_data1, aes(x = num_label, y=counts, fill=VAF)) + geom_bar(position = position_stack(), stat = "identity") + coord_flip()
g1 + theme_classic() + theme(panel.grid.major.x = element_line(colour="grey", size = 0.5), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + xlab("SNP") + scale_x_continuous("",breaks = 1:(nrow(vaf_data1)/5), labels = vaf_data1$loc[1:(nrow(vaf_data1)/5)])

g2 <- ggplot(vaf_data2, aes(x = num_label, y=counts, fill=VAF)) + geom_bar(position = position_stack(), stat = "identity") + coord_flip()
g2 + theme_classic() + theme(panel.grid.major.x = element_line(colour="grey", size = 0.5), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + xlab("SNP") + scale_x_continuous("",breaks = 1:(nrow(vaf_data2)/5), labels = vaf_data2$loc[1:(nrow(vaf_data2)/5)])


