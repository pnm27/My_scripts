library(data.table)
library(tidyverse)
library(gridExtra)
library(grid)

############################################################
#Input format for matrix should be: colnames as "SNV, sample1, sample2.....sampleN (doesn't depend on the naming style except at the simplification step
# i.e gather function)"
#with SNV having a format like 1:10000_A>T
############################################################


setwd('D:/Virtualbox/Ubuntu/Shared Folder/ReadCounts_subtask/Single_Cell/')



#Input VAF file, split snv, deselect ref and alt and simplify colnames
vaf_data <- fread('SRR11680223_10.vaf.matrix.tsv') %>% separate(SNV, into = c("chr", "loc", "ref", "alt"))
vaf_data$loc <- as.integer(vaf_data$loc)



# Annotate the VAFs with customised gene list, order by chr and loc, create a new "formatted" location column, and rename it, keep chr column for grouping
setorder(vaf_data, chr, loc)
vaf_data <- vaf_data %>% unite("loc1", chr:loc, remove = T)
vaf_data$loc1 <- sapply(vaf_data$loc1, function (x) paste0("chr", x))
names(vaf_data)[1] <- "loc"


#Tidy the data: count for each bins, group by snvs
# If the bins have changed from <0.2, 0.2-0.4....then replace '4' with (# of bins-1)
vaf_data <- vaf_data %>% group_by(loc)
vaf_data$`<0.2` <- rowSums(vaf_data[, 4:(ncol(vaf_data)-1)] < 0.2, na.rm = TRUE)
vaf_data$`0.2-0.4` <- rowSums(vaf_data[, 4:(ncol(vaf_data)-2)] >= 0.2 & vaf_data[, 4:(ncol(vaf_data)-2)] < 0.4 , na.rm = TRUE)
vaf_data$`0.4-0.6` <- rowSums(vaf_data[, 4:(ncol(vaf_data)-3)] >= 0.4 & vaf_data[, 4:(ncol(vaf_data)-3)] < 0.6, na.rm = TRUE)
vaf_data$`0.6-0.8` <- rowSums(vaf_data[, 4:(ncol(vaf_data)-4)] >= 0.6 & vaf_data[, 4:(ncol(vaf_data)-4)] < 0.8, na.rm = TRUE)
vaf_data$`>0.8` <- rowSums(vaf_data[, 4:(ncol(vaf_data)-5)] >= 0.8, na.rm = TRUE)

# For sort-order
vaf_data$`perc_0` <- vaf_data$`<0.2`/rowSums(vaf_data[, (ncol(vaf_data)-4):ncol(vaf_data)], na.rm = TRUE)
vaf_data$`perc_0.2` <- vaf_data$`0.2-0.4`/rowSums(vaf_data[, (ncol(vaf_data)-5):(ncol(vaf_data)-1)], na.rm = TRUE)
vaf_data$`perc_0.4` <- vaf_data$`0.4-0.6`/rowSums(vaf_data[, (ncol(vaf_data)-6):(ncol(vaf_data)-2)], na.rm = TRUE)
vaf_data$`perc_0.6` <- vaf_data$`0.6-0.8`/rowSums(vaf_data[, (ncol(vaf_data)-7):(ncol(vaf_data)-3)], na.rm = TRUE)
vaf_data$`perc_0.8` <- vaf_data$`>0.8`/rowSums(vaf_data[, (ncol(vaf_data)-8):(ncol(vaf_data)-4)], na.rm = TRUE)
# For removing SNVs with 0 exp (including NAs)
vaf_data$zero_exp <- rowSums(vaf_data[, 4:(ncol(vaf_data)-10)] > 0, na.rm = TRUE)
vaf_data <- vaf_data %>% filter(zero_exp != 0)
vaf_data <- vaf_data %>% select(-zero_exp)

vaf_data <- vaf_data %>% unite("loc1", loc:alt, sep = ">", remove = T)
colnames(vaf_data)[1] <- "loc"
vaf_data$loc <- sub(">", ":", vaf_data$loc)



# get into "tidy form" for ggplot, select only those SNVs that you want to visualize (enter the threshold value obtained from previous histogram)
# Split datasets into 2 parts: 1) where the sum of cells is above threshold and 2) where the sum is below or equal to threshold
# If the bins have changed from <0.2, 0.2-0.4....then replace '5' with (# of bins)
vaf_data <- vaf_data %>% select(1, tail(seq_along(vaf_data), 10))
vaf_data <- vaf_data %>% gather("VAF", "counts", -loc, -`perc_0.8`, -`perc_0.6`, -`perc_0.4`, -`perc_0.2`, -`perc_0`)

# Prepare for labels (if required) and facettting, create stacked ggpot, flip coord, remove unwanted gridlines
# If the bins have changed from <0.2, 0.2-0.4....then replace all three '5's with (# of bins)
vaf_data$VAF <- factor(vaf_data$VAF, levels = unique(vaf_data$VAF))  # if the levels of factoring is weird, reorder using the levels option
setorder(vaf_data, VAF, `perc_0.8`, `perc_0.6`, `perc_0.4`, `perc_0.2`, `perc_0`)
# setorder(vaf_data, VAF, `perc_0`, `perc_0.2`, `perc_0.4`, `perc_0.6`, `perc_0.8`)
# setorder(vaf_data, VAF, `perc_0.2`)
vaf_data$num_label <- c(rep(c(1:(nrow(vaf_data)/5)), times=5))
#vaf_data$chr <- factor(vaf_data$chr, levels = c(1:22, "X"))


vaf_data <- vaf_data %>% select(-`perc_0.8`, -`perc_0.6`, -`perc_0.4`, -`perc_0.2`, -`perc_0`)


# To reverse the ordering of labels just add reverse=TRUE to position_fill() as an argument
# png("new.png", width = 800, height = 773)
g <- ggplot(vaf_data, aes(x = num_label, y=counts, fill=VAF, width=1)) + geom_bar(position = position_fill(), stat = "identity") + coord_flip()

# Here is the make-up artist for the graph
g + theme_classic() + ylab("") + scale_x_continuous("", breaks = 1:(nrow(vaf_data)/5), labels = vaf_data$loc[1:(nrow(vaf_data)/5)]) + 
  theme(axis.text.x =  element_text(size = 16, face = "bold"), axis.ticks.x = element_line(size = 2), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 14, face = "bold")
        , legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 12, face = "bold")) # Controls legend
# dev.off()


