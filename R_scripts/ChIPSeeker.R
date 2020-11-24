###
#Wherever possible the set of commands that were run for individual files will, at the end, include a combination of all!
###

library(ChIPseeker)
library(data.table)
library(gtools)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler) # GO and KEGG enrichment analysis
library(ReactomePA) # Reactome Pathway Analysis
library(DOSE)   # Disease Ontonlogy
library(stringr)  # Needed if the dotplots have longer y labels (Descriptions)
library(msigdbr)  # Msig db
library(fgsea)   # For GSEA
library(plyr)
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
require(ggplot2)
library(RColorBrewer)
library(gplots)



txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

setwd("D:/Virtualbox/Ubuntu/Shared Folder/Bukrinsky/ATAC-Seq_First set/All peaks/")

# Get a gene set for narrowing down the analysis (The inflammatory genes)
# gene_set1 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" ) %>% dplyr::filter(gs_id == "M13657") %>% dplyr::select(7, 8)
# gs_for_GO <- as.vector(gene_set[,7])
# colnames(gene_set1)[1] <- "ENTREZID"
# gene_set1$ENTREZID <- as.character(gene_set1$ENTREZID)


file_path <- "D:/Virtualbox/Ubuntu/Shared Folder/Bukrinsky/ATAC-Seq_First set/All peaks/"
ls <- mixedsort(list.files(path = file_path, pattern = "*.narrowPeak"), decreasing = T) # List files in decreasing order so that matrix is formed in increasing order

named_list <- list("2_1" = paste0(file_path, ls[[1]]), "2_2" = paste0(file_path, ls[[2]]), "2_3" = paste0(file_path, ls[[3]]), 
                   "2_4" = paste0(file_path, ls[[4]]), "1_1" = paste0(file_path, ls[[5]]), "1_2" = paste0(file_path, ls[[6]])
                   , "1_3" = paste0(file_path, ls[[7]]), "1_4" = paste0(file_path, ls[[8]])) # If multiple file needs to be included please create a named list of all those files


# For each set of file
peak1 <- readPeakFile(named_list[[1]], header=F)   
peak2 <- readPeakFile(named_list[[2]], header=F)
peak3 <- readPeakFile(named_list[[3]], header=F)
peak4 <- readPeakFile(named_list[[4]], header=F)
peak5 <- readPeakFile(named_list[[5]], header=F)   
peak6 <- readPeakFile(named_list[[6]], header=F)
peak7 <- readPeakFile(named_list[[7]], header=F)
peak8 <- readPeakFile(named_list[[8]], header=F)



# Create dummy GRanges object with only necessary chromosomes
gr <- GRanges(seqnames = c(paste0('chr', 1:22), 'chrX', 'chrY'), ranges = IRanges(start = 1, end = 50)) 

# Filter out unwanted chromosomes from each file using the dummy GRanges object and relevel the factors
peak1_f <- peak1[-(union(grep(".*_+", seqnames(peak1)), grep("^chrM", seqnames(peak1))))@values, ]
seqlevels(peak1_f) = as.character((seqnames(gr)))
peak2_f <- peak2[-(union(grep(".*_+", seqnames(peak2)), grep("^chrM", seqnames(peak2))))@values, ]
seqlevels(peak2_f) = as.character((seqnames(gr)))
peak3_f <- peak3[-(union(grep(".*_+", seqnames(peak3)), grep("^chrM", seqnames(peak3))))@values, ]
seqlevels(peak3_f) = as.character((seqnames(gr)))
peak4_f <- peak4[-(union(grep(".*_+", seqnames(peak4)), grep("^chrM", seqnames(peak4))))@values, ]
seqlevels(peak4_f) = as.character((seqnames(gr)))
peak5_f <- peak5[-(union(grep(".*_+", seqnames(peak5)), grep("^chrM", seqnames(peak5))))@values, ]
seqlevels(peak5_f) = as.character((seqnames(gr)))
peak6_f <- peak6[-(union(grep(".*_+", seqnames(peak6)), grep("^chrM", seqnames(peak6))))@values, ]
seqlevels(peak6_f) = as.character((seqnames(gr)))
peak7_f <- peak7[-(union(grep(".*_+", seqnames(peak7)), grep("^chrM", seqnames(peak7))))@values, ]
seqlevels(peak7_f) = as.character((seqnames(gr)))
peak8_f <- peak8[-(union(grep(".*_+", seqnames(peak8)), grep("^chrM", seqnames(peak8))))@values, ]
seqlevels(peak8_f) = as.character((seqnames(gr)))


# Combine all files into one GRangelist
peak <- GRangesList("2_1" = peak1_f, "2_2" = peak2_f, "2_3" = peak3_f, 
                    "2_4" = peak4_f, "1_1" = peak5_f, "1_2" = peak3_f, 
                    "1_3" = peak4_f, "1_4" = peak5_f)

# remove the unnecessary files
remove(peak1, peak2, peak3, peak4, peak5, peak6, peak7, peak8)

# Coverage plots
covplot(peak1_f, weightCol = 5)
covplot(peak2_f, weightCol = 5)
covplot(peak3_f, weightCol = 5)
covplot(peak4_f, weightCol = 5)
covplot(peak5_f, weightCol = 5)
covplot(peak6_f, weightCol = 5)
covplot(peak7_f, weightCol = 5)
covplot(peak8_f, weightCol = 5)

# covplot(peak, weightCol = 5) # Overall

##### prepare the TSS regions, which are defined as the flanking sequence of the TSS sites #####

## WE CAN MAKE CUSTOM SITES FOR VISUALIZATION

# Create a tagmatrix that contains those peaks that align to the TSS flanking regions
# We can provide the 'windows' tag with the 
# promoter1 <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# tagMatrix1 <- getTagMatrix(peak1_f, windows=promoter1)
# tagMatrix2 <- getTagMatrix(peak2_f, windows=promoter1)
# tagMatrix3 <- getTagMatrix(peak3_f, windows=promoter1)
# tagMatrix4 <- getTagMatrix(peak4_f, windows=promoter1)
# tagMatrix5 <- getTagMatrix(peak5_f, windows=promoter1)


# tagMatrixList <- lapply(named_list, getTagMatrix, windows=promoter1)

# Plot Heatmap of ATACSeq binding to TSS regions from TxDb data
# tagHeatmap(tagMatrix1, xlim = c(-3000, 3000), color="red")
# tagHeatmap(tagMatrix2, xlim = c(-3000, 3000), color="red")
# tagHeatmap(tagMatrix3, xlim = c(-3000, 3000), color="red")
# tagHeatmap(tagMatrix4, xlim = c(-3000, 3000), color="red")
# tagHeatmap(tagMatrix5, xlim = c(-3000, 3000), color="red")


# tagHeatmap(tagMatrixList, xlim = c(-3000, 3000), color="red")

##For the Profile of ChIP peaks binding to start site of Exon/Intron use getBioRegion##########

# Annotate peaks for pairs
peakAnno1 <- annotatePeak(peak1_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno2 <- annotatePeak(peak2_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno3 <- annotatePeak(peak3_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno4 <- annotatePeak(peak4_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno5 <- annotatePeak(peak5_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno6 <- annotatePeak(peak6_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno7 <- annotatePeak(peak7_f, tssRegion=c(-3000, 3000),
TxDb=txdb)
peakAnno8 <- annotatePeak(peak8_f, tssRegion=c(-3000, 3000),
TxDb=txdb)


# Annotate peaks for all pairs
peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

# Find overlap of peaks and annotated genes
genes_o= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes_o)

#Get intersecting genes between all samples and all except the clinical sample
common_1_2 <- intersect(genes_o[[1]], genes_o[[2]])
common_3_4 <- intersect(genes_o[[3]], genes_o[[5]])
common_5_6 <- intersect(genes_o[[6]], genes_o[[7]])
common_5 <- intersect(common_1_2, common_3_4)
common_6 <- intersect(common_5, common_5_6)
common <- intersect(common_5_6, genes_o[[8]])
# q1_q4 <- intersect(genes_o[[4]], genes_o[[5]])

common_df <- bitr(common, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = "org.Hs.eg.db")
# fwrite(common_df, "Common genes_No_AnH2_4.txt", sep = "\t")


# make df from anno file and select the columns: 
vis_1 <- as.data.frame(peakAnno1) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_1 <- vis_1[grepl("Promoter*", vis_1$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_1) <- 1:nrow(vis_1)
colnames(vis_1) <- c("2_1", "genes")
vis_1$genes <- as.numeric(vis_1$genes)
vis_1 <- ddply(vis_1,"genes",numcolwise(sum))
vis_2 <- as.data.frame(peakAnno2) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_2 <- vis_2[grepl("Promoter*", vis_2$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_2) <- 1:nrow(vis_2)
colnames(vis_2) <- c("2_2", "genes")
vis_2$genes <- as.numeric(vis_2$genes)
vis_2 <- ddply(vis_2,"genes",numcolwise(sum))
vis_3 <- as.data.frame(peakAnno3) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_3 <- vis_3[grepl("Promoter*", vis_3$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_3) <- 1:nrow(vis_3)
colnames(vis_3) <- c("2_3", "genes")
vis_3$genes <- as.numeric(vis_3$genes)
vis_3 <- ddply(vis_3,"genes",numcolwise(sum))
vis_4 <- as.data.frame(peakAnno4) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_4 <- vis_4[grepl("Promoter*", vis_4$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_4) <- 1:nrow(vis_4)
colnames(vis_4) <- c("2_4","genes")
vis_4$genes <- as.numeric(vis_4$genes)
vis_4 <- ddply(vis_4,"genes",numcolwise(sum))
vis_5 <- as.data.frame(peakAnno5) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_5 <- vis_5[grepl("Promoter*", vis_5$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_5) <- 1:nrow(vis_5)
colnames(vis_5) <- c("1_1", "genes")
vis_5$genes <- as.numeric(vis_5$genes)
vis_5 <- ddply(vis_5,"genes",numcolwise(sum))
vis_6 <- as.data.frame(peakAnno6) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_6 <- vis_6[grepl("Promoter*", vis_6$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_6) <- 1:nrow(vis_6)
colnames(vis_6) <- c("1_2", "genes")
vis_6$genes <- as.numeric(vis_6$genes)
vis_6 <- ddply(vis_6,"genes",numcolwise(sum))
vis_7 <- as.data.frame(peakAnno7) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_7 <- vis_7[grepl("Promoter*", vis_7$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_7) <- 1:nrow(vis_7)
colnames(vis_7) <- c("1_3","genes")
vis_7$genes <- as.numeric(vis_7$genes)
vis_7 <- ddply(vis_7,"genes",numcolwise(sum))
vis_8 <- as.data.frame(peakAnno8) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_8 <- vis_8[grepl("Promoter*", vis_8$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_8) <- 1:nrow(vis_8)
colnames(vis_8) <- c("1_4", "genes")
vis_8$genes <- as.numeric(vis_8$genes)
vis_8 <- ddply(vis_8,"genes",numcolwise(sum))


# Create a smaller list by getting only those genes that were present in both RNASeq and ATACSeq
vis <- list(vis_1, vis_2, vis_3, vis_4, vis_5, vis_6, vis_7, vis_8) %>% purrr::reduce(full_join, by = "genes") %>% distinct_at(1, .keep_all = T)
colnames(vis)[1] <- "ENTREZID"
vis$ENTREZID <- as.character(vis$ENTREZID)
vis <- inner_join(vis, common_df, by = "ENTREZID")
vis <- vis[, c(1, ncol(vis)-1, ncol(vis), 2:(ncol(vis)-2))] %>% distinct()

fwrite(vis, "Common genes_All_AnH.txt", sep = "\t")


# vis_bk = vis
# 
# vis_q1_4 <- list(vis_4, vis_5) %>% purrr::reduce(full_join, by = "genes") %>% distinct_at(1, .keep_all = T)
# colnames(vis_q1_4)[1] <- "ENTREZID"
# vis_q1_4$ENTREZID <- as.character(vis_q1_4$ENTREZID)
# vis_q1_4 <- inner_join(vis_q1_4, common_df, by = "ENTREZID")
# vis_q1_4 <- vis_q1_4[, c(1, 4, 5, 2, 3)] %>% distinct()
# 
# 
# # intersect with inflammatory genes from msigdb
# vis <- inner_join(vis, gene_set1, by = "ENTREZID")
# vis <- vis %>% dplyr::select(-9)
# 
# # Prepare for clustering
# mydata1 <- as.matrix(vis[,c(4:8)])
# 
# # Rownames as gene symbols   
# rownames(mydata1) <- vis[, 3]
# 
# 
# # Find the range for the numeric value to plot and then divide it into parts for each color
# # Plot Heatmaps for the common genes
# split1 <- colorRampPalette(c("red"))(n=49)
# split2 <- colorRampPalette(c("green", "blue"), bias = 1.5)(n=250)
# my_palette <- c(rev(split2), split1)
# colors1 = unique(c(seq(14, 100, length=100), seq(100, 400, length=100), seq(400, 1076, length=102)))
# colors2 = unique(c(seq(14, 100, length=100), seq(100, 450, length=100), seq(450, 644, length=102)))
# 
# 
# # Prepare for clustering
# mydata3 <- mydata1[rowSums(is.na(mydata1)) == 0, ]
# distance1 = dist(mydata3, method = "euclidean")   # Manhattan distance for clustering
# cluster1 = hclust(distance1, method = "ward.D2")    
# 
# 
# 
# # Plot the Heatmap
# heatmap.2(mydata3,
#           main = "Comparison", # heat map title
#           xlab = "", # X-axis title
#           ylab = "Gene Names", # Y-axis title
#           notecol="black",      # change font color of cell labels to black
#           density.info="density",  # draw a density plot inside the colour key
#           key.title = "Color Key", # title for the colour key
#           key.xlab = "Peaks", # X-axis label for colour key
#           symkey = F, # Colour key is not symmetric about 0
#           symm = F,
#           symbreaks = T,
#           sepwidth = c(0.07, 0.07),
#           srtCol = 45, # Angle in degrees for column labels
#           cexCol = 1.2, # Column label font size
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(5,6),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           breaks = colors1, # define color range in breaks
#           dendrogram="row",     # dendrograms considerng only row
#           Rowv = as.dendrogram(cluster1), # apply default clustering method
#           Colv = FALSE) # No re-ordering for columns
# 
# 
# # Prepare for clustering
# mydata2 <- as.matrix(vis_bk[,c(4:8)])
# 
# # Rownames as gene symbols   
# rownames(mydata2) <- vis_bk[, 3]
# 
# # Prepare for clustering
# mydata4 <- mydata2[rowSums(is.na(mydata2)) == 0, ]
# distance2 = dist(mydata4, method = "euclidean")   # Manhattan distance for clustering
# cluster2 = hclust(distance2, method = "ward.D2")    
# 
# 
# 
# # Plot the Heatmap
# heatmap.2(mydata4,
#           main = "Comparison", # heat map title
#           xlab = "", # X-axis title
#           ylab = "Gene Names", # Y-axis title
#           notecol="black",      # change font color of cell labels to black
#           density.info="density",  # draw a density plot inside the colour key
#           key.title = "Color Key", # title for the colour key
#           key.xlab = "Peaks", # X-axis label for colour key
#           symkey = F, # Colour key is not symmetric about 0
#           symm = F,
#           symbreaks = T,
#           sepwidth = c(0.07, 0.07),
#           srtCol = 45, # Angle in degrees for column labels
#           cexCol = 1.2, # Column label font size
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(5,6),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           breaks = colors1, # define color range in breaks
#           dendrogram="row",     # dendrograms considerng only row
#           Rowv = as.dendrogram(cluster2), # apply default clustering method
#           Colv = FALSE) # No re-ordering for columns
# 
# # Prepare for clustering
# mydata5 <- as.matrix(vis_q1_4[,c(4:5)])
# 
# # Rownames as gene symbols   
# rownames(mydata5) <- vis_q1_4[, 3]
# 
# # Prepare for clustering
# mydata6 <- mydata5[rowSums(is.na(mydata5)) == 0, ]
# distance3 = dist(mydata6, method = "euclidean")   # Manhattan distance for clustering
# cluster3 = hclust(distance3, method = "ward.D2")    
# 
# 
# 
# # Plot the Heatmap
# heatmap.2(mydata6,
#           main = "Comparison", # heat map title
#           xlab = "", # X-axis title
#           ylab = "Gene Names", # Y-axis title
#           notecol="black",      # change font color of cell labels to black
#           density.info="density",  # draw a density plot inside the colour key
#           key.title = "Color Key", # title for the colour key
#           key.xlab = "Peaks", # X-axis label for colour key
#           symkey = F, # Colour key is not symmetric about 0
#           symm = F,
#           symbreaks = T,
#           sepwidth = c(0.07, 0.07),
#           srtCol = 45, # Angle in degrees for column labels
#           cexCol = 1.2, # Column label font size
#           trace="none",         # turns off trace lines inside the heat map
#           margins =c(5,6),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier
#           breaks = colors1, # define color range in breaks
#           dendrogram="row",     # dendrograms considerng only row
#           Rowv = as.dendrogram(cluster3), # apply default clustering method
#           Colv = FALSE) # No re-ordering for columns

