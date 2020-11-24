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

setwd("../Bukrinsky/ATAC-Seq_Feb_28/macs_new/")

# Get a gene set for narrowing down the analysis (The inflammatory genes)
gene_set1 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" ) %>% dplyr::filter(gs_id == "M13657") %>% dplyr::select(7, 8)
# gs_for_GO <- as.vector(gene_set[,7])
colnames(gene_set1)[1] <- "ENTREZID"
gene_set1$ENTREZID <- as.character(gene_set1$ENTREZID)


file_path <- "D:/Virtualbox/Ubuntu/Shared Folder/Bukrinsky/ATAC-Seq_Feb_28/macs_new/"
ls <- mixedsort(list.files(path = file_path, pattern = "*.narrowPeak"), decreasing = T) # List files in decreasing order so that matrix is formed in increasing order

named_list <- list("U-7_U-8" = paste0(file_path, ls[[1]]), "U-9_U-10" = paste0(file_path, ls[[2]]), "U-11_U-12" = paste0(file_path, ls[[3]]), "S_all" = paste0(file_path, ls[[4]]), "S-13_S-14" = paste0(file_path, ls[[4]]),
                   "S-15_S-16" = paste0(file_path, ls[[6]]), "Q-1_Q-2" = paste0(file_path, ls[[7]]), "Q-3_Q-4" = paste0(file_path, ls[[8]]), "Q-5_Q-6" = paste0(file_path, ls[[9]])) # If multiple file needs to be included please create a named list of all those files


# For each set of file
peak1 <- readPeakFile(named_list[[1]], header=F)   
peak2 <- readPeakFile(named_list[[2]], header=F)
peak3 <- readPeakFile(named_list[[3]], header=F)
peak4 <- readPeakFile(named_list[[4]], header=F)
peak5 <- readPeakFile(named_list[[5]], header=F)   
peak6 <- readPeakFile(named_list[[6]], header=F)
peak7 <- readPeakFile(named_list[[7]], header=F)
peak8 <- readPeakFile(named_list[[8]], header=F)
peak9 <- readPeakFile(named_list[[9]], header=F)


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
peak9_f <- peak9[-(union(grep(".*_+", seqnames(peak9)), grep("^chrM", seqnames(peak8))))@values, ]
seqlevels(peak9_f, pruning.mode="coarse") = as.character((seqnames(gr)))

# Combine all files into one GRangelist
peak <- GRangesList("S-13_S-14" = peak5_f, "S-15_S-16" = peak6_f, "Q-1_Q-2" = peak7_f, "Q-3_Q-4" = peak8_f, "U-9_U-10" = peak2_f)

# remove the unnecessary files
remove(peak1, peak2, peak3, peak4, peak5, peak6, peak7, peak8, peak9)

# Coverage plots
covplot(peak1_f, weightCol = 5)
covplot(peak2_f, weightCol = 5)
covplot(peak3_f, weightCol = 5)
covplot(peak4_f, weightCol = 5)
covplot(peak5_f, weightCol = 5)
covplot(peak6_f, weightCol = 5)
covplot(peak7_f, weightCol = 5)
covplot(peak8_f, weightCol = 5)
covplot(peak, weightCol = 5) # Overall

##### prepare the TSS regions, which are defined as the flanking sequence of the TSS sites #####

## WE CAN MAKE CUSTOM SITES FOR VISUALIZATION

# Create a tagmatrix that contains those peaks that align to the TSS flanking regions
# We can provide the 'windows' tag with the 
promoter1 <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix1 <- getTagMatrix(peak1_f, windows=promoter1)
tagMatrix2 <- getTagMatrix(peak2_f, windows=promoter1)
tagMatrix3 <- getTagMatrix(peak3_f, windows=promoter1)
tagMatrix4 <- getTagMatrix(peak4_f, windows=promoter1)
tagMatrix5 <- getTagMatrix(peak5_f, windows=promoter1)
tagMatrix6 <- getTagMatrix(peak6_f, windows=promoter1)
tagMatrix7 <- getTagMatrix(peak7_f, windows=promoter1)
tagMatrix8 <- getTagMatrix(peak8_f, windows=promoter1)
tagMatrix9 <- getTagMatrix(peak9_f, windows=promoter1)

# tagMatrixList <- lapply(named_list, getTagMatrix, windows=promoter1)

# Plot Heatmap of ATACSeq binding to TSS regions from TxDb data
tagHeatmap(tagMatrix1, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix2, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix3, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix4, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix5, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix6, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix7, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix8, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix9, xlim = c(-3000, 3000), color="red")

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
peakAnno9 <- annotatePeak(peak9_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)

# Annotate peaks for all pairs
peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

# Find overlap of peaks and annotated genes
genes_o= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes_o)

#Get intersecting genes between all samples and all except the clinical sample
common_1_2 <- intersect(genes_o[[1]], genes_o[[2]])
common_3_4 <- intersect(genes_o[[3]], genes_o[[4]])
common_5 <- intersect(common_1_2, common_3_4)
common <- intersect(common_5, genes_o[[5]])

common_df <- bitr(common, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = "org.Hs.eg.db")
fwrite(common_df, "Common genes_Q1_4_S_all_U_9.txt", sep = "\t")


# make df from anno file and select the columns: 
vis_1 <- as.data.frame(peakAnno5) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_1 <- vis_1[grepl("Promoter*", vis_1$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_1) <- 1:nrow(vis_1)
colnames(vis_1) <- c("S13_S14", "genes")
vis_1$genes <- as.numeric(vis_1$genes)
vis_1 <- ddply(vis_1,"genes",numcolwise(sum))
vis_2 <- as.data.frame(peakAnno6) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_2 <- vis_2[grepl("Promoter*", vis_2$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_2) <- 1:nrow(vis_2)
colnames(vis_2) <- c("S15_S16", "genes")
vis_2$genes <- as.numeric(vis_2$genes)
vis_2 <- ddply(vis_2,"genes",numcolwise(sum))
vis_3 <- as.data.frame(peakAnno7) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_3 <- vis_3[grepl("Promoter*", vis_3$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_3) <- 1:nrow(vis_3)
colnames(vis_3) <- c("Q1_Q2","genes")
vis_3$genes <- as.numeric(vis_3$genes)
vis_3 <- ddply(vis_3,"genes",numcolwise(sum))
vis_4 <- as.data.frame(peakAnno8) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_4 <- vis_4[grepl("Promoter*", vis_4$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_4) <- 1:nrow(vis_4)
colnames(vis_4) <- c("Q3_Q4", "genes")
vis_4$genes <- as.numeric(vis_4$genes)
vis_4 <- ddply(vis_4,"genes",numcolwise(sum))
vis_5 <- as.data.frame(peakAnno2) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_5 <- vis_5[grepl("Promoter*", vis_5$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_5) <- 1:nrow(vis_5)
colnames(vis_5) <- c("U9_U10", "genes")
vis_5$genes <- as.numeric(vis_5$genes)
vis_5 <- ddply(vis_5,"genes",numcolwise(sum))

# Create a smaller list by getting only those genes that were present in both RNASeq and ATACSeq
vis <- list(vis_1, vis_2, vis_3, vis_4, vis_5) %>% purrr::reduce(full_join, by = "genes") %>% distinct_at(1, .keep_all = T)
colnames(vis)[1] <- "ENTREZID"
vis$ENTREZID <- as.character(vis$ENTREZID)
vis <- inner_join(vis, common_df, by = "ENTREZID")
vis <- vis[, c(1, 7, 8, 2, 3, 4, 5, 6)] %>% distinct()

# Simplify colnames for vis and vis_3pairs
# colnames(vis)[3:5] <- c("MFa2-2", "MFa1-2", "MFa1-4")

# intersect with inflammatory genes from msigdb
vis <- inner_join(vis, gene_set1, by = "ENTREZID")
vis <- vis %>% dplyr::select(-9)

# Prepare for no clustering
mydata1 <- as.matrix(vis[,c(4:8)])

# Rownames as gene symbols   
rownames(mydata1) <- vis[, 3]


# Find the range for the numeric value to plot and then divide it into parts for each color
# Plot Heatmaps for the common genes
split1 <- colorRampPalette(c("red"))(n=49)
split2 <- colorRampPalette(c("green", "blue"), bias = 1.5)(n=250)
my_palette <- c(rev(split2), split1)
colors1 = unique(c(seq(14, 100, length=100), seq(100, 400, length=100), seq(400, 1076, length=102)))
colors2 = unique(c(seq(14, 100, length=100), seq(100, 450, length=100), seq(450, 644, length=102)))


# Prepare for clustering
mydata3 <- mydata1[rowSums(is.na(mydata1)) == 0, ]
distance1 = dist(mydata3, method = "euclidean")   # Manhattan distance for clustering
cluster1 = hclust(distance1, method = "ward.D2")    



# Plot the Heatmap
heatmap.2(mydata3,
          main = "Comparison", # heat map title
          xlab = "", # X-axis title
          ylab = "Gene Names", # Y-axis title
          notecol="black",      # change font color of cell labels to black
          density.info="density",  # draw a density plot inside the colour key
          key.title = "Color Key", # title for the colour key
          key.xlab = "Peaks", # X-axis label for colour key
          symkey = F, # Colour key is not symmetric about 0
          symm = F,
          symbreaks = T,
          sepwidth = c(0.07, 0.07),
          srtCol = 45, # Angle in degrees for column labels
          cexCol = 1.2, # Column label font size
          trace="none",         # turns off trace lines inside the heat map
          margins =c(5,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks = colors1, # define color range in breaks
          dendrogram="row",     # dendrograms considerng only row
          Rowv = as.dendrogram(cluster1), # apply default clustering method
          Colv = FALSE) # No re-ordering for columns


