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

setwd("../Bukrinsky/ATAC-Seq_Feb_28/macs_new/indiv/")

# Get a gene set for narrowing down the analysis (The inflammatory genes)
gene_set1 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" ) %>% dplyr::filter(gs_id == "M13657") %>% dplyr::select(7, 8)
# gs_for_GO <- as.vector(gene_set[,7])
colnames(gene_set1)[1] <- "ENTREZID"
gene_set1$ENTREZID <- as.character(gene_set1$ENTREZID)


file_path <- "D:/Virtualbox/Ubuntu/Shared Folder/Bukrinsky/ATAC-Seq_Feb_28/macs_new/indiv/"
ls <- mixedsort(list.files(path = file_path, pattern = "*.narrowPeak"), decreasing = T) # List files in decreasing order so that matrix is formed in increasing order

named_list <- list("U-7" = paste0(file_path, ls[[1]]), "U-8" = paste0(file_path, ls[[2]]), "U-9" = paste0(file_path, ls[[3]]), "U-10" = paste0(file_path, ls[[4]]), "U-11" = paste0(file_path, ls[[4]]),
                   "U-12" = paste0(file_path, ls[[6]]), "S-16" = paste0(file_path, ls[[7]]), "S-14" = paste0(file_path, ls[[8]]), "S-13" = paste0(file_path, ls[[9]]),
                   "S-15" = paste0(file_path, ls[[10]]), "Q-1" = paste0(file_path, ls[[11]]), "Q-2" = paste0(file_path, ls[[12]]), "Q-3" = paste0(file_path, ls[[13]]),
                   "Q-4" = paste0(file_path, ls[[14]]), "Q-5" = paste0(file_path, ls[[15]]), "Q-6" = paste0(file_path, ls[[16]])) # If multiple file needs to be included please create a named list of all those files


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
peak10 <- readPeakFile(named_list[[10]], header=F)
peak11 <- readPeakFile(named_list[[11]], header=F)
peak12 <- readPeakFile(named_list[[12]], header=F)   
peak13 <- readPeakFile(named_list[[13]], header=F)
peak14 <- readPeakFile(named_list[[14]], header=F)
peak15 <- readPeakFile(named_list[[15]], header=F)
peak16 <- readPeakFile(named_list[[16]], header=F)


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
peak9_f <- peak9[-(union(grep(".*_+", seqnames(peak9)), grep("^chrM", seqnames(peak9))))@values, ]
seqlevels(peak9_f, pruning.mode="coarse") = as.character((seqnames(gr)))
peak10_f <- peak10[-(union(grep(".*_+", seqnames(peak10)), grep("^chrM", seqnames(peak10))))@values, ]
seqlevels(peak10_f) = as.character((seqnames(gr)))
peak11_f <- peak11[-(union(grep(".*_+", seqnames(peak11)), grep("^chrM", seqnames(peak11))))@values, ]
seqlevels(peak11_f) = as.character((seqnames(gr)))
peak12_f <- peak12[-(union(grep(".*_+", seqnames(peak12)), grep("^chrM", seqnames(peak12))))@values, ]
seqlevels(peak12_f) = as.character((seqnames(gr)))
peak13_f <- peak13[-(union(grep(".*_+", seqnames(peak13)), grep("^chrM", seqnames(peak13))))@values, ]
seqlevels(peak13_f) = as.character((seqnames(gr)))
peak14_f <- peak14[-(union(grep(".*_+", seqnames(peak14)), grep("^chrM", seqnames(peak14))))@values, ]
seqlevels(peak14_f) = as.character((seqnames(gr)))
peak15_f <- peak15[-(union(grep(".*_+", seqnames(peak15)), grep("^chrM", seqnames(peak15))))@values, ]
seqlevels(peak15_f) = as.character((seqnames(gr)))
peak16_f <- peak16[-(union(grep(".*_+", seqnames(peak16)), grep("^chrM", seqnames(peak16))))@values, ]
seqlevels(peak16_f, pruning.mode="coarse") = as.character((seqnames(gr)))

# Combine all files into one GRangelist
peak <- GRangesList("U-7" = peak1_f, "U-8" = peak2_f, "U-9" = peak3_f, "U-10" = peak4_f, "U-11" = peak5_f,
                    "U-12" = peak6_f, "S-16" = peak7_f, "S-14" = peak8_f, "S-13" = peak9_f,
                    "S-15" = peak10_f, "Q-1" = peak11_f, "Q-2" = peak12_f, "Q-3" = peak13_f,
                    "Q-4" = peak14_f, "Q-5" = peak15_f, "Q-6" = peak16_f)

# remove the unnecessary files
remove(peak1, peak2, peak3, peak4, peak5, peak6, peak7, peak8, peak9, peak10, peak11, peak12, peak13, peak14, peak15, peak16)

# Coverage plots
covplot(peak1_f, weightCol = 5)
covplot(peak2_f, weightCol = 5)
covplot(peak3_f, weightCol = 5)
covplot(peak4_f, weightCol = 5)
covplot(peak5_f, weightCol = 5)
covplot(peak6_f, weightCol = 5)
covplot(peak7_f, weightCol = 5)
covplot(peak8_f, weightCol = 5)
covplot(peak9_f, weightCol = 5)
covplot(peak10_f, weightCol = 5)
covplot(peak11_f, weightCol = 5)
covplot(peak12_f, weightCol = 5)
covplot(peak13_f, weightCol = 5)
covplot(peak14_f, weightCol = 5)
covplot(peak15_f, weightCol = 5)
covplot(peak16_f, weightCol = 5)
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
tagMatrix10 <- getTagMatrix(peak10_f, windows=promoter1)
tagMatrix11 <- getTagMatrix(peak11_f, windows=promoter1)
tagMatrix12 <- getTagMatrix(peak12_f, windows=promoter1)
tagMatrix13 <- getTagMatrix(peak13_f, windows=promoter1)
tagMatrix14 <- getTagMatrix(peak14_f, windows=promoter1)
tagMatrix15 <- getTagMatrix(peak15_f, windows=promoter1)
tagMatrix16 <- getTagMatrix(peak16_f, windows=promoter1)

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
tagHeatmap(tagMatrix10, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix11, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix12, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix13, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix14, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix15, xlim = c(-3000, 3000), color="red")
tagHeatmap(tagMatrix16, xlim = c(-3000, 3000), color="red")

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
peakAnno10 <- annotatePeak(peak10_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno11 <- annotatePeak(peak11_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno12 <- annotatePeak(peak12_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno13 <- annotatePeak(peak13_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno14 <- annotatePeak(peak14_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno15 <- annotatePeak(peak15_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)
peakAnno16 <- annotatePeak(peak16_f, tssRegion=c(-3000, 3000),
                          TxDb=txdb)

# Annotate peaks for all pairs
peakAnnoList <- lapply(peak, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

# Find overlap of peaks and annotated genes
genes_o= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes_o)

#Get intersecting genes between all samples and all except the clinical sample
common_1_2 <- intersect(genes_o[[1]], genes_o[[2]])
common_3_4 <- intersect(genes_o[[3]], genes_o[[4]])
common_5_6 <- intersect(genes_o[[5]], genes_o[[6]])
common_7_8 <- intersect(genes_o[[7]], genes_o[[8]])
common_9_10 <- intersect(genes_o[[9]], genes_o[[10]])
common_11_12 <- intersect(genes_o[[11]], genes_o[[12]])
common_13_14 <- intersect(genes_o[[13]], genes_o[[14]])
common_15_16 <- intersect(genes_o[[15]], genes_o[[16]])
com_1 <- intersect(common_1_2, common_3_4)
com_2 <- intersect(common_5_6, common_7_8)
com_3 <- intersect(common_9_10, common_11_12)
com_4 <- intersect(common_13_14, common_15_16)
co_1 <- intersect(com_1, com_2)
co_2 <- intersect(com_3, com_4)
common <- intersect(co_1, co_2)

common_df <- bitr(common, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = "org.Hs.eg.db")
fwrite(common_df, "Common_genes_16_samples.txt", sep = "\t")


# make df from anno file and select the columns: 
vis_1 <- as.data.frame(peakAnno1) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_1 <- vis_1[grepl("Promoter*", vis_1$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_1) <- 1:nrow(vis_1)
colnames(vis_1) <- c("U-7", "genes")
vis_1$genes <- as.numeric(vis_1$genes)
vis_1 <- ddply(vis_1,"genes",numcolwise(sum))
vis_2 <- as.data.frame(peakAnno2) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_2 <- vis_2[grepl("Promoter*", vis_2$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_2) <- 1:nrow(vis_2)
colnames(vis_2) <- c("U-8", "genes")
vis_2$genes <- as.numeric(vis_2$genes)
vis_2 <- ddply(vis_2,"genes",numcolwise(sum))
vis_3 <- as.data.frame(peakAnno3) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_3 <- vis_3[grepl("Promoter*", vis_3$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_3) <- 1:nrow(vis_3)
colnames(vis_3) <- c("U-9","genes")
vis_3$genes <- as.numeric(vis_3$genes)
vis_3 <- ddply(vis_3,"genes",numcolwise(sum))
vis_4 <- as.data.frame(peakAnno4) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_4 <- vis_4[grepl("Promoter*", vis_4$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_4) <- 1:nrow(vis_4)
colnames(vis_4) <- c("U-10", "genes")
vis_4$genes <- as.numeric(vis_4$genes)
vis_4 <- ddply(vis_4,"genes",numcolwise(sum))
vis_5 <- as.data.frame(peakAnno5) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_5 <- vis_5[grepl("Promoter*", vis_5$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_5) <- 1:nrow(vis_5)
colnames(vis_5) <- c("U-11", "genes")
vis_5$genes <- as.numeric(vis_5$genes)
vis_5 <- ddply(vis_5,"genes",numcolwise(sum))
vis_6 <- as.data.frame(peakAnno6) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_6 <- vis_6[grepl("Promoter*", vis_6$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_6) <- 1:nrow(vis_6)
colnames(vis_6) <- c("U-12", "genes")
vis_6$genes <- as.numeric(vis_6$genes)
vis_6 <- ddply(vis_6,"genes",numcolwise(sum))
vis_7 <- as.data.frame(peakAnno7) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_7 <- vis_7[grepl("Promoter*", vis_7$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_7) <- 1:nrow(vis_7)
colnames(vis_7) <- c("S-16", "genes")
vis_7$genes <- as.numeric(vis_7$genes)
vis_7 <- ddply(vis_7,"genes",numcolwise(sum))
vis_8 <- as.data.frame(peakAnno8) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_8 <- vis_8[grepl("Promoter*", vis_8$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_8) <- 1:nrow(vis_8)
colnames(vis_8) <- c("S-14","genes")
vis_8$genes <- as.numeric(vis_8$genes)
vis_8 <- ddply(vis_8,"genes",numcolwise(sum))
vis_9 <- as.data.frame(peakAnno9) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_9 <- vis_9[grepl("Promoter*", vis_9$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_9) <- 1:nrow(vis_9)
colnames(vis_9) <- c("S-13", "genes")
vis_9$genes <- as.numeric(vis_9$genes)
vis_9 <- ddply(vis_9,"genes",numcolwise(sum))
vis_10 <- as.data.frame(peakAnno10) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_10 <- vis_10[grepl("Promoter*", vis_10$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_10) <- 1:nrow(vis_10)
colnames(vis_10) <- c("S-15", "genes")
vis_10$genes <- as.numeric(vis_10$genes)
vis_10 <- ddply(vis_10,"genes",numcolwise(sum))
vis_11 <- as.data.frame(peakAnno11) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_11 <- vis_11[grepl("Promoter*", vis_11$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_11) <- 1:nrow(vis_11)
colnames(vis_11) <- c("Q-1", "genes")
vis_11$genes <- as.numeric(vis_11$genes)
vis_11 <- ddply(vis_11,"genes",numcolwise(sum))
vis_12 <- as.data.frame(peakAnno12) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_12 <- vis_12[grepl("Promoter*", vis_12$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_12) <- 1:nrow(vis_12)
colnames(vis_12) <- c("Q-2", "genes")
vis_12$genes <- as.numeric(vis_12$genes)
vis_12 <- ddply(vis_12,"genes",numcolwise(sum))
vis_13 <- as.data.frame(peakAnno13) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_13 <- vis_13[grepl("Promoter*", vis_13$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_13) <- 1:nrow(vis_13)
colnames(vis_13) <- c("Q-3","genes")
vis_13$genes <- as.numeric(vis_13$genes)
vis_13 <- ddply(vis_13,"genes",numcolwise(sum))
vis_14 <- as.data.frame(peakAnno14) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_14 <- vis_14[grepl("Promoter*", vis_14$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_14) <- 1:nrow(vis_14)
colnames(vis_14) <- c("Q-4", "genes")
vis_14$genes <- as.numeric(vis_14$genes)
vis_14 <- ddply(vis_14,"genes",numcolwise(sum))
vis_15 <- as.data.frame(peakAnno15) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_15 <- vis_15[grepl("Promoter*", vis_15$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_15) <- 1:nrow(vis_15)
colnames(vis_15) <- c("Q-5", "genes")
vis_15$genes <- as.numeric(vis_15$genes)
vis_15 <- ddply(vis_15,"genes",numcolwise(sum))
vis_16 <- as.data.frame(peakAnno16) %>% dplyr::select(1, 2, 3, 7, 13, 19)
vis_16 <- vis_16[grepl("Promoter*", vis_16$annotation),] %>% dplyr::select(4, 6) %>% dplyr::distinct()
row.names(vis_16) <- 1:nrow(vis_16)
colnames(vis_16) <- c("Q-6", "genes")
vis_16$genes <- as.numeric(vis_16$genes)
vis_16 <- ddply(vis_16,"genes",numcolwise(sum))




# Create a smaller list by getting only those genes that were present in both RNASeq and ATACSeq
vis <- list(vis_1, vis_2, vis_3, vis_4, vis_5,
            vis_6, vis_7, vis_8, vis_9, vis_10,
            vis_11, vis_12, vis_13, vis_14, vis_15,
            vis_16) %>% purrr::reduce(full_join, by = "genes") %>% distinct_at(1, .keep_all = T)
colnames(vis)[1] <- "ENTREZID"
vis$ENTREZID <- as.character(vis$ENTREZID)
vis <- inner_join(vis, common_df, by = "ENTREZID")
vis <- vis[, c(1, 18, 19, 2:17)] %>% distinct()

# Simplify colnames for vis and vis_3pairs
# colnames(vis)[3:5] <- c("MFa2-2", "MFa1-2", "MFa1-4")

# intersect with inflammatory genes from msigdb
vis <- inner_join(vis, gene_set1, by = "ENTREZID")
vis <- vis %>% dplyr::select(-ncol(vis))



# Find the range for the numeric value to plot and then divide it into parts for each color
# Plot Heatmaps for the common genes
split1 <- colorRampPalette(c("red"))(n=49)
split2 <- colorRampPalette(c("green", "blue"), bias = 1.5)(n=250)
my_palette <- c(rev(split2), split1)
colors1 = unique(c(seq(14, 100, length=100), seq(100, 400, length=100), seq(400, 1076, length=102)))
colors2 = unique(c(seq(14, 100, length=100), seq(100, 450, length=100), seq(450, 644, length=102)))


# Prepare for clustering
mydata1 <- as.matrix(vis[,c(4:ncol(vis))])
rownames(mydata1) <- vis[, 3]     # Rownames as gene symbols   
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


