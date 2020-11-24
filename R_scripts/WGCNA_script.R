library(tidyverse)
library(data.table)
library(WGCNA)
library(org.Hs.eg.db)
library(msigdbr)  # Msig db


setwd("D:/Virtualbox/Ubuntu/Shared Folder/Entcheva/")

gene_exp <- fread("Full-data.txt", na.strings = NULL)
inp <- as.data.table(t(gene_exp))
inp <- inp[-1, ]
rownames(inp) <- colnames(gene_exp)[-1]
colnames(inp) <- unlist(gene_exp[, 1])
inp <- sapply(inp, as.numeric)
inp <- as.data.frame(inp)

rm(gene_exp)

gsg = goodSamplesGenes(inp, verbose = 3)
gsg$allOK
#----------------------------------------------------------------------------------------
# Filter genes that have lots of missing values
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(inp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(inp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  inp = inp[gsg$goodSamples, gsg$goodGenes]
}
#----------------------------------------------------------------------------------------
# Cluster samples to find outliers

sampleTree = hclust(dist(inp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#----------------------------------------------------------------------------------------
# Pick a branch cut corresponding to the outlier samples by picking a sufficient height from the previous tree
# Plot a line to show the cut
branch_cut = 75

abline(h = branch_cut, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = branch_cut, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = inp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#----------------------------------------------------------------------------------------


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#----------------------------------------------------------------------------------------

net = blockwiseModules(datExpr, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
#----------------------------------------------------------------------------------------

# Get a gene set for narrowing down the analysis (The inflammatory genes)
gene_set <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" ) %>% dplyr::filter(gs_id == "M13657") %>% dplyr::select(7, 8)
# gs_for_GO <- as.vector(gene_set[,7])
colnames(gene_set)[1] <- "ENTREZID"
gene_set$ENTREZID <- as.character(gene_set$ENTREZID)
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10)
