# list.of.packages <- c("ggplot2", "Rcpp", "BiocManager" ,"devtools", "data.table", "gtools", "optparse")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

# suppressPackageStartupMessages(library(BiocManager))
# suppressPackageStartupMessages(library(devtools))
# devtools::install_github("hadley/lazyeval")
# devtools::install_github("hadley/dplyr")
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(msigdbr))

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="You will see progress bar now! :)"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Now you won't see the progress bar! :("),
  make_option(c("-d", "--dir"), action="store", type = "character", default=NULL,
              help="provide the directory where the featureCount file is"),
  make_option(c("-c", "--coldata"), action="store", type = "character", default=NULL,
              help="provide the directory where the coldata file is"),
  make_option(c("-n", "--namecontrol"), action="store", type = "character", default="control",
              help="How you want the \"control\" (default: control) to be named"),
  make_option(c("-f", "--fdr"), action="store", type = "numeric", default=0.05,
              help="The FDR cut-off for significance. DEFAULT: 0.05")
)

opt <- parse_args(OptionParser(option_list=option_list, description = "-d and -c options are necessary!!!!"))

if (any(is.null(c(opt$dir, opt$coldata)))) {
  stop("Please make sure you have provided all the files")
}

# Check if file exists
res <- tryCatch({mixedsort(list.files(path = opt$dir, pattern = "*.txt$"), decreasing = F)}, error=function(cond){
  message("No file exists. Please check the address!")
})


# create progress bar
if (opt$verbose) {
  total <- 6
  pb <- txtProgressBar(min = 0, max = total, style = 3)
}
ifelse(opt$verbose, setTxtProgressBar(pb, 1))

# Get a gene set for narrowing down the analysis (The inflammatory genes)
# gene_set1 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP" ) %>% dplyr::filter(gs_id == "M13657") %>% dplyr::select(7, 8)
# colnames(gene_set1)[1] <- "ENTREZID"
# gene_set1$ENTREZID <- as.character(gene_set1$ENTREZID)
# fwrite(gene_set1, "Inflammatory_gens_msgib.txt", sep = "\t")

#File path where there are txt files from featureCounts
file_path <- opt$dir

ls <- list.files(path = file_path, pattern = "*.txt$") # List files in decreasing order so that matrix is formed in increasing order

s <- fread(paste0(file_path,ls))

fn <- list.files(path = opt$coldata, pattern = "*.txt$")
print(paste0(opt$coldata, fn))
cd <- fread(paste0(opt$coldata, fn), select = 1:2,
            col.names = c("sample", "condition"), header = F)

ifelse(opt$verbose, setTxtProgressBar(pb, 2))

coldata <- as.matrix(cd)

s <- na.omit(s)
#cts <- as.matrix(s[,-1], rownames = FALSE) #Make Geneids as the row names

# Here, according to how your counts file looks you should change the marked values
cts <- as.matrix(s[,7:ncol(s)], rownames = FALSE) #
rownames(cts) <- as.character(unlist(s[,1])) # column that has gene_id/gene_name
colnames(cts) <- colnames(s)[7:ncol(s)] #
rownames(coldata) <- as.character(unlist(coldata[, 1]))

# Check the columns vs rows equivalence. Don't go ahead if both of them aren't TRUE!!!!!
print("If any of the following 2 statements is false then please check your inputs!!!!")
# rownames(coldata)
# head(colnames(cts))
# coldata[1:5,]
cts <- cts[,rownames(coldata)]
print(all(rownames(coldata) %in% colnames(cts)))
print(all(rownames(coldata) == colnames(cts)))

# Creating DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

ifelse(opt$verbose, setTxtProgressBar(pb, 3))


#Pre- Filtering
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# Relevel the factors (Set which one is the control sample)
dds$condition <- relevel(dds$condition, ref = opt$namecontrol)

#DE Analysis
dds <- DESeq(dds) # Add "minReplicatesForReplace = Inf" for single-cell Analysis
res <- results(dds)
res_all <- results(dds, alpha = 0.99)

ifelse(opt$verbose, setTxtProgressBar(pb, 4))

# MA plot
plotMA(res, ylim=c(-2,2))

#Log fold change shrinkage
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm") # Set "coef" to the particular comparison you are interested in. coef=2 is equivalent to coef=resultsNames(dds)[2]
resLFC
resOrdered <- res[order(res$pvalue),]

ifelse(opt$verbose, setTxtProgressBar(pb, 5))

# Summary
print("SUMMARY------------------")
print(summary(res))

print("-------------------------")
# q < 0.1
print(sum(res$padj < opt$fdr, na.rm=TRUE))

vsd <- vst(dds, blind=FALSE) # Will remove batch effects
corrs <- cor(assay(vsd), method="spearman")
corr.dists <- as.dist(1 - corrs)
colors <- colorRampPalette(c("red","white","blue"))(99)
png("HM_plot.png")
pheatmap(corrs, breaks=seq(from=-1,to=1,length=100),
         clustering_distance_rows=corr.dists,
         clustering_distance_cols=corr.dists,
         col=colors)
dev.off()

# MA plot after removing genes with low counts
png("MA_plot.png")
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# png("HM_top_var_genes_sample.png", width = 764, height = 800)
rld <- rlog( dds )
# dimnames(assay(rld, withDimnames = T))
topVarGenes <- head(order(-rowVars(assay(rld))),50)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="both",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

dev.off()

# Plot normalized counts after adding a pseudocount of 1/2
# Plot the gene which had the smallest p-value from the results table
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", # "intgroup" specifies the grouping variables(s)
                returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))


# Subset genes with padj < 0.1
resSig <- subset(resOrdered, padj < opt$fdr)
final_out <- as.data.table(resSig, keep.rownames = T)
colnames(final_out)[1] <- "Gene_ID"
data.table::fwrite(x=final_out, "DESEq2_results.txt", sep = "\t")
Allres <- as.data.table(res_all, keep.rownames = T)
data.table::fwrite(Allres, "Pos_v_Neg_fC_gn_all.txt", sep = '\t')
ifelse(opt$verbose, setTxtProgressBar(pb, 6))

# head(counts(dds, normalized=T))
indiv_norm_values <- as.matrix(counts(dds, normalized=T))
temp <- as.data.frame(Gene_ID = row.names(indiv_norm_values), indiv_norm_values)
data.table::fwrite(x=temp, "Indiv_norm_results.txt", sep = "\t", row.names = T)

close(pb)