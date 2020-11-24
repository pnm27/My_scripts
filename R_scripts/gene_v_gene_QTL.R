

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("MatrixEQTL"); load_package("tidyverse"); load_package("data.table")

setwd("D:/Virtualbox/Ubuntu/Shared Folder/")

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel <<- modelLINEAR

# Genotype and gene location file names
SNV_file_name <<- "Heart_LV_GE_q_transf.txt"
snvs_location_file_name <<- "gene_loc_OK.txt"

# Gene expression file name
expression_file_name <<- "Heart_LV_GE_q_transf.txt"
gene_location_file_name <<- "gene_loc_OK.txt"

# Covariates file name
covariates_file_name <<- "Heart_LV_cov.txt"

# Whether to split into cis and trans
split_cis_trans <<- "F"

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <<- numeric()

# Filename for qqplot
qqplot_filename <<-  "Heart_LV_female_ReQTL_combo.tiff"

if(split_cis_trans == "T") {
  # Output file name
  output_file_name_cis <<- "Heart_LV_female_cis_ReQTL.txt"
  output_file_name_tra <<- "Heart_LV_female_trans_ReQTL.txt"
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis <<- 1e-27
  pvOutputThreshold_tra <<- 1e-27
  
  # Distance for local gene-SNV pairs
  cisDist <<- 1e6
  
} else {
  # Output file name
  output_file_name <<- "Heart_LV_female_combo_ReQTL.txt"
  
  # Only associations significant at this level will be saved
  pvOutputThreshold <<- 1e-27
}


## Load genotype data
snvs <- SlicedData$new()
snvs$fileDelimiter = "\t"      # the TAB character
snvs$fileOmitCharacters = "NA" # denote missing values
snvs$fileSkipRows = 1          # one row of column labels
snvs$fileSkipColumns = 1       # one column of row labels
snvs$fileSliceSize = 2000      # read file in slices of 2,000 rows
snvs$LoadFile(SNV_file_name)

## Load gene expression data

gene <- SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt <- SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name) > 1) {
  cvrt$LoadFile(covariates_file_name)
}

## Run the analysis
# snvspos <- read.table(snvs_location_file_name, header = TRUE, stringsAsFactors = FALSE)
# genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
snvspos <- fread(snvs_location_file_name, header = T, select = 1:3)
genepos <- fread(gene_location_file_name, header = T)


if(split_cis_trans == "T") {
  me <- Matrix_eQTL_main(
    snps = snvs, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold  = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snvspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
} else {
  me <- Matrix_eQTL_main(
    snps = snvs, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = output_file_name,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    pvOutputThreshold= pvOutputThreshold,
    snpspos = snvspos, 
    genepos = genepos,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
}

tiff(filename = qqplot_filename)
plot(me)
dev.off()

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
