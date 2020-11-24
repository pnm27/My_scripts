# PLOT_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019

# install missing required packages and load packages
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))

setwd("D:/Virtualbox/Ubuntu/Shared Folder/")

# load in the required matrix files
gene_exp_matrix <<- fread("Heart_LV_GE_q_transf.txt") %>% gather(sample, TPM, -gene_id) %>% drop_na()


# load in a results file
results <<- fread("Heart_LV_female_combo_ReQTL_filt.txt") %>% 
  as.data.frame() %>%
  arrange(FDR) %>%
  mutate(pair = paste0(SNP, "_", gene))


# specify plotting mode
plot_mode <<- "bulk" # or single

# if bulk plotting, specify number of top correlations to plot
n_to_plot <<- 200
n_to_plot <<- ifelse(n_to_plot > nrow(results), nrow(results), n_to_plot)


# specify output prefix
output_prefix <<- "Heart_LV_female_filt_red_combo"
  



if(plot_mode == "bulk") {
  

  plot_df <- left_join(left_join(results, gene_exp_matrix, by = c("SNP" = "gene_id")), gene_exp_matrix, by = c("gene" = "gene_id", "sample"))
  # Remove duplicate pairs (present as reversed pairs)
  a <- duplicated(t(apply(plot_df[, c(1:2, 8)], 1, sort)))
  plot_df2 <- plot_df[!a, ]
  ### FOR PLOTTING BULK SNVs
  uniq_pairs <- head(unique(plot_df2$pair), n_to_plot)
  to_plot <- plot_df2[which(plot_df2$pair %in% uniq_pairs)]
  setorder(plot_df3, beta)
  p <- ggscatter(to_plot, x = "TPM.x", y = "TPM.y",
                 fill = "lightsteelblue", color = "lightsteelblue", shape = 21, size = 1.5, # Points color, shape and size
                 add = "reg.line",  # Add regression line
                 add.params = list(color = "dodgerblue4", fill = alpha("dodgerblue4"), 0.5), # Customize reg. line
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coef.size = 5,
                 cor.coeff.args = list(method = "spearman", label.sep = "\n")) +
    labs(x = "Gene1", y = "Gene2")
  
  pdf(paste0("Heart_LV_ReQTL_combo_plots_top", n_to_plot, ".pdf"), width = ifelse(n_to_plot / 4 >= 8, n_to_plot / 4, 8), height = ifelse(n_to_plot / 4 >= 8, n_to_plot / 4, 8))
  suppressWarnings(print(facet(p, facet.by = "pair", ncol = round(sqrt(n_to_plot)), scales = "free_y")))
  garbage <- dev.off()
  cat(paste0("Bulk plot saved to saved to output/", "Heart_LV_ReQTL_combo_plots_top", n_to_plot, ".pdf\n"))
} else if (plot_mode == "single") {
  ### FOR PLOTTING INDIVIDUAL SNVs
  plot_df <- left_join(vaf_matrix %>% filter(SNV == !!my_snv), 
                       splicing_matrix %>% filter(intron == !!my_intron), 
                       by = 'sample')
  
  p <- ggscatter(plot_df, x = "vaf", y = "ratio",
                 fill = "lightsteelblue", color = "lightsteelblue", shape = 21, size = 1.5, # Points color, shape and size
                 add = "reg.line",  # Add regression line
                 add.params = list(color = "dodgerblue4", fill = alpha("dodgerblue4"), 0.5), # Customize reg. line
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coef.size = 5,
                 cor.coeff.args = list(method = "spearman", label.sep = "\n")) +
    labs(x = "Variant allele fraction", y = "Proportion of read spanning intron junction", title = paste0(my_snv, "-", my_intron))
  
  pdf(paste0("output/", paste0(gsub(":", "_", paste0("RsQTL_plot_", my_snv, "-", my_intron, ".pdf")))))
  suppressWarnings(print(p))
  garbage <- dev.off()
  cat(paste0("Single plot saved to output/", gsub(":", "_", paste0("RsQTL_plot_", my_snv, "-", my_intron, ".pdf\n"))))
  
} else stop("ERROR: Invalid mode specified!\n")
