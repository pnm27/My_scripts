

a <- fread("../common_ReQTL/sc_ol_GTEx_sQTL_43tissues/ol_10/N7_10_VAF_matrix_bar_harmonized.txt")

b <- fread("N7_adip_10_bc.txt", header = F)
c <- fread("N7_ery_10_bc.txt", header = F)
d <- fread("N7_naive_10_bc.txt", header = F)

b_h <- colnames(a) %in% b$V1
c_h <- colnames(a) %in% c$V1
d_h <- colnames(a) %in% d$V1
b_h[1] <- T
c_h[1] <- T
d_h[1] <- T

b_a <- a[, c(b_h)]
c_a <- a[, c(c_h)]
d_a <- a[, c(d_h)]


fwrite(b_a, "N7_adip_10_VAF_matrix_GTEX.txt", sep = "\t")
fwrite(c_a, "N7_ery_10_VAF_matrix_GTEX.txt", sep = "\t")
fwrite(d_a, "N7_naive_10_VAF_matrix_GTEX.txt", sep = "\t")
