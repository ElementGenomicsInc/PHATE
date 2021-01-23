#library(reticulate)
#virtualenv_create("r-reticulate", "/usr/bin/python3", packages="pip")
#virtualenv_install("r-reticulate", "phate")
#virtualenv_install("r-reticulate", "packaging")
#virtualenv_install("r-reticulate", "magic-impute")
#install.packages("phateR")

library(phateR)
library(ggplot2)
library(Seurat)

get_dists = function(in_path, out_path)
{
  data = Read10X(in_path)
  data = as.data.frame(t(as.matrix(data)))
  
  # keep genes expressed in at least 10 cells
  keep_cols <- colSums(data > 0) > 10
  data <- data[,keep_cols]
  
  # look at the distribution of library sizes
  #ggplot() +
  #  geom_histogram(aes(x=rowSums(data)), bins=50) +
  #  geom_vline(xintercept = 1500, color='red')
  
  # keep cells with at least 1500 UMIs
  keep_rows <- rowSums(data) > 1500
  data <- data[keep_rows,]
  
  data <- library.size.normalize(data)
  data <- sqrt(data)
  
  data_transpose = as.data.frame(t(as.matrix(data)))
  dists <- phate(data_transpose)
  write.table(dists, out_path, sep="\t", row.names=FALSE, col.names=FALSE)
}

get_dists("/data/frozen_pbmc_donor_a/filtered_matrices_mex/hg19", "/data/dists1.tsv")
get_dists("/data/frozen_pbmc_donor_b/filtered_matrices_mex/hg19", "/data/dists2.tsv")
