library(ggplot2)
library(Seurat)

potentials_from_counts = function(data, out_path, threshold)
{
  data <- library.size.normalize(data)
  data <- sqrt(data)
  
  data_transpose = as.data.frame(t(as.matrix(data)))
  potentials <- phate(data_transpose, n.landmark = NULL)
  write.table(potentials, out_path, sep="\t") 
}

merge_counts = function(counts1, counts2, counts3){
  genes_intersection = intersect(rownames(counts1), rownames(counts2))
  genes_intersection = intersect(genes_intersection, rownames(counts3))
  return(as.data.frame(t(as.matrix(cbind(counts1[genes_intersection,], counts2[genes_intersection,], counts3[genes_intersection,])))))
}

#counts1 = Read10X("/data/frozen_pbmc_donor_a/filtered_matrices_mex/hg19")
#counts2 = Read10X("/data/frozen_pbmc_donor_b/filtered_matrices_mex/hg19")

#diabetes

#rep 1
counts1 = readRDS("/data/GSM3823942_diabetes.s1.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts1 = counts1[rowSums(counts1 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(x=colSums(counts1)), bins=50) 
threshold = 1000
# keep cells with UMIs meeting threshold
counts1 = counts1[,colSums(counts1) > threshold]

#rep 2
counts2 = readRDS("/data/GSM3823943_diabetes.s2.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts2 = counts2[rowSums(counts2 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(colSums(counts2)), bins=50) 
threshold = 500
# keep cells with UMIs meeting threshold
counts2 = counts2[,colSums(counts2) > threshold]

#rep 3
counts3 = readRDS("/data/GSM3823944_diabetes.s3.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts3 = counts3[rowSums(counts3 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(colSums(counts3)), bins=50) 
threshold = 1000 
# keep cells with UMIs meeting threshold
counts3 = counts3[,colSums(counts3) > threshold]
  
diabetes_counts = merge_counts(counts1, counts2, counts3)
write.table(diabetes_counts, "/data/diabetes_counts.tsv",sep="\t")

potentials_from_counts(diabetes_counts, threshold, "/data/diabetes_diff_potentials.tsv")

#control

#rep 1
counts1 = readRDS("/data/GSM3823939_control.s1.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts1 = counts1[rowSums(counts1 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(colSums(counts1)), bins=50) 
threshold = 2500
# keep cells with UMIs meeting threshold
counts1 = counts1[,colSums(counts1) > threshold]

#rep 2
counts2 = readRDS("/data/GSM3823940_control.s2.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts2 = counts2[rowSums(counts2 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(colSums(counts2)), bins=50) 
threshold = 500
# keep cells with UMIs meeting threshold
counts2 = counts2[,colSums(counts2) > threshold]

#rep 3
counts3 = readRDS("/data/GSM3823941_control.s3.dgecounts.rds")$umicount$exon$all
# keep genes expressed in at least 10 cells
counts3 = counts3[rowSums(counts3 > 0) > 10,]
#find threshold
ggplot() + geom_histogram(aes(colSums(counts3)), bins=50) 
threshold = 1000
# keep cells with UMIs meeting threshold
counts3 = counts3[,colSums(counts3) > threshold]

control_counts = merge_counts(counts1, counts2, counts3)
write.table(control_counts, "/data/control_counts.tsv",sep="\t")
  
potentials_from_counts(control_counts, threshold, "/data/diabetes_diff_potentials.tsv")