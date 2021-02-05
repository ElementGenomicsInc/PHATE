#library(reticulate)
#virtualenv_create("r-reticulate", "/usr/bin/python3")
#bash commands:
#git clone --recursive https://github.com/ElementGenomicsInc/PHATE.git
#cd PHATE/Python
#sudo /home/ruser/.virtualenvs/r-reticulate/bin/python setup.py install
#cd ../phateR
#sudo R CMD INSTALL .
.libPaths("/usr/local/lib/R/site-library")
library(phateR)

args = commandArgs(trailingOnly = TRUE)
in_path = args[1]
out_path = args[2]

data = read.table(in_path, sep="\t")
data <- library.size.normalize(data)
data <- sqrt(data)

data_transpose = as.data.frame(t(as.matrix(data)))
potentials <- phate(data_transpose, n.landmark = NULL)
write.table(potentials, out_path, sep="\t") 