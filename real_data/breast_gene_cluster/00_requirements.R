
## MODIFY ROOT PATH.
ROOT_PATH =  "/breast_gene_expression_cluster/"

if (!require(mvtnorm, lib = paste0(ROOT_PATH, "req_lib"))) {
  .libPaths(paste0(ROOT_PATH, "req_lib"))
  install.packages("mvtnorm")
  library(mvtnorm)
}

if(!require(MASS, lib = paste0(ROOT_PATH, "req_lib"))) {
  .libPaths(paste0(ROOT_PATH, "req_lib"))
  install.packages("MASS")
  library(MASS)
}


if(!require(e1071, lib = paste0(ROOT_PATH, "req_lib"))) {
  .libPaths(paste0(ROOT_PATH, "req_lib"))
  install.packages("e1071")
  library(e1071)
}

if (!require(equalCovs, lib = paste0(ROOT_PATH, "req_lib"))) {
  .libPaths(paste0(ROOT_PATH, "req_lib"))
  install.packages('equalCovs_1.0.tar.gz', repos=NULL, type='source')
  library(equalCovs,lib = paste0(ROOT_PATH, "req_lib"))
}

if (!require(e1071,lib = paste0(ROOT_PATH, "req_lib"))) {
  .libPaths(paste0(ROOT_PATH, "req_lib"))
  install.packages('equalCovs_1.0.tar.gz', repos=NULL, type='source')
  library(e1071,lib = paste0(ROOT_PATH, "req_lib"))
}

