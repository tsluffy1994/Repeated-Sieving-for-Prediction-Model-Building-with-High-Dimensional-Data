# Real data analysis, binary, repeated sieving
library(ROCR)
library(modelr)
library(purrr)
library(glmnet)
library(stringr)
source("Function_rda.R")

# Read in data
data <- readRDS('GSE2034.rds')
data <- data[,-c(1,2)]
block_size_sieve <- 50
set.seed(53)

# Analysis
out_bin <- rda_bin_50hold(data,block_size_sieve,0.000001,0.000002,0.0025,0.005,0.0025,0.005)
out_bin <- data.frame(matrix(out_bin,ncol=4))
colnames(out_bin) <- c('lasso','en','l_svs_sieve_100','l_svs')
rownames(out_bin) <- c('ns','auc','auc.val','sel','time','n_can')
write.table(out_bin,file='rda_bin_50hold.txt')
