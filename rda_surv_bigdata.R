# Real data analysis, survival, repeated sieving
library(survival)
library(modelr)
library(purrr)
library(glmnet)
library(stringr)
source("Function_rda.R")

# Read in data
data <- readRDS('GSE2034.rds')
data <- data[,-3]
block_size_sieve <- 50
set.seed(619)

# Analysis
out_surv <- rda_surv_50hold(data,block_size_sieve,0.0001,0.0002,0.0005,0.001,0.0005,0.001)
out_surv <- data.frame(matrix(out_surv,ncol=4))
colnames(out_surv) <- c('lasso','en','c_svs_sieve_100','c_svs')
rownames(out_surv) <- c('ns','nlp','nlp.val','cin','cin.val','sel','time','n_can')
write.table(out_surv,file='rda_surv_50hold.txt')
