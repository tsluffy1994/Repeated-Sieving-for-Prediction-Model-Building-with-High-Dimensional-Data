# Simulation study for binary data, repeated sieving
library(doParallel)
library(foreach)
library(ROCR)
library(modelr)
library(purrr)
library(glmnet)
library(caret)
source("Function_simu.R")

# Set cluster
ncl <- 20
cl <- makeCluster(ncl)
registerDoParallel(cl)
getDoParWorkers()

# Parameter setting
n <- 400
m <- 10000
m_sig <- 6
m_nonsig <- m-m_sig
block_size_gen <- 20
r <- 0.2
block_size_sieve <- 50
beta <- (-1)^(1:(m_sig+1))*2
iteration <- 100

# Simulation
OUT <- foreach(i=1:20, .combine=rbind) %dopar% {
  library(ROCR)
  library(modelr)
  library(purrr)
  library(glmnet)
  library(caret)
  source("Function_simu.R")
  
  RES <- NULL
  for (h in (5*(i-1)+1):(5*(i-1)+5)){
    set.seed(h)
    out <- simulation_single_bin(n,m,m_sig,m_nonsig,block_size_gen,r,beta,block_size_sieve,0.01,0.02,0.0025,0.005)
    RES <- rbind(RES, out) 
  }
  RES
}

# Summary table
ta <- matrix(colMeans(OUT),ncol=4)
colnames(ta) <- c('lasso','en','l_svs_sieve_10','l_svs_sieve_100')
rownames(ta) <- c('ncov','nts','auc','auc.val','time')
write.table(ta,file='simu_bin.txt')
