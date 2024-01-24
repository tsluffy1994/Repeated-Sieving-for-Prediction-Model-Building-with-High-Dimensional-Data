# Simulation study for survival data, repeated sieving
library(doParallel)
library(foreach)
library(simsurv)
library(survival)
library(modelr)
library(purrr)
library(glmnet)
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
beta <- c(sig1=0.7,sig2=-0.7,sig3=0.7,sig4=-0.7,sig5=0.7,sig6=-0.7)
T1 <- 50
T2 <- 40
iteration <- 100

# Simulation
OUT <- foreach(i=1:20, .combine=rbind) %dopar% {
  library(simsurv)
  library(survival)
  library(modelr)
  library(purrr)
  library(glmnet)
  source("Function_simu.R")
  
  RES <- NULL
  for (h in (5*(i-1)+1):(5*(i-1)+5)){
    set.seed(h)
    out <- simulation_single_surv(n,m,m_sig,m_nonsig,block_size_gen,r,beta,T1,T2,block_size_sieve,0.01,0.02,0.0005,0.001)
    RES <- rbind(RES, out) 
  }
  RES
}

# Summary table
ta <- matrix(colMeans(OUT,na.rm=TRUE),ncol=8)
colnames(ta) <- paste0(c('lasso','en','c_svs_sieve_10','c_svs_sieve_100'),rep(c('_30%','_10%'),each=4))
rownames(ta) <- c('ncov','nts','nlp','nlp.val','cin','cin.val','time')
write.table(ta,file='simu_surv.txt')
