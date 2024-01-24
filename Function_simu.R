# Functions for repeated sieving prediction, 6 true covs

######################
# Forward stepwise selection based on p-value
# Binary outcome
fit_for1 <- function(y,x,i){
  fit <- glm(y~x[,i],family = "binomial")
  p <- coef(summary(fit))[2,4]
  return(p)
}
fit_for2 <- function(y,x,enter_index,i){
  data_for <- cbind(y,x[,c(enter_index,i)])
  fit_for <- glm(y~.,data=data_for,family = "binomial")
  p <- tail(coef(summary(fit_for)),1)[4]
  return(p)
}
fit_for3 <- function(y,x,enter_index,i){
  data_for2 <- cbind(y,x[,c(enter_index,i)])
  fit_for2 <- glm(y~.,data=data_for2,family = "binomial")
  z <- abs(tail(coef(summary(fit_for2)),1)[3])
  return(z)
}

stepwiseF <- function(data,alpha_E,alpha_R){
  y <- data[,1]
  x <- data[,-1]
  enter_index <- NULL
  out_index <- 1:dim(x)[2]
  
  # First Step
  p <- sapply(out_index,fit_for1,y=y,x=x)
  if (min(p)<alpha_E){
    ind <- which(p==min(p))
    enter_index <- c(enter_index,ind)
    out_index <- out_index[out_index!=ind]
    print(paste0('covariate ',colnames(x)[ind],' added, with p-value ', min(p)))
  } else {
    ind <- which(p==min(p))
    #print(paste0('No covariate pass alpha_entry, select the covariate with smallest p-value, covariate ',colnames(x)[ind],', with p-value ', min(p),' and stop'))
    print('No covariate pass alpha_entry')
    return(NULL)
    #return(colnames(x)[ind])
  }
  
  # Continue
  while (min(p)<alpha_E & min(p)>0){
    # Forward
    p <- sapply(out_index,fit_for2,y=y,x=x,enter_index=enter_index)
    if (min(p)<alpha_E){
      ind <- out_index[which(p==min(p))]
      if (length(ind)>1){
        z <- sapply(ind,fit_for3,y=y,x=x,enter_index=enter_index)
        ind <- ind[z==max(z)]
      }
      enter_index <- c(enter_index,ind)
      out_index <- out_index[out_index!=ind]
      print(paste0('covariate ',colnames(x)[ind],' added, with p-value ', min(p)))
      
      # Backward check
      data_bac <- cbind(y,x[,enter_index])
      fit_bac <- glm(y~.,data_bac,family = "binomial")
      p_bac <- coef(summary(fit_bac))[-1,4]
      if (max(p_bac)>=alpha_R){
        ind_out <- enter_index[p_bac >= alpha_R]
        enter_index <- enter_index[!enter_index%in%ind_out]
        out_index <- c(out_index,ind_out)
        for (ind in ind_out){
          print(paste0('covariate ',colnames(x)[ind],' excluded, with p-value ', p_bac[names(p_bac)==colnames(x)[ind]]))
        }
      }
    } 
  }
  
  cov <- colnames(x)[enter_index]
  
  return(cov)
}

# Survival outcome
fit_for1_surv <- function(y,x,i){
  fit <- coxph(Surv(y$eventtime,y$status)~x[,i])
  p <- coef(summary(fit))[5]
  return(p)
}
fit_for2_surv <- function(y,x,enter_index,i){
  data_for <- cbind(y,x[,c(enter_index,i)])
  fit_for <- coxph(Surv(eventtime,status)~.,data=data_for)
  p <- tail(coef(summary(fit_for)),1)[5]
  return(p)
}
fit_for3_surv <- function(y,x,enter_index,i){
  data_for2 <- cbind(y,x[,c(enter_index,i)])
  fit_for2 <- coxph(Surv(eventtime,status)~.,data=data_for2)
  z <- abs(tail(coef(summary(fit_for2)),1)[4])
  return(z)
}

stepwiseF_surv <- function(data,alpha_E,alpha_R){
  y <- data[,c(1,2)]
  x <- data[,-c(1,2)]
  enter_index <- NULL
  out_index <- 1:dim(x)[2]
  
  # First Step
  p <- sapply(out_index,fit_for1_surv,y=y,x=x)
  if (min(p)<alpha_E){
    ind <- which(p==min(p))
    enter_index <- c(enter_index,ind)
    out_index <- out_index[out_index!=ind]
    print(paste0('covariate ',colnames(x)[ind],' added, with p-value ', min(p)))
  } else {
    ind <- which(p==min(p))
    #print(paste0('No covariate pass alpha_entry, select the covariate with smallest p-value, covariate ',colnames(x)[ind],', with p-value ', min(p),' and stop'))
    print('No covariate pass alpha_entry')
    return(NULL)
    #return(colnames(x)[ind])
  }
  
  # Continue
  while (min(p)<alpha_E & min(p)>0){
    # Forward
    p <- sapply(out_index,fit_for2_surv,y=y,x=x,enter_index=enter_index)
    if (min(p)<alpha_E){
      ind <- out_index[which(p==min(p))]
      if (length(ind)>1){
        z <- sapply(ind,fit_for3_surv,y=y,x=x,enter_index=enter_index)
        ind <- ind[z==max(z)]
      }
      enter_index <- c(enter_index,ind)
      out_index <- out_index[out_index!=ind]
      print(paste0('covariate ',colnames(x)[ind],' added, with p-value ', min(p)))
      
      # Backward check
      data_bac <- cbind(y,x[,enter_index])
      fit_bac <- coxph(Surv(eventtime,status)~.,data=data_bac)
      p_bac <- coef(summary(fit_bac))[,5]
      if (max(p_bac)>=alpha_R){
        ind_out <- enter_index[p_bac >= alpha_R]
        enter_index <- enter_index[!enter_index%in%ind_out]
        out_index <- c(out_index,ind_out)
        for (ind in ind_out){
          print(paste0('covariate ',colnames(x)[ind],' excluded, with p-value ', p_bac[names(p_bac)==colnames(x)[ind]]))
        }
      }
    } 
  }
  
  cov <- colnames(x)[enter_index]
  
  return(cov)
}

######################
# Binary Outcome
# Data generation_binary
data_gen_bin <- function(n,m,m_sig,m_nonsig,block_size,r,beta){
  # Simulate covariates and outcomes
  x <- matrix(0,nrow=n,ncol=m)
  x_block <- matrix(0,nrow=n,ncol=block_size+1)
  for (i in 1:(m/block_size)){
    for (j in 1:ncol(x_block)){
      x_block[,j] <- rnorm(n)
      if (j>1){
        x[,block_size*(i-1)+j-1] <- sqrt(r)*x_block[,1]+sqrt(1-r)*x_block[,j]
      }
    }
  }
  # Scenario 3: Select each 2 covariates from the same block
  ind <- c(seq(1,1+(m_sig/2-1)*block_size,block_size),seq(2,2+(m_sig/2-1)*block_size,block_size))
  ind <- ind[order(ind)]
  x_sig <- x[,ind]
  colnames(x) <- 1:ncol(x)
  colnames(x)[ind] <- paste0("sig",1:m_sig)
  colnames(x)[-ind] <- paste0("nonsig",1:m_nonsig)
  
  # Generate outcomes
  logit <- c(cbind(1,x_sig) %*% beta)
  prob <- exp(logit)/(1+exp(logit))
  y <- rbinom(n,1,prob)
  data <- data.frame(cbind(y,x))
  
  return(data)
}

# Permutation and partition
permutation_partition <- function(data,block_size){
  y <- data[,1]
  x <- data[,-1]
  block <- sample(rep(1:(ncol(x)/block_size),each=block_size))
  x_df <- data.frame(t(rbind(x,block)))
  x_pp <- split(x_df,x_df[,ncol(x_df)])
  data_pp <- lapply(x_pp,function(df) data.frame(cbind(y,t(df[,-ncol(df)]))))
  
  return(data_pp)
}

# R-SVS with Sieving
sieve <- function(data.train,data.test,block_size,P,m_sig,alpha1,alpha2,alpha3,alpha4){
  # Sieve with permutation
  cov_sieve <- NULL
  t0 <- Sys.time()
  for (p in 1:P){
    data_pp <- permutation_partition(data.train,block_size)
    cov <- unlist(map(data_pp, stepwiseF, alpha_E=alpha1, alpha_R=alpha2))
    cov_sieve <- unique(c(cov_sieve,cov))
  }
  n_can <- length(cov_sieve)
  cov_final <- stepwiseF(data.train[,c('y',cov_sieve)],alpha3,alpha4)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  
  # Selections
  sig <- paste0("sig",1:m_sig)
  n_can_sig <- sum(cov_sieve %in% sig)
  ns <- length(cov_final)
  nts <- sum(cov_final %in% sig)
  
  # AUC
  fit <- glm(y~.,data=data.train[,c('y',cov_final)],family='binomial')
  auc.train <- performance(prediction(predict(fit),data.train$y),measure='auc')@y.values[[1]]
  auc.val <- performance(prediction(predict(fit,newdata=data.test),data.test$y),measure='auc')@y.values[[1]]
  
  return(c(ns,nts,auc.train,auc.val,time.step,n_can,n_can_sig))
}

# LASSO
lasso <- function(data.train,data.test,m_sig){
  y.train <- data.train[,1]
  x.train <- data.train[,-1]
  y.test <- data.test[,1]
  x.test <- data.test[,-1]
  
  # LASSO
  t0 <- Sys.time()
  model <- model.matrix(~., data = x.train)[,-1]
  lasFit <- cv.glmnet(x=model, y=y.train, family = "binomial",alpha=1,nfolds = 10)
  model.min <- glmnet(x=model, y=y.train, family = "binomial", alpha = 1,lambda=lasFit$lambda.min)
  time.lasso <- as.numeric(Sys.time()-t0,units='secs')
  cov.lasso <- rownames(coef(model.min))[coef(model.min)[,1]!=0][-1]
  ncov.lasso <- length(cov.lasso)
  sig <- paste0("sig",1:m_sig)
  nts.lasso <- sum(cov.lasso %in% sig)
  auc.lasso <- performance(prediction(predict(model.min,as.matrix(x.train)),y.train),measure='auc')@y.values[[1]]
  auc.val.lasso <- performance(prediction(predict(model.min,as.matrix(x.test)),y.test),measure='auc')@y.values[[1]]
  
  return (c(ncov.lasso,nts.lasso,auc.lasso,auc.val.lasso,time.lasso))
}

# EN
en <- function(data.train,data.test,m_sig){
  y.train <- data.train[,1]
  x.train <- data.train[,-1]
  y.test <- data.test[,1]
  x.test <- data.test[,-1]
  
  # EN
  data.train$y <- as.factor(data.train$y)
  t0 <- Sys.time()
  model.en <- train(
    y~., data = data.train, method = "glmnet",
    trControl = trainControl("cv", number = 10),
    tuneLength = 10
  )
  time.en <- as.numeric(Sys.time()-t0,units='secs')
  cov.en <- rownames(coef(model.en$finalModel, model.en$bestTune$lambda))[coef(model.en$finalModel, model.en$bestTune$lambda)[,1]!=0][-1]
  ncov.en <- length(cov.en)
  sig <- paste0("sig",1:m_sig)
  nts.en <- sum(cov.en %in% sig)
  auc.en <- performance(prediction(predict(model.en,type='prob')[,2],y.train),measure='auc')@y.values[[1]]
  auc.val.en <- performance(prediction(predict(model.en,newdata=x.test,type='prob')[,2],y.test),measure='auc')@y.values[[1]]
  
  return (c(ncov.en,nts.en,auc.en,auc.val.en,time.en))
}

# Single simulation/iteration binary
simulation_single_bin <- function(n,m,m_sig,m_nonsig,block_size_gen,r,beta,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Data generation
  data <- data_gen_bin(n,m,m_sig,m_nonsig,block_size_gen,r,beta)
  
  # 50-50 hold-out
  index <- sample(1:n,size=n*0.5)
  data.test <- data[index,]
  data.train <- data[-index,]
  
  # Models
  out_step_10 <- sieve(data.train,data.test,block_size_sieve,10,m_sig,alpha1,alpha2,alpha3,alpha4)
  out_step_10 <- out_step_10[-c(6,7)]
  out_step_100 <- sieve(data.train,data.test,block_size_sieve,100,m_sig,alpha1,alpha2,alpha3,alpha4)
  out_step_100 <- out_step_100[-c(6,7)]
  out_lasso <- lasso(data.train,data.test,m_sig)
  out_en <- en(data.train,data.test,m_sig)

  # Combine results
  res <- c(out_lasso,out_en,out_step_10,out_step_100)
  return(res)
}

# Single simulation/iteration binary, sieve only
simulation_single_bin_sieve <- function(n,m,m_sig,m_nonsig,block_size_gen,r,beta,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Data generation
  data <- data_gen_bin(n,m,m_sig,m_nonsig,block_size_gen,r,beta)
  
  # 50-50 hold-out
  index <- sample(1:n,size=n*0.5)
  data.test <- data[index,]
  data.train <- data[-index,]
  
  # Models
  out_step_10 <- sieve(data.train,data.test,block_size_sieve,10,m_sig,alpha1,alpha2,alpha3,alpha4)
  out_step_100 <- sieve(data.train,data.test,block_size_sieve,100,m_sig,alpha1,alpha2,alpha3,alpha4)
  
  # Combine results
  res <- c(out_step_10,out_step_100)
  return(res)
}

######################
# Survival Outcome
# Data generation_survival
data_gen_surv <- function(n,m,m_sig,m_nonsig,block_size,r,beta,T1,T2){
  # Simulate covariates and outcomes
  x <- matrix(0,nrow=n,ncol=m)
  x_block <- matrix(0,nrow=n,ncol=block_size+1)
  for (i in 1:(m/block_size)){
    for (j in 1:ncol(x_block)){
      x_block[,j] <- rnorm(n)
      if (j>1){
        x[,block_size*(i-1)+j-1] <- sqrt(r)*x_block[,1]+sqrt(1-r)*x_block[,j]
      }
    }
  }
  # Scenario 3: Select each 2 covariates from the same block
  ind <- c(seq(1,1+(m_sig/2-1)*block_size,block_size),seq(2,2+(m_sig/2-1)*block_size,block_size))
  ind <- ind[order(ind)]
  x_sig <- x[,ind]
  colnames(x) <- 1:ncol(x)
  colnames(x)[ind] <- paste0("sig",1:m_sig)
  colnames(x)[-ind] <- paste0("nonsig",1:m_nonsig)
  
  # Generate outcomes
  colnames(x_sig) <- paste0("sig",1:m_sig)
  covs <- data.frame(id=1:n,x_sig)
  y <- simsurv(dist='exponential',lambdas=0.1,x=covs,betas=beta)[,2]
  cen1 <- runif(n,0,T1)
  cen2 <- runif(n,T2,T1+T2)
  y1 <- data.frame(eventtime=pmin(y,cen1),status=ifelse(y<cen1,1,0))
  y2 <- data.frame(eventtime=pmin(y,cen2),status=ifelse(y<cen2,1,0))
  data1 <- data.frame(cbind(y1,x))
  data2 <- data.frame(cbind(y2,x))
  data <- list(data1,data2)
  
  return(data)
}

# Permutation and partition
permutation_partition_surv <- function(data,block_size){
  y <- data[,c(1,2)]
  x <- data[,-c(1,2)]
  block <- sample(rep(1:(ncol(x)/block_size),each=block_size))
  x_df <- data.frame(t(rbind(x,block)))
  x_pp <- split(x_df,x_df[,ncol(x_df)])
  data_pp <- lapply(x_pp,function(df) data.frame(cbind(y,t(df[,-ncol(df)]))))
  
  return(data_pp)
}

# R-SVS with Sieving
sieve_surv <- function(data.train,data.test,block_size,P,m_sig,alpha1,alpha2,alpha3,alpha4){
  # Sieve with permutation
  cov_sieve <- NULL
  t0 <- Sys.time()
  for (p in 1:P){
    data_pp <- permutation_partition_surv(data.train,block_size)
    cov <- unlist(map(data_pp, stepwiseF_surv, alpha_E=alpha1, alpha_R=alpha2))
    cov_sieve <- unique(c(cov_sieve,cov))
  }
  n_can <- length(cov_sieve)
  cov_final <- stepwiseF_surv(data.train[,c('eventtime','status',cov_sieve)],alpha3,alpha4)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  
  # Selections
  sig <- paste0("sig",1:m_sig)
  n_can_sig <- sum(cov_sieve %in% sig)
  ns <- length(cov_final)
  nts <- sum(cov_final %in% sig)
  
  # C-index and -log10 p-value
  fit <- coxph(Surv(eventtime,status)~.,data=data.train[,c('eventtime','status',cov_final)])
  linear.step.x <- predict(fit)
  linear.step.x.val <- predict(fit,newdata=data.test)
  s.step <- summary(coxph(Surv(data.train$eventtime,data.train$status)~linear.step.x))
  s.step.val <- summary(coxph(Surv(data.test$eventtime,data.test$status)~linear.step.x.val))
  nlp.step <- -log10(coef(s.step)[5])
  nlp.val.step <- -log10(coef(s.step.val)[5])
  cin.step <- s.step$concordance[1]
  cin.val.step <- s.step.val$concordance[1]
  
  return(c(ns,nts,nlp.step,nlp.val.step,cin.step,cin.val.step,time.step,n_can,n_can_sig))
}

# LASSO
lasso_surv <- function(data.train,data.test,m_sig){
  y.train <- data.train[,c(1,2)]
  x.train <- data.train[,-c(1,2)]
  y.test <- data.test[,c(1,2)]
  x.test <- data.test[,-c(1,2)]
  
  # LASSO
  t0 <- Sys.time()
  lasFit <- cv.glmnet(x=as.matrix(x.train), y=Surv(y.train$eventtime,y.train$status),family="cox",alpha=1,nfolds = 10)
  model.min <- glmnet(x=as.matrix(x.train), y=Surv(y.train$eventtime,y.train$status),family="cox",alpha = 1,lambda=lasFit$lambda.min)
  time.lasso <- as.numeric(Sys.time()-t0,units='secs')
  cov.lasso <- rownames(coef(model.min))[coef(model.min)[,1]!=0]
  ncov.lasso <- length(cov.lasso)
  sig <- paste0("sig",1:m_sig)
  nts.lasso <- sum(cov.lasso %in% sig)
  linear.x <- predict(model.min,as.matrix(x.train))
  linear.x.val <- predict(model.min,as.matrix(x.test))
  s <- summary(coxph(Surv(y.train$eventtime,y.train$status)~linear.x))
  s.val <- summary(coxph(Surv(y.test$eventtime,y.test$status)~linear.x.val))
  nlp.lasso <- -log10(coef(s)[5])
  nlp.val.lasso <- -log10(coef(s.val)[5])
  cin.lasso <- s$concordance[1]
  cin.val.lasso <- s.val$concordance[1]
  
  return (c(ncov.lasso,nts.lasso,nlp.lasso,nlp.val.lasso,cin.lasso,cin.val.lasso,time.lasso))
}

# EN
en_surv <- function(data.train,data.test,m_sig){
  y.train <- data.train[,c(1,2)]
  x.train <- data.train[,-c(1,2)]
  y.test <- data.test[,c(1,2)]
  x.test <- data.test[,-c(1,2)]
  
  # EN
  alpha <- seq(0,1,by=0.1)
  lambda <- rep(NA,length(alpha))
  cve <- rep(NA,length(alpha))
  t0 <- Sys.time()
  for (i in 1:length(alpha)){
    model.en <- cv.glmnet(x=as.matrix(x.train), y=Surv(y.train$eventtime,y.train$status),family="cox",alpha=alpha[i],nfolds = 10)
    lambda[i] <- model.en$lambda.min
    cve[i] <- min(model.en$cvm)
  }
  model.min.en <- glmnet(x=as.matrix(x.train), y=Surv(y.train$eventtime,y.train$status),family="cox",alpha=alpha[cve==min(cve)],lambda=lambda[cve==min(cve)])
  time.en <- as.numeric(Sys.time()-t0,units='secs')
  cov.en <- rownames(coef(model.min.en))[coef(model.min.en)[,1]!=0]
  ncov.en <- length(cov.en)
  sig <- paste0("sig",1:m_sig)
  nts.en <- sum(cov.en %in% sig)
  linear.en.x <- predict(model.min.en,as.matrix(x.train))
  linear.en.x.val <- predict(model.min.en,as.matrix(x.test))
  s.en <- summary(coxph(Surv(y.train$eventtime,y.train$status)~linear.en.x))
  s.en.val <- summary(coxph(Surv(y.test$eventtime,y.test$status)~linear.en.x.val))
  nlp.en <- -log10(coef(s.en)[5])
  nlp.val.en <- -log10(coef(s.en.val)[5])
  cin.en <- s.en$concordance[1]
  cin.val.en <- s.en.val$concordance[1]
  
  return (c(ncov.en,nts.en,nlp.en,nlp.val.en,cin.en,cin.val.en,time.en))
}

# Single simulation/iteration survival
simulation_single_surv <- function(n,m,m_sig,m_nonsig,block_size_gen,r,beta,T1,T2,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Data generation
  data <- data_gen_surv(n,m,m_sig,m_nonsig,block_size_gen,r,beta,T1,T2)
  
  # 50-50 hold-out
  index <- sample(1:n,size=n*0.5)
  data.test <- lapply(data,function(da) da[index,])
  data.train <- lapply(data,function(da) da[-index,])
  
  # Models
  out_step_10 <- pmap(list(data.train,data.test),sieve_surv,block_size=block_size_sieve,P=10,m_sig=m_sig,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4)
  out_step_10 <- map(out_step_10,function(da) da[-c(8,9)])
  out_step_100 <- pmap(list(data.train,data.test),sieve_surv,block_size=block_size_sieve,P=100,m_sig=m_sig,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4)
  out_step_100 <- map(out_step_100,function(da) da[-c(8,9)])
  out_lasso <- pmap(list(data.train,data.test),lasso_surv,m_sig=m_sig)
  out_en <- pmap(list(data.train,data.test),en_surv,m_sig=m_sig)
  
  # Combine results
  res <- unlist(Map(c,out_lasso,out_en,out_step_10,out_step_100))
  return(res)
}

# Single simulation/iteration survival, sieve only
simulation_single_surv_sieve <- function(n,m,m_sig,m_nonsig,block_size_gen,r,beta,T1,T2,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Data generation
  data <- data_gen_surv(n,m,m_sig,m_nonsig,block_size_gen,r,beta,T1,T2)
  
  # 50-50 hold-out
  index <- sample(1:n,size=n*0.5)
  data.test <- lapply(data,function(da) da[index,])
  data.train <- lapply(data,function(da) da[-index,])
  
  # Models
  out_step_10 <- pmap(list(data.train,data.test),sieve_surv,block_size=block_size_sieve,P=10,m_sig=m_sig,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4)
  out_step_100 <- pmap(list(data.train,data.test),sieve_surv,block_size=block_size_sieve,P=100,m_sig=m_sig,alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,alpha4=alpha4)
  
  # Combine results
  res <- unlist(Map(c,out_step_10,out_step_100))
  return(res)
}