# Functions for repeated sieving prediction, real data analysis

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
  while (min(p)<alpha_E & min(p)>0 & length(out_index)>0){
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
  while (min(p)<alpha_E & min(p)>0 & length(out_index)>0){
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
# Permutation and partition
permutation_partition <- function(data,block_size){
  y <- data[,1]
  x <- data[,-1]
  block <- sample(rep(1:ceiling(ncol(x)/block_size),each=block_size,length.out=ncol(x)))
  x_df <- data.frame(t(rbind(x,block)))
  x_pp <- split(x_df,x_df[,ncol(x_df)])
  data_pp <- lapply(x_pp,function(df) data.frame(cbind(y,t(df[,-ncol(df)]))))
  
  return(data_pp)
}

# Pure R-SVS
l_svs <- function(data.train,data.test,alpha5,alpha6){
  t0 <- Sys.time()
  cov.step <- stepwiseF(data.train,alpha5,alpha6)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  
  # Selections
  ns <- length(cov.step)
  
  # AUC
  fit <- glm(y~.,data=data.train[,c('y',cov.step)],family='binomial')
  auc.train <- performance(prediction(predict(fit),data.train$y),measure='auc')@y.values[[1]]
  auc.val <- performance(prediction(predict(fit,newdata=data.test),data.test$y),measure='auc')@y.values[[1]]
  
  return(c(ns,auc.train,auc.val,str_c(cov.step, collapse = ', '),time.step,0))
}

# R-SVS with Sieving
sieve <- function(data.train,data.test,block_size,P,alpha1,alpha2,alpha3,alpha4){
  # Sieve with permutation
  cov_sieve <- NULL
  t0 <- Sys.time()
  for (p in 1:P){
    data_pp <- permutation_partition(data.train,block_size)
    cov <- unlist(map(data_pp, stepwiseF, alpha_E=alpha1, alpha_R=alpha2))
    cov_sieve <- unique(c(cov_sieve,cov))
  }
  n_can <- length(cov_sieve)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  cov_sieve[!cov_sieve%in%colnames(data)] <- str_sub(cov_sieve[!cov_sieve%in%colnames(data)],2)
  cov_final <- stepwiseF(data.train[,c('y',cov_sieve)],alpha3,alpha4)

  # Selections
  ns <- length(cov_final)

  # AUC
  fit <- glm(y~.,data=data.train[,c('y',cov_final)],family='binomial')
  auc.train <- performance(prediction(predict(fit),data.train$y),measure='auc')@y.values[[1]]
  auc.val <- performance(prediction(predict(fit,newdata=data.test),data.test$y),measure='auc')@y.values[[1]]
  
  return(c(ns,auc.train,auc.val,str_c(cov_final, collapse = ', '),time.step,n_can))
}

# LASSO
lasso <- function(data.train,data.test){
  y.train <- data.train[,1]
  x.train <- data.train[,-1]
  y.test <- data.test[,1]
  x.test <- data.test[,-1]
  
  # LASSO
  t0 <- Sys.time()
  lasFit <- cv.glmnet(x=as.matrix(x.train), y=y.train, family = "binomial",alpha=1,nfolds = 10)
  model.min <- glmnet(x=as.matrix(x.train), y=y.train, family = "binomial", alpha = 1,lambda=lasFit$lambda.min)
  time.lasso <- as.numeric(Sys.time()-t0,units='secs')
  cov.lasso <- rownames(coef(model.min))[coef(model.min)[,1]!=0][-1]
  ncov.lasso <- length(cov.lasso)
  auc.lasso <- performance(prediction(predict(model.min,as.matrix(x.train)),y.train),measure='auc')@y.values[[1]]
  auc.val.lasso <- performance(prediction(predict(model.min,as.matrix(x.test)),y.test),measure='auc')@y.values[[1]]
  
  return (c(ncov.lasso,auc.lasso,auc.val.lasso,str_c(cov.lasso, collapse = ', '),time.lasso,0))
}

# EN
en <- function(data.train,data.test){
  y.train <- data.train[,1]
  x.train <- data.train[,-1]
  y.test <- data.test[,1]
  x.test <- data.test[,-1]
  
  # EN
  alpha <- seq(0,1,by=0.1)
  lambda <- rep(NA,length(alpha))
  cve <- rep(NA,length(alpha))
  t0 <- Sys.time()
  for (i in 1:length(alpha)){
    model.en <- cv.glmnet(x=as.matrix(x.train), y=y.train,family = "binomial",alpha=alpha[i],nfolds = 10)
    lambda[i] <- model.en$lambda.min
    cve[i] <- min(model.en$cvm)
  }
  model.min.en <- glmnet(x=as.matrix(x.train), y=y.train,family = "binomial",alpha=alpha[cve==min(cve)],lambda=lambda[cve==min(cve)])
  time.en <- as.numeric(Sys.time()-t0,units='secs')
  cov.en <- rownames(coef(model.min.en))[coef(model.min.en)[,1]!=0][-1]
  ncov.en <- length(cov.en)
  auc.en <- performance(prediction(predict(model.min.en,as.matrix(x.train)),y.train),measure='auc')@y.values[[1]]
  auc.val.en <- performance(prediction(predict(model.min.en,as.matrix(x.test)),y.test),measure='auc')@y.values[[1]]
  
  return (c(ncov.en,auc.en,auc.val.en,str_c(cov.en, collapse = ', '),time.en,0))
}

# Real data analysis binary
rda_bin <- function(data,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Hold-out 100 samples for validation
  index <- sample(1:nrow(data),size=100)
  data.test <- data[index,]
  data.train <- data[-index,]
  
  # Models
  out_step_100 <- sieve(data.train,data.test,block_size_sieve,100,alpha1,alpha2,alpha3,alpha4)
  out_step_1000 <- sieve(data.train,data.test,block_size_sieve,1000,alpha1,alpha2,alpha3,alpha4)
  out_lasso <- lasso(data.train,data.test)
  out_en <- en(data.train,data.test)
  
  # Combine results
  res <- c(out_lasso,out_en,out_step_100,out_step_1000)
  return(res)
}

# Real data analysis binary_50-50 hold-out
rda_bin_50hold <- function(data,block_size_sieve,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6){
  # 50-50 hold-out by case
  sam0 <- sample(which(data$y==0),ceiling(length(which(data$y==0))/2))
  sam1 <- sample(which(data$y==1),ceiling(length(which(data$y==1))/2))
  data.train <- data[c(sam0,sam1),]
  data.test <- data[-c(sam0,sam1),]
  
  # Models
  out_step <- l_svs(data.train,data.test,alpha5,alpha6)
  out_step_100 <- sieve(data.train,data.test,block_size_sieve,100,alpha1,alpha2,alpha3,alpha4)
  out_lasso <- lasso(data.train,data.test)
  out_en <- en(data.train,data.test)
  
  # Combine results
  res <- c(out_lasso,out_en,out_step_100,out_step)
  return(res)
}

######################
# Survival Outcome
# Permutation and partition
permutation_partition_surv <- function(data,block_size){
  y <- data[,c(1,2)]
  x <- data[,-c(1,2)]
  block <- sample(rep(1:ceiling(ncol(x)/block_size),each=block_size,length.out=ncol(x)))
  x_df <- data.frame(t(rbind(x,block)))
  x_pp <- split(x_df,x_df[,ncol(x_df)])
  data_pp <- lapply(x_pp,function(df) data.frame(cbind(y,t(df[,-ncol(df)]))))
  
  return(data_pp)
}

# Pure R-SVS
c_svs <- function(data.train,data.test,alpha5,alpha6){
  t0 <- Sys.time()
  cov.step <- stepwiseF_surv(data.train,alpha5,alpha6)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  
  # Selections
  ns <- length(cov.step)
  
  # C-index and -log10 p-value
  fit <- coxph(Surv(eventtime,status)~.,data=data.train[,c('eventtime','status',cov.step)])
  linear.step.x <- predict(fit)
  linear.step.x.val <- predict(fit,newdata=data.test)
  s.step <- summary(coxph(Surv(data.train$eventtime,data.train$status)~linear.step.x))
  s.step.val <- summary(coxph(Surv(data.test$eventtime,data.test$status)~linear.step.x.val))
  nlp.step <- -log10(coef(s.step)[5])
  nlp.val.step <- -log10(coef(s.step.val)[5])
  cin.step <- s.step$concordance[1]
  cin.val.step <- s.step.val$concordance[1]
  
  return(c(ns,nlp.step,nlp.val.step,cin.step,cin.val.step,str_c(cov.step,collapse = ', '),time.step,0))
}

# R-SVS with Sieving
sieve_surv <- function(data.train,data.test,block_size,P,alpha1,alpha2,alpha3,alpha4){
  # Sieve with permutation
  cov_sieve <- NULL
  t0 <- Sys.time()
  for (p in 1:P){
    data_pp <- permutation_partition_surv(data.train,block_size)
    cov <- unlist(map(data_pp, stepwiseF_surv, alpha_E=alpha1, alpha_R=alpha2))
    cov_sieve <- unique(c(cov_sieve,cov))
  }
  n_can <- length(cov_sieve)
  time.step <- as.numeric(Sys.time()-t0,units='secs')
  cov_sieve[!cov_sieve%in%colnames(data)] <- str_sub(cov_sieve[!cov_sieve%in%colnames(data)],2)
  cov_final <- stepwiseF_surv(data.train[,c('eventtime','status',cov_sieve)],alpha3,alpha4)

  # Selections
  ns <- length(cov_final)

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
  
  return(c(ns,nlp.step,nlp.val.step,cin.step,cin.val.step,str_c(cov_final,collapse = ', '),time.step,n_can))
}

# LASSO
lasso_surv <- function(data.train,data.test){
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
  linear.x <- predict(model.min,as.matrix(x.train))
  linear.x.val <- predict(model.min,as.matrix(x.test))
  s <- summary(coxph(Surv(y.train$eventtime,y.train$status)~linear.x))
  s.val <- summary(coxph(Surv(y.test$eventtime,y.test$status)~linear.x.val))
  nlp.lasso <- -log10(coef(s)[5])
  nlp.val.lasso <- -log10(coef(s.val)[5])
  cin.lasso <- s$concordance[1]
  cin.val.lasso <- s.val$concordance[1]
  
  return (c(ncov.lasso,nlp.lasso,nlp.val.lasso,cin.lasso,cin.val.lasso,str_c(cov.lasso,collapse = ', '),time.lasso,0))
}

# EN
en_surv <- function(data.train,data.test){
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
  linear.en.x <- predict(model.min.en,as.matrix(x.train))
  linear.en.x.val <- predict(model.min.en,as.matrix(x.test))
  s.en <- summary(coxph(Surv(y.train$eventtime,y.train$status)~linear.en.x))
  s.en.val <- summary(coxph(Surv(y.test$eventtime,y.test$status)~linear.en.x.val))
  nlp.en <- -log10(coef(s.en)[5])
  nlp.val.en <- -log10(coef(s.en.val)[5])
  cin.en <- s.en$concordance[1]
  cin.val.en <- s.en.val$concordance[1]
  
  return (c(ncov.en,nlp.en,nlp.val.en,cin.en,cin.val.en,str_c(cov.en,collapse = ', '),time.en,0))
}

# Real data analysis survival
rda_surv <- function(data,block_size_sieve,alpha1,alpha2,alpha3,alpha4){
  # Hold-out 100 samples for validation
  index <- sample(1:nrow(data),size=100)
  data.test <- data[index,]
  data.train <- data[-index,]
  
  # Models
  out_step_100 <- sieve_surv(data.train,data.test,block_size_sieve,100,alpha1,alpha2,alpha3,alpha4)
  out_step_1000 <- sieve_surv(data.train,data.test,block_size_sieve,1000,alpha1,alpha2,alpha3,alpha4)
  out_lasso <- lasso_surv(data.train,data.test)
  out_en <- en_surv(data.train,data.test)
  
  # Combine results
  res <- c(out_lasso,out_en,out_step_100,out_step_1000)
  return(res)
}

# Real data analysis survival_50-50 hold-out
rda_surv_50hold <- function(data,block_size_sieve,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6){
  # 50-50 hold-out by case
  sam0 <- sample(which(data$status==0),ceiling(length(which(data$status==0))/2))
  sam1 <- sample(which(data$status==1),ceiling(length(which(data$status==1))/2))
  data.train <- data[c(sam0,sam1),]
  data.test <- data[-c(sam0,sam1),]
  
  # Models
  out_step <- c_svs(data.train,data.test,alpha5,alpha6)
  out_step_100 <- sieve_surv(data.train,data.test,block_size_sieve,100,alpha1,alpha2,alpha3,alpha4)
  out_lasso <- lasso_surv(data.train,data.test)
  out_en <- en_surv(data.train,data.test)
  
  # Combine results
  res <- c(out_lasso,out_en,out_step_100,out_step)
  return(res)
}
