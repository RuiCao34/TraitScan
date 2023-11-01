#From Multi_trait_Scan_individual_v1.R
#extend to summary level data

# Function to simulate summary level data
# b1 is effect size, 
# d is the number of trait asspciated with the SNP
pleio.data <- function(n,p,d,maf,rho,b1){
  # library(MASS)
  #randomly select d trait to be associate with the SNP
  if (d > 0) beta <- sample(c(rep(b1,d), rep(0, p-d)), p, replace=F)
  if (d == 0) beta <- rep(0,p)
  truebeta=(beta>0)
  # Simulate dosage of minor allele
  x <- rbinom(n, 2, maf) # 1 by n
  
  # Covariance matrix
  Sigma <- matrix(rho, nrow=p, ncol=p) # p by p
  diag(Sigma) <- 1
  
  
  # Simulate traits for each individual
  epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    epsilon[i,] <- MASS::mvrnorm(1, rep(0,p), Sigma)
    y[i,] <- x[i]*beta + epsilon[i,]
  }
  
  est.epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  z <- rep(NA, p)    
  beta <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    #standardized x and y
    x <- scale(x, center=T, scale=T)
    y[,j] <- scale(y[,j], center=T, scale=T)
    model <- lm(y[,j] ~ x)
    est.epsilon[,j] <- model$resid
    z[j] <- summary(model)$coefficient[2,3]
    # z[j] <- qnorm(pt(summary(model)$coefficient[2,3],n-2))
    beta[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  
  sigma <- cov(est.epsilon)
  output <- cbind(beta,se,z,n,p)
  return(list("y"=y,"x"=x, "sigma"=sigma,"summary"=output,"truebeta"=truebeta))
}

pleio.data.het <- function(n,p1,p2,d1,d2,maf,rho,b1,b2){
  # library(MASS)
  p = p1 + p2
  d = d1 + d2
  #randomly select d trait to be associate with the SNP
  if (d > 0) beta <- sample(c(rep(b1,d1), rep(b2,d2),rep(0, p-d)), p, replace=F)
  if (d == 0) beta <- rep(0,p)
  truebeta=(!beta==0)
  # Simulate dosage of minor allele
  x <- rbinom(n, 2, maf) # 1 by n
  
  # Covariance matrix
  Sigma <- matrix(rho, nrow=p, ncol=p) # p by p
  diag(Sigma) <- 1
  
  
  # Simulate traits for each individual
  epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    epsilon[i,] <- MASS::mvrnorm(1, rep(0,p), Sigma)
    y[i,] <- x[i]*beta + epsilon[i,]
  }
  
  est.epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  z <- rep(NA, p)    
  beta <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    #standardized x and y
    x <- scale(x, center=T, scale=T)
    y[,j] <- scale(y[,j], center=T, scale=T)
    model <- lm(y[,j] ~ x)
    est.epsilon[,j] <- model$resid
    z[j] <- summary(model)$coefficient[2,3]
    # z[j] <- qnorm(pt(summary(model)$coefficient[2,3],n-2))
    beta[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  
  sigma <- cov(est.epsilon)
  output <- cbind(beta,se,z,n,p)
  return(list("y"=y,"x"=x, "sigma"=sigma,"summary"=output,"truebeta"=truebeta))
}


pleio.data.N <- function(n1,n2,p1,p2,d1,d2,maf,rho,b1,n_effect){
  # default n1>n2
  n <- n1
  p <- p1+p2
  d <- d1+d2
  
  # library(MASS)
  #randomly select d trait to be associate with the SNP
  
  if (d1 > 0) beta_part1 <- sample(c(rep(b1,d1), rep(0, p1-d1)), p1, replace=F)
  if (d1 == 0) beta_part1 <- rep(0,p1)
  if (d2 > 0) beta_part2 <- sample(c(rep(b1,d2), rep(0, p2-d2)), p2, replace=F)
  if (d2 == 0) beta_part2 <- rep(0,p2)
  
  beta=c(beta_part1, beta_part2)
  truebeta=c(beta_part1>0,beta_part2>0)
  
  # Simulate dosage of minor allele
  x <- rbinom(n1, 2, maf) # 1 by n
  
  # Covariance matrix
  Sigma <- matrix(rho, nrow=p, ncol=p) # p by p
  diag(Sigma) <- 1
  
  
  # Simulate traits for each individual
  epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    epsilon[i,] <- MASS::mvrnorm(1, rep(0,p), Sigma)
    y[i,] <- x[i]*beta + epsilon[i,]
  }
  
  est.epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  z <- rep(NA, p)    
  beta <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    #standardized x and y
    x <- scale(x, center=T, scale=T)
    y[,j] <- scale(y[,j], center=T, scale=T)
    if(j <= p1){
      model <- lm(y[,j] ~ x)
      z[j] <- qnorm(pt(summary(model)$coefficient[2,3],n1-2))
    }else{
      model <- lm(y[1:n2,j] ~ x[1:n2])
      z[j] <- qnorm(pt(summary(model)$coefficient[2,3],n2-2))
    }
    # est.epsilon[,j] <- model$resid
    beta[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  
  # sigma <- cov(est.epsilon)
  n = n_effect
  output <- cbind(beta,se,z,n,p)
  return(list("y"=y,"x"=x, "summary"=output,"truebeta"=truebeta))
}



pleio.data.binary <- function(n,p,d,maf,rho,b1,b0){
  # library(MASS)
  # library(bindata)
  #randomly select d trait to be associate with the SNP
  if (d > 0) beta1 <- sample(c(rep(b1,d), rep(0, p-d)), p, replace=F)
  if (d == 0) beta1 <- rep(0,p)
  truebeta=(beta1>0)
  # Simulate dosage of minor allele
  x <- rbinom(n, 2, maf) # 1 by n
  pai=exp(b0+x*beta1)/(1+exp(b0+x*beta1))
  
  # Covariance matrix
  Sigma <- matrix(rho, nrow=p, ncol=p) # p by p
  diag(Sigma) <- 1
  
  # Simulate traits for each individual
  epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    pai=exp(b0+x[i]*beta1)/(1+exp(b0+x[i]*beta1))
    y[i,] <- bindata::rmvbin(1, pai, bincorr=Sigma)
  }
  
  
  z <- rep(NA, p)    
  beta <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    model <- glm(y[,j] ~ x, family="binomial")
    z[j] <- summary(model)$coefficient[2,3]
    beta[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  
  
  output <- cbind(beta,se,z,n,p)
  return(list("y"=y,"x"=x,"summary"=output,"truebeta"=truebeta))
}


pleio.data.binary.v2 <- function(n,p,d,maf,rho,b1,threshold){
  # library(MASS)
  #randomly select d trait to be associate with the SNP
  if (d > 0) beta <- sample(c(rep(b1,d), rep(0, p-d)), p, replace=F)
  if (d == 0) beta <- rep(0,p)
  truebeta=(beta>0)
  # Simulate dosage of minor allele
  x <- rbinom(n, 2, maf) # 1 by n
  
  # Covariance matrix
  Sigma <- matrix(rho, nrow=p, ncol=p) # p by p
  diag(Sigma) <- 1
  
  
  # Simulate traits for each individual
  epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    epsilon[i,] <- MASS::mvrnorm(1, rep(0,p), Sigma)
    y[i,] <- ifelse(x[i]*beta + epsilon[i,]>threshold,1,0)
  }
  
  est.epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  z <- rep(NA, p)    
  beta <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    #standardized x and y
    x <- scale(x, center=T, scale=T)
    model <- glm(y[,j] ~ x, family="binomial")
    z[j] <- summary(model)$coefficient[2,3]
    beta[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  
  output <- cbind(beta,se,z,n,p)
  return(list("y"=y,"x"=x,"summary"=output,"truebeta"=truebeta))
}


pleio.data.sim <- function(n,beta_cont, beta_bin,maf,variance,threshold){
  # library(MASS)
  #randomly select d trait to be associate with the SNP
  p_cont <- length(beta_cont)
  p_bin <- length(beta_bin)
  p <- p_cont + p_bin
  # if(length(p_cont)==0){
  #   beta <- sample(beta_bin, length(beta_bin), replace=F)
  # }else if(length(p_bin)==0){
  #   beta <- sample(beta_cont, length(beta_cont), replace=F)
  # }else{
  #   beta <- c(sample(beta_cont, length(beta_cont), replace=F), sample(beta_bin, length(beta_bin), replace=F))
  # }
  beta <- c(beta_cont, beta_bin)
  
  truebeta=!(beta==0)
  # Simulate dosage of minor allele
  x <- rbinom(n, 2, maf) # 1 by n
  
  # Covariance matrix
  Sigma <- variance
  
  # Simulate traits for each individual
  epsilon <- MASS::mvrnorm(n, rep(0,p), Sigma) # n by p
  y <- matrix(NA, nrow=n, ncol=p) # n by p
  for (i in 1:n){
    y[i,] <- x[i]*beta + epsilon[i,]
    if(p_bin>0){
      y[i,(p-p_bin+1):p] <- ifelse(y[i,(p-p_bin+1):p]>threshold,1,0)
    }
  }
  
  est.epsilon <- matrix(NA, nrow=n, ncol=p) # n by p
  z <- rep(NA, p)    
  beta_hat <- rep(NA, p)  
  se <- rep(NA,p)
  for (j in 1:p){
    #standardized x and y
    x <- scale(x, center=T, scale=T)
    if(j <= p_cont){
      model <- lm(y[,j] ~ x)
    }else{
      model <- glm(y[,j] ~ x, family="binomial")
    }
    z[j] <- summary(model)$coefficient[2,3]
    beta_hat[j] <- summary(model)$coefficient[2,1]
    se[j] <- summary(model)$coefficient[2,2]
  }
  output <- cbind(beta_hat,se,z,n,p)
  return(list("y"=y,"x"=x,"summary"=output,"truebeta"=truebeta))
}

# get a summary statistics matrix for q SNPs, p traits, and 5 summary statistics (beta, se_beta, et al)
# from the second variant, it is all null variants (for all the p-1 traits) and will be used for the purpose of estimating the var-cov matrix (simulation purpose)
# maf is the maf of the truly associated snp, other SNPs have their maf drawed randomly

#seperate null and interest snp
interest <- function(n, p, d, maf, rho, b1){
  temp <- pleio.data(n=n, p, d, maf=maf, rho, b1=b1)	
  return(temp)
}

null <- function(n,p,rho,q,maf){
  summary_data <- array(dim=c(q+1,p,5))
  maf_all <- rep(maf,q+1)
  for (i in 1:q+1){
    temp <- pleio.data(n=n, p, d=0, maf=maf_all[i], rho, b1=0)
    summary_data[i,,] <- temp$summary
  }
  Zscore.null <- summary_data[1:q+1,,3]
  se.null <- summary_data[1:q+1,,2]
  beta.null <- summary_data[1:q+1,,1]
  comp1 <- cor(Zscore.null) * n
  # comp1 <- cor(Zscore.null) * (n-p)
  return(comp1)
}

null.binary <- function(n,p,rho,q,maf,b0){
  summary_data <- array(dim=c(q+1,p,5))
  maf_all <- rep(maf,q+1)
  for (i in 1:q+1){
    temp <- pleio.data.binary(n=n, p, d=0, maf=maf_all[i], rho, b1=0, b0=b0)
    summary_data[i,,] <- temp$summary
  }
  Zscore.null <- summary_data[1:q+1,,3]
  se.null <- summary_data[1:q+1,,2]
  beta.null <- summary_data[1:q+1,,1]
  comp1 <- cor(Zscore.null) * n
  return(comp1)
}

null.binary.v2 <- function(n,p,rho,q,maf,threshold){
  summary_data <- array(dim=c(q+1,p,5))
  maf_all <- rep(maf,q+1)
  for (i in 1:q+1){
    temp <- pleio.data.binary.v2(n=n, p, d=0, maf=maf_all[i], rho, b1=0,threshold=threshold)
    summary_data[i,,] <- temp$summary
  }
  Zscore.null <- summary_data[1:q+1,,3]
  se.null <- summary_data[1:q+1,,2]
  beta.null <- summary_data[1:q+1,,1]
  comp1 <- cor(Zscore.null) * n
  return(comp1)
}

null.binary.v3 <- function(n,p,rho,q,maf,threshold){
  summary_data <- array(dim=c(q+1,p,5))
  maf_all <- rep(maf,q+1)
  for (i in 1:q+1){
    temp <- pleio.data.binary.v2(n=n, p, d=0, maf=maf_all[i], rho, b1=0,threshold=threshold)
    summary_data[i,,] <- temp$summary
  }
  Zscore.null <- summary_data[1:q+1,,3]
  se.null <- summary_data[1:q+1,,2]
  beta.null <- summary_data[1:q+1,,1]
  return(cor(Zscore.null))
}

null.mixed <- function(n,p_cont, p_bin,maf,variance,threshold,q){
  p <- p_cont + p_bin
  summary_data <- array(dim=c(q,p,5))
  maf_all <- rep(maf,q)
  for (i in 1:q){
    summary_data[i,,] <- pleio.data.sim(n,beta_cont=rep(0,p_cont),beta_bin=rep(0,p_bin),maf,variance,threshold)$summary
  }
  Zscore.null <- summary_data[1:(q),,3]
  se.null <- summary_data[1:(q),,2]
  beta.null <- summary_data[1:(q),,1]
  return(cor(Zscore.null))
}

null.mixed.truncated <- function(n,p_cont, p_bin,maf,variance,threshold,q){
  # library(qlcMatrix)
  p <- p_cont + p_bin
  summary_data <- array(dim=c(q,p,5))
  maf_all <- rep(maf,q)
  for (i in 1:q){
    summary_data[i,,] <- pleio.data.sim(n,beta_cont=rep(0,p_cont),beta_bin=rep(0,p_bin),maf,variance,threshold)$summary
  }
  Zscore.null <- summary_data[1:(q),,3]
  Zscore.null[Zscore.null > -qnorm(0.025) | Zscore.null < qnorm(0.025)] = NA
  # Zscore.null <- Zscore.null[rowMax(Zscore.null) < -qnorm(0.025) & rowMin(Zscore.null) > qnorm(0.025),]
  # se.null <- summary_data[1:q+1,,2]
  # beta.null <- summary_data[1:q+1,,1]
  return(cor(Zscore.null,use="pairwise.complete.obs"))
}

# null.mixed.truncated2nontrunctated <- function(n,p_cont,p_bin,rho,q,maf,threshold){
#   library(qlcMatrix)
#   library(tmvtnorm)
#   p <- p_cont + p_bin
#   summary_data <- array(dim=c(q+1,p,5))
#   maf_all <- rep(maf,q+1)
#   for (i in 1:q+1){
#     data_cont <- pleio.data(n=n, p_cont, d=0, maf=maf_all[i], rho, b1=0)
#     data_bin <- pleio.data.binary.v2(n=n, p_bin, d=0, maf=maf_all[i], rho, b1=0,threshold=threshold)
#     summary_data[i,,] <- rbind(data_cont$summary, data_bin$summary)
#   }
#   Zscore.null <- summary_data[1:q+1,,3]
#   Zscore.null <- Zscore.null[rowMax(Zscore.null) < -qnorm(0.025) & rowMin(Zscore.null) > qnorm(0.025),]
#   # se.null <- summary_data[1:q+1,,2]
#   # beta.null <- summary_data[1:q+1,,1]
#   mle.fit1 <- mle.tmvnorm(Zscore.null, lower=rep(qnorm(0.025),p), upper=rep(-qnorm(0.025),p))
#   cor.matrix <- matrix(rep(NA,p*p), nrow=p)
#   cor.matrix[lower.tri(cor.matrix, diag = T)] <- mle.fit1@coef[(p+1):length(mle.fit1@coef)]
#   cor.matrix <- forceSymmetric(cor.matrix, uplo="L")
#   cor.matrix <- as.matrix(cor.matrix)
#   return(cor.matrix)
# }

null.mixed.assosnp <- function(n,p_cont,p_bin,rho,q,maf,threshold,snp.percen,trait.percen,b1){
  # snp.percen: percentage of snps truly associated with traits
  # trait.percen: percentage of traits truly associated with the associated snps
  p <- p_cont + p_bin
  summary_data <- array(dim=c(q+1,p,5))
  maf_all <- rep(maf,q+1)
  for (i in 1:q+1){
    data_cont <- pleio.data(n=n, p_cont, d=0, maf=maf_all[i], rho, b1=0)
    summary_data[i,,] <- data_cont$summary
    # data_bin <- pleio.data.binary.v2(n=n, p_bin, d=0, maf=maf_all[i], rho, b1=0,threshold=threshold)
    # summary_data[i,,] <- rbind(data_cont$summary, data_bin$summary)
    
  }
  
  summary_data_2 <- array(dim=c(round(q*snp.percen),p,5))
  maf_all_2 <- rep(maf,round(q*snp.percen))
  for (i in 1:round(q*snp.percen)){
    data_cont <- pleio.data(n=n, p_cont, d=round(p_cont*trait.percen), maf=maf_all[i], rho, b1=b1)
    summary_data_2[i,,] <- data_cont$summary
    # data_bin <- pleio.data.binary.v2(n=n, p_bin, d=round(p_bin*trait.percen), maf=maf_all[i], rho, b1=b1,threshold=threshold)
    # summary_data_2[i,,] <- rbind(data_cont$summary, data_bin$summary)
  }
  
  Zscore.null <- summary_data[1:q+1,,3]
  se.null <- summary_data[1:q+1,,2]
  beta.null <- summary_data[1:q+1,,1]
  return(cor(Zscore.null))
}



# null_v2 <- function(n,p,rho,q,maf,output){
#   summary_data <- array(dim=c(q+1,p,5))
#   maf_all <- rep(maf,q+1)
#   for (i in 1:q+1){
#     temp <- pleio.data(n=n, p, d=0, maf=maf_all[i], rho, b1=0)
#     summary_data[i,,] <- temp$summary
#   }
#   Zscore.null <- summary_data[1:q+1,,3]
#   se.null <- summary_data[1:q+1,,2]
#   beta.null <- summary_data[1:q+1,,1]
#   cor_z <- cor(Zscore.null)
#   
#   x <- rbinom(n, 2, maf)
#   x=x-mean(x)
#   # x=x/sqrt(sum(x^2))
#   Y_sqaure=sum(x^2)*(output[,1])^2+(n-1)*sum(x^2)*(output[,2])^2
#   comp1=diag(sqrt(Y_sqaure)) %*% cor_z %*% diag(sqrt(Y_sqaure))
#   
#   return(comp1)
# }

#function to estimate the cov(epislon) from summary level data; 
#here by default X and Y are standardized
#comp1 and 2, and 3 correspond to derivation for fun, Page 18

sigma_fun <- function(interest,null.cov,q){
  #assume sample size are equal
  
  #interest SNP
  n <- interest$summary[1,4]
  beta <- interest$summary[,1]
  se <- interest$summary[,2]
  p <- interest$summary[1,5]
  
  comp32 <- -beta %*% t(beta) * n 
  comp1 <- null.cov
  
  sigma_sum <- (comp1 + comp32)/(n-p)
  return(sigma_sum)
}

#sigma.summary <- sigma_fun(sdata,q=100)

fgss.pleio <- function(interest, null.cov, q, alpha_grid=10, test=c('BJ','HC')){
  # library(whitening)
  #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  w <- whitening::whiteningMatrix(epsilon, method='ZCA-cor')
  
  z.decor <- z %*% t(w) 
  
  pval.de <- (1-pnorm(abs(z.decor)))*2
  
  # Perform pleio.scan() with whitened after removing j'
  if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC(pval.de,alpha_grid)
  
  return(obj)
}

fgss.pleio.binary <- function(interest, null.cov, q, alpha_grid=10, test=c('BJ','HC')){
  # #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  # epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  w <- whitening::whiteningMatrix(null.cov, method='ZCA-cor')
  
  z.decor <- z %*% t(w) 
  
  pval.de <- (1-pnorm(abs(z.decor)))*2
  
  # Perform pleio.scan() with whitened after removing j'
  if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC(pval.de,alpha_grid)
  
  return(obj)
}

fgss.pleio.binary.whitened <- function(interest, w, q, alpha_grid=10, test=c('BJ','HC')){
  
  # #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  # epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  # w <- whiteningMatrix(null.cov, method='ZCA-cor')
  
  z.decor <- z %*% t(w) 
  
  pval.de <- (1-pnorm(abs(z.decor)))*2
  
  # Perform pleio.scan() with whitened after removing j'
  if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC(pval.de,alpha_grid)
  
  return(obj)
}


fgss.pleio.alphainf.binary <- function(interest, null.cov, min.alpha, max.alpha, test=c('BJ','HC')){
  # library(whitening)
  # #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  # epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  w <- whitening::whiteningMatrix(null.cov, method='ZCA-cor')
  
  z.decor <- z %*% t(w) 
  
  pval.de <- (1-pnorm(abs(z.decor)))*2
  
  # Perform pleio.scan() with whitened after removing j'
  # if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC.alphainf(pval.de, min.alpha, max.alpha)
  
  return(obj)
}

pleio.scan.BJ <- function(pval, alpha_grid){

  # alpha <- seq(0, 1, length.out=alpha_grid)
  alpha <- alpha_grid
  todo <- sort(pval) # order of pvalues, small to large
  s <- length(todo)

  # Define BJ function
  BJ <- function(x, y){x * log(x/y) + (1-x)*log((1-x)/(1-y)) }

  # Run subset scan
  max.f <- 0
  max.subset <- NULL
  max.alpha <- 0
  l <- matrix(NA, nrow=length(alpha_grid), ncol=s)
  # Scan over a grid of {s} and alpha values
  if(min(todo) > max(alpha)){
    max.f <- -Inf
    max.alpha <- integer(0)
    #max.subset <- which(pval < max.alpha)
    max.size <- 0
    max.subset <- integer(0)
  }else{
    for (ss in seq(1,s)){
      for (aa in 1:length(alpha_grid)){
        a <- alpha[aa]
        l[aa,ss] <- ss * BJ(sum(todo[seq(1,ss)] < a)/ss, a)
      }
    }
    l[is.nan(l)] = -Inf
    
    # Calculate max F
    max.f <- max(l,na.rm=T)
    temp <- which(l==max.f, arr.ind=T)
    max.alpha <- alpha[temp[1]]
    #max.subset <- which(pval < max.alpha)
    max.size <- temp[2]
    max.subset <- which(pval %in% todo[1:temp[2]])
  }
  
  

  return(list("subset"=max.size, "statistic"=max.f, "max.alpha"=max.alpha, "select"=max.subset))
}

pleio.scan.HC <- function(pval, alpha_grid){
  
  #alpha <- seq(0, 1, length.out=alpha_grid) 
  # alpha <- alpha_grid
  
  
  todo <- sort(pval) # order of pvalues, small to large
  s <- length(todo)
  
  # Run subset scan
  max.f <- 0
  max.subset <- NULL
  max.alpha <- 0
  l <- matrix(NA, nrow=length(alpha_grid), ncol=s)
  # Scan over a grid of {s} and alpha values
  for (ss in seq(1,s)){
    for (aa in 1:length(alpha_grid)){ 
      a <- alpha_grid[aa]	  	   
      l[aa,ss] <- (sum(todo[seq(1,ss)] < a) - ss * a)/sqrt(ss * a * (1-a))
    }
  }
  
  # Calculate max F
  max.f <- max(l,na.rm=T)
  temp <- which(l==max.f, arr.ind=T)
  max.alpha <- alpha_grid[temp[1]]
  #max.subset <- which(pval < max.alpha)
  max.size <- temp[2]
  max.subset <- which(pval %in% todo[1:temp[2]])
  
  return(list("subset"=max.size, "statistic"=max.f, "max.alpha"=max.alpha, "select"=max.subset))
}


pleio.scan.HC.alphainf <- function(pval, min.alpha, max.alpha){
  
  #alpha <- seq(0, 1, length.out=alpha_grid) 
  # alpha <- alpha_grid
  pval[pval < min.alpha] = min.alpha
  pval[pval > max.alpha] = max.alpha
  
  todo <- sort(pval) # order of pvalues, small to large
  s <- length(todo)
  
  # Run subset scan
  max.f <- 0
  max.subset <- NULL
  max.alpha <- 0
  l <- unlist(lapply(1:s, function(i){sqrt(i * (1 - todo[i]) / todo[i])}))
  
  
  # Calculate max F
  max.f <- max(l,na.rm=T)
  temp <- which(l==max.f)
  # max.alpha <- alpha_grid[temp[1]]
  #max.subset <- which(pval < max.alpha)
  max.size <- temp
  max.subset <- which(pval %in% todo[1:temp])
  
  return(list("subset"=max.size, "statistic"=max.f, "select"=max.subset))
}



# Function to simulate draws from the null distribution of F(S), using Monte Carlo simulations
null.stat <- function(p, alpha_grid, test){
  z.m <- rnorm(p,0,1)
  pval.de <- (1-pnorm(abs(z.m)))*2
  if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC(pval.de,alpha_grid)
  #obj <- pleio.scan(pval.de, alpha_grid)
  # Return the statistic
  return(obj$statistic)
}

null.stat.alphainf <- function(p, min.alpha, max.alpha, test){
  z.m <- rnorm(p,0,1)
  pval.de <- (1-pnorm(abs(z.m)))*2
  # if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC.alphainf(pval.de, min.alpha, max.alpha)
  #obj <- pleio.scan(pval.de, alpha_grid)
  # Return the statistic
  return(obj$statistic)
}

null.stat.mclapply <- function(sim, p, alpha_grid, test){
  z.m <- rnorm(p,0,1)
  pval.de <- (1-pnorm(abs(z.m)))*2
  if (test=='BJ') obj <- pleio.scan.BJ(pval.de,alpha_grid)
  if (test=='HC') obj <- pleio.scan.HC(pval.de,alpha_grid)
  #obj <- pleio.scan(pval.de, alpha_grid)
  # Return the statistic
  return(obj$statistic)
}

null.stat.theoretical <- function(n.p.values, alpha_grid, test.statistic){
  m = length(alpha_grid)
  threshold = alpha_grid/(1-alpha_grid)*test.statistic^2
  threshold = ceiling(threshold) - 1
  threshold[threshold > n.p.values] = n.p.values
  conditional.prob = vector("list", length = m)
  # temp = lapply(0:threshold[m], function(x){choose(n.p.values,x) * (alpha_grid[m])^(x) * (1-alpha_grid[m])^(n.p.values-x)})
  temp = lapply(0:threshold[m], function(x){dbinom(x,n.p.values,alpha_grid[m])})
  conditional.prob[[m]] = unlist(temp)
  for(i in (m-1):1){
    temp.matrix = matrix(0, nrow=1+threshold[i], ncol=length(conditional.prob[[i+1]]))
    ## column: (value of m-1)-1; row: given (value in m)-1
    for(j in 1:ncol(temp.matrix)){
      # temp = lapply(0:min(threshold[i],j-1), function(x){choose(j-1,x) * (alpha_grid[i]/alpha_grid[i+1])^(x) * (1-alpha_grid[i]/alpha_grid[i+1])^(j-1-x)})
      temp = lapply(0:min(threshold[i],j-1), function(x){dbinom(x,j-1,alpha_grid[i]/alpha_grid[i+1])})
      temp.matrix[1:min(threshold[i]+1,j),j] = unlist(temp) * conditional.prob[[i+1]][j]
    }
    
    ## joint prob/marginal prob
    if(sum(conditional.prob[[i+1]]) ==0 ){
      conditional.prob[1:i] = 0
      break
    }else{
      conditional.prob[[i]] = as.vector(rowSums(temp.matrix))/sum(conditional.prob[[i+1]])
    }
  }
  return(1 - prod(unlist(lapply(conditional.prob, sum))))
}


# null.stat.alphainf.theoretical <- function(n.p.values, test.statistic, min.alpha, max.alpha){
#   m = n.p.values
#   alpha_grid = (1:n.p.values) / ((1:n.p.values) + test.statistic^2)
#   threshold = (1:m)-1
#   # threshold[alpha_grid < min.alpha] = threshold[min(which(alpha_grid >= min.alpha))]
#   # threshold[alpha_grid > max.alpha] = threshold[max(which(alpha_grid <= max.alpha))]
#   # alpha_grid[alpha_grid < min.alpha] = min.alpha
#   # alpha_grid[alpha_grid > max.alpha] = max.alpha
#   
#   threshold <- threshold[alpha_grid > min.alpha & alpha_grid < max.alpha]
#   alpha_grid <- alpha_grid[alpha_grid > min.alpha & alpha_grid < max.alpha]
#   m = length(threshold)
#   
#   
#   conditional.prob = vector("list", length = m)
#   temp = lapply(0:threshold[m], function(x){choose(n.p.values,x) * (alpha_grid[m])^(x) * (1-alpha_grid[m])^(n.p.values-x)})
#   conditional.prob[[m]] = unlist(temp)
#   if(m > 1){
#     for(i in (m-1):1){
#       temp.matrix = matrix(0, nrow=1+threshold[i], ncol=length(conditional.prob[[i+1]]))
#       ## column: (value of m-1)-1; row: given (value in m)-1
#       for(j in 1:ncol(temp.matrix)){
#         temp = lapply(0:min(threshold[i],j-1), function(x){choose(j-1,x) * (alpha_grid[i]/alpha_grid[i+1])^(x) * (1-alpha_grid[i]/alpha_grid[i+1])^(j-1-x)})
#         temp.matrix[1:min(threshold[i]+1,j),j] = unlist(temp) * conditional.prob[[i+1]][j]
#       }
#       
#       ## joint prob/marginal prob
#       conditional.prob[[i]] = as.vector(rowSums(temp.matrix))/sum(conditional.prob[[i+1]])
#     }
#   }
#   return(1 - prod(unlist(lapply(conditional.prob, sum))))
# }



## return G value given decor z, alpha and gamma
truncated.G <- function(z.scores.decor, gamma, power =2){
  ## alpha is a number
  z.scores.decor <- abs(z.scores.decor)
  if(max(abs(z.scores.decor)) <= gamma){
    test.statistics = 0
    select = c()
    return(result=list(test.statistics = test.statistics, select = select))
  }else{
    test.statistics = sum((z.scores.decor[z.scores.decor > gamma])^power)
    select = which(z.scores.decor > gamma)
    return(result=list(test.statistics = test.statistics, select = select))
  }
}

## null dist given alpha and gamma
null.dist.truncated.G <- function(N, p, alpha=0, gamma){
  vector = rep(0, N)
  for(i in 1:N){
    z.m <- rnorm(p,0,1)
    obj <- truncated.G(z.m, alpha, gamma)
    vector[i] <- obj$test.statistics
  }
  return(vector)
}

truncated.power.test <- function(interest, null.cov, cutoff, null.dist, gamma){
  ## cutoff is vector of alpha
  ## null.dist is a list of null dist, length = length(threshold)
  # library(expm)
  # z.scores <- interest$summary[,3]
  # whiten.matrix <- sqrtm(solve(null.cov))
  # z.scores.decor <- t(z.scores) %*% whiten.matrix
  
  # #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  # epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  w <- whitening::whiteningMatrix(null.cov, method='ZCA-cor')
  
  z.scores.decor <- z %*% t(w) 
  
  
  
  # threshold <- lapply(abs(z.scores.decor)[abs(z.scores.decor) > min(threshold)], 
  #                     function(z){max(threshold[threshold < z])})
  # threshold <- unlist(threshold)
  # threshold <- threshold[!duplicated(threshold)]
  
  m = length(cutoff)
  G.statistic <- vector("list", length = m)
  H.statistic <- rep(0, m)
  for(i in 1:m){
    G.statistic[[i]] <- truncated.G(z.scores.decor, cutoff[i], gamma)
    H.statistic[i] <- sum(null.dist[[i]] >= G.statistic[[i]]$test.statistics)/length(null.dist[[i]])
  }
  
  select.alpha = which(H.statistic == min(H.statistic))
  select.alpha = select.alpha[1]
  test.statistics = H.statistic[[select.alpha]]
  select = G.statistic[[select.alpha]]$select
  return(result=list(test.statistics = test.statistics, select = select))
  
}


null.dist.power.test <- function(N, p, cutoff, null.dist, gamma){
  vector = rep(0, N)
  for(i in 1:N){
    z.m <- rnorm(p,0,1)
    obj <- truncated.power.test(interest=list(summary=cbind(z.m,z.m,z.m)), 
                                null.cov = diag(rep(1,p)), 
                                cutoff, 
                                null.dist, 
                                gamma)
    vector[i] <- obj$test.statistics
  }
  return(vector)
}

# gamma.test.theor.p <- function(gamma, p, test.statistic){
#   mu_i_gamma_even <- function(gamma){factorial(gamma)/(2^(gamma/2))*
#   do.call("sum",lapply(0:(gamma/2), function(d) 1/factorial(d)/factorial(gamma/2-d)))/
#   (2^(gamma/2))}
#   
#   
#   if(gamma == 1){
#     mu_gamma <- NA
#     sigma2_gamma <- NA
#   }else if(gamma %% 2 ==0){
#     ## let n1=n2=1, sigma_ii=0.5
#     mu_i_gamma <- mu_i_gamma_even(gamma)
#     
#     mu_gamma <- mu_i_gamma * p
#     sigma2_gamma <- mu_i_gamma_even(gamma*2)*p - (mu_i_gamma_even(gamma)^2)*p
#   }else{
#     mu_gamma <- NA
#     sigma2_gamma <- NA
#   }
#   return(2*(1-pnorm((test.statistic-mu_gamma)/sqrt(sigma2_gamma))))
# 
# }


combined.test <- function(interest, null.cov, null.dist, alpha_grid, null.dist.power.test.gamma2){
  ## threshold is vector of alpha
  ## null.dist is a list of null dist, length = length(threshold)
  cutoff = sapply(alpha_grid, function(x){
    abs(qnorm(x/2, 0, 1))
  })
  
  obj_hc <- fgss.pleio.binary(interest, null.cov, q=0, alpha_grid=alpha_grid, test='HC')
  pvalue_hc <- null.stat.theoretical(dim(null.cov)[1], alpha_grid=alpha_grid, obj_hc$statistic)
  
  obj_gamma2 = truncated.power.test(interest, null.cov, cutoff, null.dist, gamma=2)
  pvalue_gamma2 <- sum(obj_gamma2$test.statistics >= null.dist.power.test.gamma2)/length(null.dist.power.test.gamma2)
  
  test.statistics = min(pvalue_hc,pvalue_gamma2)
  if(pvalue_hc < pvalue_gamma2){
    select <- obj_hc$select
    minp.test = "HC"
  }else{
    select <- obj_gamma2$select
    minp.test = "TC"
    }
  
  return(result=list(test.statistics = test.statistics, select = select, 
                     select.HC =obj_hc$select, select.TC = obj_gamma2$select, minp.test = minp.test))
  
}

null.dist.combined.test <- function(N, p, alpha_grid, null.dist, null.dist.power.test.gamma2){
  vector = rep(0, N)
  for(i in 1:N){
    z.m <- rnorm(p,0,1)
    obj <- combined.test(interest=list(summary=cbind(z.m,z.m,z.m)),
                         null.cov = diag(rep(1,p)),
                         null.dist, alpha_grid, null.dist.power.test.gamma2)
    vector[i] <- obj$test.statistics
  }
  return(vector)
}

TC.null.dist.asy <- function(p, cutoff, nsim){
  m = length(cutoff)
  para.asy <- vector("list", length = m)
  z.matrix = matrix(rep(0, nsim*p), nrow=nsim)
  for(i in 1:nsim){
    z.matrix[i,] = rnorm(p,0,1)
  }
  for(i in 1:m){
    zs = z.matrix
    zs[abs(zs) < cutoff[i]] = 0
    tc = rowSums(zs^2)
    
    K1 = mean(tc)
    K2 = mean((tc-K1)^2)
    K3 = mean((tc-K1)^3)
    a = K3/(4*K2)
    b = K1-2*(K2^2)/K3
    d = 8*K2^3/K3^2
    
    para.asy[[i]] = c(a,b,d)
  }
  return(para.asy)
}

TC.test.asy <- function(interest, null.cov, cutoff, para.asy){
  ## threshold is vector of alpha
  ## null.dist is a list of null dist, length = length(threshold)
  # library(expm)
  # z.scores <- interest$summary[,3]
  # whiten.matrix <- sqrtm(solve(null.cov))
  # z.scores.decor <- t(z.scores) %*% whiten.matrix
  
  # library(whitening)
  # #decorrelate the z score, individual level data
  z <- interest$summary[,3]
  # epsilon <- sigma_fun(interest,null.cov, q)
  
  # Find W matrix
  
  w <- whitening::whiteningMatrix(null.cov, method='ZCA-cor')
  
  z.scores.decor <- z %*% t(w) 
  
  # threshold <- lapply(abs(z.scores.decor)[abs(z.scores.decor) > min(threshold)], 
  #                     function(z){max(threshold[threshold < z])})
  # threshold <- unlist(threshold)
  # threshold <- threshold[!duplicated(threshold)]
  
  m = length(cutoff)
  G.statistic <- vector("list", length = m)
  p.TC <- rep(0, m)
  
  
  for(i in 1:m){
    G.statistic[[i]] <- truncated.G(z.scores.decor, alpha=cutoff[i], gamma=2)
    
    a = para.asy[[i]][1]
    b = para.asy[[i]][2]
    d = para.asy[[i]][3]
    
    p.TC[i] <- 1 - pchisq((G.statistic[[i]]$test.statistics-b)/a, df = d)
  }
  
  select.alpha = which(p.TC == min(p.TC))
  select.alpha = select.alpha[1]
  test.statistics = p.TC[[select.alpha]]
  select = G.statistic[[select.alpha]]$select
  return(result=list(test.statistics = test.statistics, select = select))
  
}


null.dist.TC.test.asy <- function(nsim, p, cutoff, para.asy){
  vector = rep(0, nsim)
  for(i in 1:nsim){
    z.m <- rnorm(p,0,1)
    obj <- TC.test.asy(interest=list(summary=cbind(z.m,z.m,z.m)), 
                       null.cov = diag(rep(1,p)), 
                       cutoff, para.asy)
    
    vector[i] <- obj$test.statistics
  }
  return(vector)
}

null.dist.combined.test.v2 <- function(nsim,p,alpha_grid){
  cutoff = sapply(alpha_grid, function(x){
    abs(qnorm(x/2, 0, 1))
  })
  
  z.matrix = matrix(rnorm(nsim*p), nrow=nsim) # each row is one iteration
  null.cov = diag(1, nrow=p, ncol=p)
  p.value.HC <- p.value.TC <- p.value.combined <- rep(0,nsim)
  tc.matrix <- p.value.TC.cutoff <- matrix(0, nrow=nsim, ncol = length(alpha_grid)) # row is iteration and column is alpha grid
  tc.chi.para.matrix = matrix(0, nrow=3, ncol = length(alpha_grid))
  
  for(i in 1:length(alpha_grid)){
    z.temp <- z.matrix 
    z.temp[ abs(z.temp) < cutoff[i] ] = 0
    tc.temp <- rowSums(z.temp)
    
    K1 = mean(tc.temp)
    K2 = mean((tc.temp-K1)^2)
    K3 = mean((tc.temp-K1)^3)
    tc.chi.para.matrix[1,i] = K3/(4*K2)
    tc.chi.para.matrix[2,i] = K1-2*(K2^2)/K3
    tc.chi.para.matrix[3,i] = 8*K2^3/K3^2
    
    tc.matrix[,i] = rowSums(z.temp)
  }
  
  
  
  for(i in 1:nsim){
    # print(i)
    z <- z.matrix[i,]
    obj <- fgss.pleio.binary(list(summary=cbind(z,z,z)), null.cov, alpha_grid=alpha_grid, test='HC')
    p.value.HC[i] = null.stat.theoretical(p, alpha_grid, obj$statistic)
    
    for(j in 1:length(alpha_grid)){
      p.value.TC.cutoff[i,j] = 1 - pchisq((tc.matrix[i,j]-tc.chi.para.matrix[2,j])/tc.chi.para.matrix[1,j], df = tc.chi.para.matrix[3,j])
    }
    
  }
  TC.cutoff.dist <- apply(p.value.TC.cutoff,1,min)
  p.value.TC <- order(TC.cutoff.dist) / length(TC.cutoff.dist)
  p.value.combined <- apply(rbind(p.value.HC, p.value.TC), 2, min)
  
  
  return(list=list(tc.chi.para.matrix=tc.chi.para.matrix, 
                   TC.cutoff.dist=TC.cutoff.dist,
                   p.value.combined=p.value.combined))
  
}


TC.test.v2 <- function(interest.data, null.cov, null.dist.TraitScan, gamma_grid){
  TC.parameters <- null.dist.TraitScan$tc.chi.para.matrix
  
  # library(whitening)
  z <- interest.data$summary[,3]
  w <- whitening::whiteningMatrix(null.cov, method='ZCA-cor')
  z.decor <- z %*% t(w) 
  
  m = length(gamma_grid)
  G.statistic <- vector("list", length = m)
  H.statistic <- rep(0, m)
  for(i in 1:m){
    G.statistic[[i]] <- truncated.G(z.decor, alpha=gamma_grid[i], gamma=2)
    H.statistic[i] <- 1 - pchisq((G.statistic[[i]]$test.statistics-TC.parameters[2,i])/TC.parameters[1,i], df = TC.parameters[3,i])
  }
  
  select.gamma = which(H.statistic == min(H.statistic))
  select.gamma = select.gamma[1]
  test.statistics = H.statistic[[select.gamma]]
  select = G.statistic[[select.gamma]]$select
  
  pval = sum(min(H.statistic) < null.dist.TraitScan$TC.cutoff.dist) / length(null.dist.TraitScan$TC.cutoff.dist)
  return(list("gamma"=gamma_grid[select.gamma], "select"=select, "pval"=pval))
}

null.dist.Monte.Carlo <- function(N, p, alpha_grid,para=F){
  gamma_grid = lapply(alpha_grid, function(x){
    abs(qnorm(x/2, 0, 1))
  })
  gamma_grid = unlist(gamma_grid)
  z.matrix <- matrix(rnorm(p*N,0,1), nrow=N)
  
  
  if(para==F){
    null_dist_TC_prior <- vector("list", length = length(alpha_grid))
    vector_G = rep(0, N)
    
    # TC priority function
    for(j in 1:length(alpha_grid)){
      for(i in 1:N){
        obj <- truncated.G(z.matrix[i,], gamma=gamma_grid[j], power=2)
        vector_G[i] <- obj$test.statistics
      }
      null_dist_TC_prior[[j]] <- vector_G
    }
    
    # TC min p
    null_dist_TC_stat <- null_dist_TC_pval <- rep(0, N)
    for(i in 1:N){
      obj <- truncated.power.test(interest=list(summary=cbind(z.matrix[i,],z.matrix[i,],z.matrix[i,])), 
                                  null.cov = diag(rep(1,p)), 
                                  cutoff = gamma_grid, 
                                  null_dist_TC_prior, 
                                  gamma = 2)
      null_dist_TC_stat[i] <- obj$test.statistics
    }
    
    null_dist_TC_pval = rank(null_dist_TC_stat)/N
  }
  
  # HC
  null_dist_HC_pval = rep(0, N)
  for(i in 1:N){
    obj <- fgss.pleio.binary(list(summary=cbind(z.matrix[i,],z.matrix[i,],z.matrix[i,])), 
                             null.cov = diag(rep(1,p)), alpha_grid=alpha_grid, test='HC')
      
    null_dist_HC_pval[i] <- null.stat.theoretical(p, alpha_grid, obj$statistic)
  }
  
  
  
  null_dist_combined_test = pmin(null_dist_TC_pval, null_dist_HC_pval)

  return(list("null_dist_TC_prior"=null_dist_TC_prior, 
              "null_dist_TC_stat"=null_dist_TC_stat, 
              "null_dist_combined_test"=null_dist_combined_test))
}

t2cor <- function (t, n, min.n = max(n, 30)/10) 
{
  # stopifnot(length(n) == 1 || length(n) == length(p))
  # stopifnot(length(p) == length(sign))
  # t <- sign(sign) * qt(p/2, df = n - 2, lower.tail = F)
  invalid <- n < min.n
  if (any(invalid)) {
    warning(paste(sum(invalid), "statistics has n < ", 
                  min.n, "and are coded as NA."))
  }
  t[invalid] <- NA
  return(t/sqrt(n - 2 + t^2))
}

null.mixed.mc <- function(n,p_cont, p_bin,maf,variance,threshold,q,mc.cores){
  p <- p_cont + p_bin
  summary_data <- array(dim=c(q,p,5))
  maf_all <- rep(maf,q)

  Zscore.null <- parallel::mclapply(1:q, function(q_temp){
    pleio.data.sim(n,beta_cont=rep(0,p_cont),beta_bin=rep(0,p_bin),maf,variance,threshold)$summary[,3]
  }, mc.cores = mc.cores)
  Zscore.null <- do.call("rbind",Zscore.null)

  return(cor(Zscore.null))
}

null.dist.Monte.Carlo.mc <- function(N, p, alpha_grid,para=F, mc.cores = NULL){
  gamma_grid = lapply(alpha_grid, function(x){
    abs(qnorm(x/2, 0, 1))
  })
  gamma_grid = unlist(gamma_grid)
  z.matrix <- matrix(rnorm(p*N,0,1), nrow=N)
  
  

  null_dist_TC_prior <- vector("list", length = length(alpha_grid))
  vector_G = rep(0, N)
    
  # TC priority function
  null_dist_TC_prior <- parallel::mclapply(1:length(alpha_grid),function(j){
    for(i in 1:N){
      obj <- truncated.G(z.matrix[i,], gamma=gamma_grid[j], power=2)
      vector_G[i] <- obj$test.statistics
    }
    return(vector_G)
  },mc.cores=mc.cores)
    
    
  # TC min p
  null_dist_TC_stat <- null_dist_TC_pval <- rep(0, N)
  null_dist_TC_stat <- parallel::mclapply(1:N,function(i){
    obj <- truncated.power.test(interest=list(summary=cbind(z.matrix[i,],z.matrix[i,],z.matrix[i,])), 
                                null.cov = diag(rep(1,p)), 
                                cutoff = gamma_grid, 
                                null_dist_TC_prior, 
                                gamma = 2)
    return(obj$test.statistics)
  },mc.cores=mc.cores)
  null_dist_TC_stat = unlist(null_dist_TC_stat)
  null_dist_TC_pval = rank(null_dist_TC_stat)/N

  
  # HC
  null_dist_HC_pval = rep(0, N)
  null_dist_HC_pval <- parallel::mclapply(1:N,function(i){
    obj <- fgss.pleio.binary(list(summary=cbind(z.matrix[i,],z.matrix[i,],z.matrix[i,])), 
                             null.cov = diag(rep(1,p)), alpha_grid=alpha_grid, test='HC')
    return(null.stat.theoretical(p, alpha_grid, obj$statistic))
  },mc.cores=mc.cores)
  null_dist_HC_pval = unlist(null_dist_HC_pval)
  
  

  null_dist_combined_test = pmin(null_dist_TC_pval, null_dist_HC_pval)
  
  return(list("null_dist_TC_prior"=null_dist_TC_prior, 
              "null_dist_TC_stat"=null_dist_TC_stat, 
              "null_dist_combined_test"=null_dist_combined_test))
}

