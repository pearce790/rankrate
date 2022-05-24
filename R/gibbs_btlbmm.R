#' Gibbs Sampling of a BTL-B Mixture Model model
#' 
#' This function samples approximately from the posterior of a Bradley-Terry-Luce-Binomial distribution with potentially K mixture components via an adaptive Metropolis-Hastings-within-Gibbs algorithm.
#' 
#' @import stats
#' @import matrixStats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param I Numeric specifying number of judges
#' @param J Numeric specifying number of objects
#' @param Pi_full Matrix of objects considered when formulating each ranking (to account for groupwise comparisons). Should have same number of rows as Pi, and first entries of each row should match Pi. For example, if a judge considers obejcts 1,2,3,4 and provides the top-2 ranking 1<2, then Pi should have a row with entries 1,2 and Pi_full should have a row with entries 1,2,3,4.
#' @param num_iter Number of Gibbs iterations before returning results (default=1000).
#' @param sigma Covariance matrix of transformed parameters to sample, for use in restarting the Gibbs sampler after some prior function call. If null, starts with small values and independent structure (default=NULL).
#' @param p0 Starting object quality parameter values for the algorithm. (default=NULL).
#' @param theta0 Starting consensus scale parameter values for the algorithm. (default=NULL)
#' @param seed Random seed for replicability
#' @param adapt_par Vector of four values: (1) first iteration on which to begin covariance adaptation, (2) how many iterations between updates, (3) proportion of previous iterations to use when updating, and (4) proportion of total iterations before stopping updating. 
#' @param thin Numeric specifying thinning (e.g., 10 means every tenth iteration is retained).
#' @param verbose Boolean specifying if iteration information should be printed while function runs.
#' @param K Integer specifying number of latent classes (default = 1, i.e., no latent classes)
#' @param alpha0 Numeric or vector specifying hyperparameters of Dirichlet prior on latent classes. If K>1 and alpha0 is a single number, alpha0 becomes rep(alpha0,K). Default is a flat (uninformative) prior that corresponds to MLE.
#' @param a Numeric specifying first hyperparameter on Beta(a,b) prior on each p_jk. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param b Numeric specifying second hyperparameter on Beta(a,b) prior on each p_jk. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param gamma1 Numeric specifying first hyperparameter on Gamma(gamma1,gamma2) prior on each theta_k. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param gamma2 Numeric specifying second hyperparameter on Gamma(gamma1,gamma2) prior on each theta_k. Default is a flat (uninformative) prior that corresponds to MLE.
#'  
#' @return List of estimated parameters (alpha,p,theta,Z).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,5,2,1,NA,NA,NA),byrow=TRUE,nrow=2)
#' Pi_full <- matrix(c(1,2,3,4,5,2,1,3,4,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,4,1,2,2,5,5),byrow=TRUE,nrow=2)
#' gibbs_btlbmm(Pi=Pi,X=X,M=6,I=2,J=5,num_iter=100)
#' gibbs_btlbmm(Pi=Pi,X=X,M=6,I=2,J=5,a=5,b=1,gamma1=10,gamma2=1,num_iter=100)
#'  
#' @export
gibbs_btlbmm <- function(Pi,X,M,I,J,Pi_full=NULL,
                         num_iter = 1000,sigma = NULL, p0 = NULL, theta0 = NULL,
                         seed = NULL, adapt_par = c(100,50,0.50,0.50), thin = 1, verbose = FALSE,
                         K=1,alpha0=1,a=1,b=1,gamma1=1,gamma2=0){
  
  #data checks
  if(any(dim(X)!=c(I,J))){stop("Requires dim(X) = c(I,J)")}
  if(nrow(Pi)!=I){stop("Requires nrow(Pi) = I")}
  if(length(alpha0)==1){alpha0 <- rep(alpha0,K)}
  
  #fix null values, if applicable
  if(is.null(seed)){seed <- round(runif(1,0,100))}
  set.seed(seed)
  if(is.null(Pi_full)){
    Pi_full <- matrix(NA,nrow=nrow(Pi),ncol=J)
    for(i in 1:nrow(Pi)){Pi_full[i,] <- c(na.exclude(Pi[i,]),setdiff(1:J,na.exclude(Pi[i,])))}
  }else{
    if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full)")}
    if(any(Pi!=Pi_full[,1:ncol(Pi)],na.rm=T)){stop("Existing entries in Pi, Pi_full must be equivalent")}
  }
  if(is.null(sigma)){
    sigma <- array(0,dim=c(J+1,J+1,K))
    for(k in 1:K){sigma[,,k] <- c(rep(.01/J,J),.1)*diag(J+1)}
  }
  if(is.null(p0)){
    Xfilled <- X
    for(j in 1:J){Xfilled[is.na(Xfilled[,j]),j] <- mean(Xfilled[,j],na.rm=T)}
    p0 <- as.matrix(t(kmeans(na.exclude(Xfilled),K)$centers/M))
    p0[p0==0] <- 0.01
    p0[p0==1] <- 0.99
  }
  if(is.null(theta0)){
    if(gamma1==1 & gamma2==0){
      theta0 <- runif(K,2,30)
    }else{
      theta0 <- rgamma(K,gamma1,gamma2)
    }
  }
  
  #calculate other useful quantities
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  alpha.out <- matrix(NA,nrow=num_iter,ncol=K)
  q.out <- array(NA,dim=c(num_iter,J,K))
  rho.out <- matrix(NA,nrow=num_iter,ncol=K)
  Z.out <- matrix(NA,nrow=num_iter,ncol=I)
  class.probs <- array(NA,dim=c(num_iter,I,K))
  accept <- rep(0,K)
  
  alpha.out[1,] <- rep(1/K,K)
  q.out[1,,] <- log(p0/(1-p0))
  rho.out[1,] <- log(theta0)
  for(i in 1:I){
    for(k in 1:K){class.probs[1,i,k] <- log(alpha.out[1,k]) + 
      dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=p0[,k],theta=theta0[k], M=M,log=TRUE,Pi_full=matrix(Pi_full[i,],nrow=1))
    }
    class.probs[1,i,] <- exp(class.probs[1,i,] - logSumExp(class.probs[1,i,]))
    Z.out[1,i] <- sample(K,1,prob = class.probs[1,i,])
  }
  
  for(iter in 2:num_iter){
    
    #print progress if applicable
    if(verbose & iter%in%(seq(0,num_iter,length=11)[-1])){
      print(paste0("Progress: ",which(iter == seq(0,num_iter,length=11)[-1])*10,"% complete"))
    }
    
    #update covariance matrix at pre-specified intervals
    if(iter >= adapt_par[1] & iter <= num_iter*adapt_par[4] & iter%%adapt_par[2]==0){
      for(k in 1:K){
        sigma[,,k] <- cov(cbind(q.out[(iter*adapt_par[3]):(iter-1),,k],rho.out[(iter*adapt_par[3]):(iter-1),k]))
        diag(sigma[,,k])[diag(sigma[,,k])==0] <- 0.00001
      }
    }
    
    #draw alpha.out[iter,] | Z.out[iter-1,]
    alpha.out[iter,] <- c(rdirichlet(1,c(alpha0 + unlist(lapply(1:K,function(k){sum(Z.out[iter-1,]==k)})))))
    
    #draw c(q.out[iter,,k],rho.out[iter,k]) | Z.out[iter-1,]
    for(k in 1:K){
      
      whichk <- Z.out[iter-1,] == k
      while(sum(whichk) == 0){whichk <- rep(TRUE,I)}
      
      curr_trans <- c(q.out[iter-1,,k],rho.out[iter-1,k])
      prop_trans <- c(rmvnorm(n=1,mean=curr_trans,sigma=2.38/sqrt(J+1)*sigma[,,k]))
      curr <- c(exp(curr_trans[1:J])/(1+exp(curr_trans[1:J])),exp(curr_trans[J+1]))
      prop <- c(exp(prop_trans[1:J])/(1+exp(prop_trans[1:J])),exp(prop_trans[J+1]))
      
      
      loglikratio <- dbtlb(Pi=matrix(Pi[whichk,],nrow=sum(whichk)),X=matrix(X[whichk,],nrow=sum(whichk)),
                           p=prop[1:J],theta=prop[J+1],M=M,log=TRUE,
                           Pi_full=matrix(Pi_full[whichk,],nrow=sum(whichk)))-
        dbtlb(Pi=matrix(Pi[whichk,],nrow=sum(whichk)),X=matrix(X[whichk,],nrow=sum(whichk)),
              p=curr[1:J],theta=curr[J+1],M=M,log=TRUE,
              Pi_full=matrix(Pi_full[whichk,],nrow=sum(whichk)))
      logqpriorratio <- sum(dbeta(prop[1:J],a,b,log=T))-sum(dbeta(curr[1:J],a,b,log=T))
      if(gamma1==1 & gamma2==0){logrhopriorratio <- 0
      }else{logrhopriorratio <- dgamma(prop[J+1],gamma1,gamma2,log=T)-dgamma(curr[J+1],gamma1,gamma2,log=T)}
      logjacobianratio <- sum(c(prop_trans[J+1],prop_trans[1:J]-2*log(1+exp(prop_trans[1:J]))))-
        sum(c(curr_trans[J+1],curr_trans[1:J]-2*log(1+exp(curr_trans[1:J]))))
      logprob <- loglikratio + logqpriorratio + logrhopriorratio + logjacobianratio
      
      u <- runif(1)
      if(log(u) < logprob){
        q.out[iter,,k] <- prop_trans[1:J]
        rho.out[iter,k] <- prop_trans[J+1]
        if(iter>(adapt_par[4]*num_iter)){accept[k] <- accept[k] + 1}
      }else{
        q.out[iter,,k] <- curr_trans[1:J]
        rho.out[iter,k] <- curr_trans[J+1]
      }
    }
    
    #draw Z.out[iter,] | alpha.out[iter,],q.out[iter,,],rho.out[iter,]
    if(K>1){
      pcurr <- exp(q.out[iter,,])/(1+exp(q.out[iter,,]))
      thetacurr <- exp(rho.out[iter,])
      for(i in 1:I){
        for(k in 1:K){class.probs[iter,i,k] <- log(alpha.out[iter,k]) + 
          dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pcurr[,k],theta=thetacurr[k], M=M,log=TRUE,Pi_full=matrix(Pi_full[i,],nrow=1))
        }
        class.probs[iter,i,] <- exp(class.probs[iter,i,] - logSumExp(class.probs[iter,i,]))
        Z.out[iter,i] <- sample(K,1,prob = class.probs[iter,i,])
      }
    }else{
      class.probs[iter,,] <- 1
      Z.out[iter,] <- 1
    }
  }
  
  alpha.out <- alpha.out[seq(adapt_par[4]*num_iter+1,num_iter,by=thin),]
  q.out <- q.out[seq(adapt_par[4]*num_iter+1,num_iter,by=thin),,]
  rho.out <- rho.out[seq(adapt_par[4]*num_iter+1,num_iter,by=thin),]
  Z.out <- Z.out[seq(adapt_par[4]*num_iter+1,num_iter,by=thin),]
  class.probs <- class.probs[seq(adapt_par[4]*num_iter+1,num_iter,by=thin),,]
  
  
  if(verbose & K>1){print("Label-swapping and assembling final outputs")}
  if(K>1){
    labels <- stephens(class.probs)$permutations
    for(iter in 1:nrow(rho.out)){
      alpha.out[iter,] <- alpha.out[iter,labels[iter,]]
      q.out[iter,,] <- q.out[iter,,labels[iter,]]
      rho.out[iter,] <- rho.out[iter,labels[iter,]]
      Z.out[iter,] <- labels[iter,][Z.out[iter,]]
      class.probs[iter,,] <- class.probs[iter,,labels[iter,]]
    }
  }
  
  if(K>1){
    p.out <- q.out
    for(k in 1:K){p.out[,,k] <- apply(q.out[,,k],2,function(x){1/(1+exp(-x))})}
    theta.out <- apply(rho.out,2,exp)
    for(k in 1:K){sigma[,,k] <- cov(cbind(q.out[,,k],rho.out[,k]))}
  }else{
    p.out <- 1/(1+exp(-q.out))
    theta.out <- matrix(exp(rho.out),ncol=1)
    sigma[,,1] <- cov(cbind(q.out,rho.out))
  }
  
  
  return(list(alpha.out = alpha.out, p.out = p.out, theta.out = theta.out, Z.out = Z.out,
              accept = accept/(thin*nrow(rho.out)),sigma = sigma))
}
