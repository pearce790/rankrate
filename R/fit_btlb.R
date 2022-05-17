#' MLE and MAP Estimation of a BTL-B model
#' 
#' This function estimates the maximum likelihood or maximum a posteriori estimates of a Bradley-Terry-Luce-Binomial distribution via an Expectation-Maximization (EM) algorithm.
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
#' @param tol Numeric specifying stopping tolerance of loglikelihood in successive EM iterations.
#' @param maxiter Maximum number of EM iterations before returning results (default=50).
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
#' fit_btlb(Pi=Pi,X=X,M=6,I=2,J=5,tol=0.01)
#' fit_btlb(Pi=Pi,X=X,M=6,I=2,J=5,tol=0.01,a=5,b=1,gamma1=10,gamma2=10)
#'  
#' @export
fit_btlb <- function(Pi,X,M,I,J,Pi_full=NULL,tol=1,maxiter=50,verbose=FALSE,
                     K=1,alpha0=1,a=1,b=1,gamma1=1,gamma2=0){
  
  if(any(dim(X)!=c(I,J))){stop("Requires dim(X) = c(I,J)")}
  if(length(alpha0)==1){alpha0 <- rep(alpha0,K)}
  
  #preprocessing
  if(is.null(Pi_full)){
    Pi_full <- matrix(NA,nrow=nrow(Pi),ncol=J)
    for(i in 1:nrow(Pi)){Pi_full[i,] <- c(na.exclude(Pi[i,]),setdiff(1:J,na.exclude(Pi[i,])))}
  }else{if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full)")}
  }
  
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  
  #initializing
  Zcurr <- Znew <- matrix(1/K,nrow=I,ncol=K)
  alphacurr <- alphanew <- rep(1/K,K)
  Xfilled <- X
  for(j in 1:J){Xfilled[is.na(Xfilled[,j]),j] <- mean(Xfilled[,j],na.rm=T)}
  pcurr <- pnew <- as.matrix(t(kmeans(na.exclude(Xfilled),K)$centers/M))
  thetacurr <- thetanew <- rep(1,K)
  
  loglikcurr <- sum(unlist(lapply(1:I,function(i){
    logSumExp(unlist(lapply(1:K,function(k){
      log(alphacurr[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pcurr[,k],theta=thetacurr[k],
                              M=M,log=TRUE,Pi_full=matrix(Pi_full[i,],nrow=1))
    })))
  })))
  diff <- tol+1
  counter <- 1
  
  while(diff>tol){
    if(counter==maxiter){
      print("Max iteration reached. Returning current solution, which may be imprecise.")
      break()}
    
    # E-Step
    for(i in 1:I){for(k in 1:K){
      Znew[i,k] <- log(alphacurr[k]) + dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pcurr[,k],theta=thetacurr[k],
                                             M=M,log=TRUE,Pi_full=matrix(Pi_full[i,],nrow=1))
    }}
    Znew <- exp(Znew-apply(Znew,1,logSumExp))
    
    
    # M-Step
    alphanew <- (alpha0 - 1 + apply(Znew,2,sum))/(I-K+sum(alpha0))
    for(k in 1:K){
      
      D <- a-1+apply(Znew[,k]*X,2,function(x){sum(x,na.rm=T)})
      E <- b-1+apply(Znew[,k]*(M-X),2,function(x){sum(x,na.rm=T)})
      
      loglik <- function(param,J,Pi,D,E,Pi_full,Ri,Z,k){
        p <- param[1:J]
        theta <- param[J+1]
        exp <- -theta*p
        omega <- exp(exp)
        
        term1 <- -theta*(gamma2+sum(Z[,k]*p[Pi],na.rm=TRUE))
        term2 <- (gamma1-1)*log(theta)
        term3 <- -sum(unlist(lapply(1:nrow(Pi),function(i){
          if(Ri[i]==0){0
          }else{
            Z[i,k]*sum(unlist(lapply(1:Ri[i],function(r){
              logSumExp(exp[Pi_full[i,r:length(na.exclude(Pi_full[i,]))]])
            })))
          }
        })))
        term4 <- sum(D*log(p)+E*log(1-p))
        
        return(-(term1+term2+term3+term4))
      }
      
      est <- optim(c(pcurr[,k],thetacurr[k]),loglik,method="L-BFGS-B",lower=rep(1e-5,J+1),upper=c(rep(1-1e-5,J),Inf),
                   J=J,Pi=Pi,D=D,E=E,Pi_full=Pi_full,Ri=Ri,Z=Znew,k=k)$par
      pnew[,k] <- est[1:J]
      thetanew[k] <- est[J+1]
    }
    
    logliknew <- sum(unlist(lapply(1:I,function(i){
      logSumExp(unlist(lapply(1:K,function(k){
        log(alphanew[k])+dbtlb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pnew[,k],theta=thetanew[k],
                               M=M,log=TRUE,Pi_full=matrix(Pi_full[i,],nrow=1))
      })))
    })))
    (diff <- abs(logliknew-loglikcurr))
    alphacurr <- alphanew
    pcurr <- pnew
    thetacurr <- thetanew
    Zcurr <- Znew
    loglikcurr <- logliknew
    if(verbose){print(paste0("Iteration No. ",counter,". Abs. diff. in loglik from previous iteration: ",diff))}
    counter <- counter + 1
  }
  
  return(list(alpha=alphacurr,p=pcurr,theta=thetacurr,Z=Zcurr))
}