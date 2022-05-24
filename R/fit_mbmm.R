#' MLE and MAP Estimation of a MB Mixture Model model
#' 
#' This function estimates the maximum likelihood or maximum a posteriori estimates of a Mallows-Binomial distribution with potentially K mixture components via an Expectation-Maximization (EM) algorithm.
#' 
#' @import stats
#' @import matrixStats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param I Numeric specifying number of judges
#' @param J Numeric specifying number of objects
#' @param Pi_full Matrix of objects considered when formulating each ranking (to account for limited object access). Should have same number of rows as Pi, and first entries of each row should match Pi. For example, if a judge considers obejcts 1,2,3,4 and provides the top-2 ranking 1<2, then Pi should have a row with entries 1,2 and Pi_full should have a row with entries 1,2,3,4.
#' @param tol Numeric specifying stopping tolerance of loglikelihood in successive EM iterations.
#' @param maxiter Maximum number of EM iterations before returning results (default=50).
#' @param scope Number of rankings to sample at each update stage in the EM algorithm. When scope is large, algorithm will be slower but more likely to obtain an exact solution.
#' @param verbose Boolean specifying if iteration information should be printed while function runs.
#' @param K Integer specifying number of latent classes (default = 1, i.e., no latent classes)
#' @param alpha0 Numeric or vector specifying hyperparameters of Dirichlet prior on latent classes. If K>1 and alpha0 is a single number, alpha0 becomes rep(alpha0,K). Default is a flat (uninformative) prior that corresponds to MLE.
#' @param a Numeric specifying first hyperparameter on Beta(a,b) prior on each p_jk. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param b Numeric specifying second hyperparameter on Beta(a,b) prior on each p_jk. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param gamma1 Numeric specifying first hyperparameter on Gamma(gamma1,gamma2) prior on each theta_k. Default is a flat (uninformative) prior that corresponds to MLE.
#' @param gamma2 Numeric specifying second hyperparameter on Gamma(gamma1,gamma2) prior on each theta_k. Default is a flat (uninformative) prior that corresponds to MLE.
#'  
#' @return List of estimated parameters (alpha,p,pi0,theta,Z).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,5,2,1,NA,NA,NA),byrow=TRUE,nrow=2)
#' Pi_full <- matrix(c(1,2,3,4,5,2,1,3,4,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,4,1,2,2,5,5),byrow=TRUE,nrow=2)
#' fit_mbmm(Pi=Pi,X=X,M=6,I=2,J=5,tol=0.01)
#' fit_mbmm(Pi=Pi,X=X,M=6,I=2,J=5,tol=0.01,a=5,b=1,gamma1=10,gamma2=10)
#'  
#' @export

fit_mbmm <- function(Pi,X,M,I,J,Pi_full=NULL,tol=1,maxiter=50,scope=max(10,log(factorial(J))),verbose=FALSE,
                     K=1,alpha0=1,a=1,b=1,gamma1=1,gamma2=0){
  
  if(any(dim(X)!=c(I,J))){stop("Requires dim(X) = c(I,J)")}
  if(length(alpha0)==1){alpha0 <- rep(alpha0,K)}
  
  #functions for updating phat and thetahat in EM algorithm
  phat_unconditional <- apply(X,2,function(x){mean(x,na.rm=T)})/M
  phat_unconditional[phat_unconditional<0.01] <- 0.01
  phat_unconditional[phat_unconditional>0.99] <- 0.99
  phat_conditional_map <- function(D,E,phat,order,J){
    
    opt_f <- function(p,D,E,order,J){
      # quantity that must be optimized with respect to p
      return(-sum(D*log(p)+E*log(1-p)))
    }
    opt_g <- function(p,D,E,order,J){
      # contraints, which are met when all quantities returned are nonpositive
      return(c(-diff(p[order]),p[order[length(order)]]-p[setdiff(1:J,order)]))
    }
    res <- nloptr(x0=phat,
                  eval_f = opt_f, eval_g_ineq = opt_g,
                  lb = rep(1e-8,J), ub = rep(1-1e-8,J),
                  D=D,E=E,order=order,J=J,
                  opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_abs"=1e-8,maxeval=1000))
    list(phat=res$solution,obj=res$objective)
  }
  theta_conditional_map <- function(G,H,Pi,Z,k,Ji,Ri){
    tmp_fn <- function(theta){
      -(-theta*G+log(theta)*H - sum(unlist(lapply(1:nrow(Pi),function(i){
        if(Ri[i]==0){return(0)
        }else{return(Z[i,k]*psi(theta,Ji[i],Ri[i],log=T))}
      }))))
    }
    res <- optim(par=1,tmp_fn,lower=1e-5,upper=Inf,method="L-BFGS-B")
    list(thetahat = res$par, obj = res$value)
  }
  
  
  #preprocessing
  if(is.null(Pi_full)){
    Pi_full <- matrix(NA,nrow=nrow(Pi),ncol=J)
    for(i in 1:nrow(Pi)){Pi_full[i,] <- c(na.exclude(Pi[i,]),setdiff(1:J,na.exclude(Pi[i,])))}
  }else{
    if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full)")}
    if(any(Pi!=Pi_full[,1:ncol(Pi)],na.rm=T)){stop("Existing entries in Pi, Pi_full must be equivalent")}
  }
  
  Ri <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  Ji <- apply(Pi_full,1,function(pi){length(na.exclude(pi))})
  
  #initializing
  Zcurr <- Znew <- matrix(1/K,nrow=I,ncol=K)
  alphacurr <- alphanew <- rep(1/K,K)
  Xfilled <- X
  for(j in 1:J){Xfilled[is.na(Xfilled[,j]),j] <- mean(Xfilled[,j],na.rm=T)}
  pcurr <- pnew <- as.matrix(t(kmeans(na.exclude(Xfilled),K)$centers/M))
  pcurr[pcurr==0] <- 1e-5
  pcurr[pcurr==1] <- 1-1e-5
  pi0curr <- pi0new <- as.matrix(apply(pcurr,2,order))
  thetacurr <- thetanew <- rep(1,K)
  
  loglikcurr <- sum(unlist(lapply(1:I,function(i){
    logSumExp(unlist(lapply(1:K,function(k){
      log(alphacurr[k])+dmb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pcurr[,k],
                            pi0=pi0curr[,k],theta=thetacurr[k],M=M,log=T,Pi_full=matrix(Pi_full[i,],nrow=1))
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
      Znew[i,k] <- log(alphacurr[k])+dmb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pcurr[,k],
                                         pi0=pi0curr[,k],theta=thetacurr[k],M=M,log=T,Pi_full=matrix(Pi_full[i,],nrow=1))
    }}
    Znew <- exp(Znew-apply(Znew,1,logSumExp))
    
    
    # M-Step
    alphanew <- (alpha0 - 1 + apply(Znew,2,sum))/(I-K+sum(alpha0))
    for(k in 1:K){
      
      D <- a-1+apply(Znew[,k]*X,2,function(x){sum(x,na.rm=T)})
      E <- b-1+apply(Znew[,k]*(M-X),2,function(x){sum(x,na.rm=T)})
      
      pi0_options <- unique(rbind(pi0curr[,k],order(D/(D+E)),
                                  t(replicate(scope/2,order(pcurr[,k]+rnorm(J,mean=0,sd=sd(diff(sort(pcurr[,k])))/3)))),
                                  t(replicate(scope/2,order(D/(D+E)+rnorm(J,mean=0,sd=sd(diff(sort(D/(D+E))))/3))))
      ))
      objk <- Inf
      
      for(ind in 1:nrow(pi0_options)){
        G <- gamma2 + sum(Znew[,k]*unlist(lapply(1:I,function(i){
          if(length(na.exclude(Pi[i,]))==0){return(0)
          }else{
            removed <- setdiff(1:J,Pi_full[i,])
            return(kendall(Pi[i,],c(setdiff(pi0_options[ind,],removed),removed)))
            }
        })))
        H <- gamma1-1
        ptmp <- phat_conditional_map(D=D,E=E,phat=phat_unconditional,order=pi0_options[ind,],J)
        thetatmp <- theta_conditional_map(G=G,H=H,Pi=Pi,Z=Znew,k=k,Ji=Ji,Ri=Ri)
        if(ptmp$obj+thetatmp$obj < objk){
          pnew[,k] <- ptmp$phat
          thetanew[k] <- thetatmp$thetahat
          pi0new[,k] <- pi0_options[ind,]
          objk <- ptmp$obj+thetatmp$obj
        }
      }
    }
    
    logliknew <- sum(unlist(lapply(1:I,function(i){
      logSumExp(unlist(lapply(1:K,function(k){
        log(alphanew[k])+dmb(Pi=matrix(Pi[i,],nrow=1),X=matrix(X[i,],nrow=1),p=pnew[,k],
                             pi0=pi0new[,k],theta=thetanew[k],M=M,log=T,Pi_full=matrix(Pi_full[i,],nrow=1))
      })))
    })))
    (diff <- abs(logliknew-loglikcurr))
    alphacurr <- alphanew
    pcurr <- pnew
    pi0curr <- pi0new
    thetacurr <- thetanew
    Zcurr <- Znew
    loglikcurr <- logliknew
    if(verbose){print(paste0("Iteration No. ",counter,". Abs. diff. in loglik from previous iteration: ",diff))}
    counter <- counter + 1
  }
  
  return(list(alpha=alphacurr,p=pcurr,pi0=pi0curr,theta=thetacurr,Z=Zcurr))
}
