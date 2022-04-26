#' Calculate Confidence Intervals for Mallows-Binomial parameters.
#' 
#' This function calculates confidence intervals for parameters in a Mallows-Binomial model using the nonparametric bootstrap.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param interval Numeric between 0 and 1 specifying the confidence interval (e.g., .90 indicates a 90% confidence interval). Defaults to 0.90.
#' @param nsamples Numeric indicating desired number of bootstrap samples to be used when calculating confidence intervals. Defaults to 100.
#' @param all Boolean indicating if estimated parameters from all bootstrap samples should be returned. Defaults to FALSE.
#' @param method String indicating which estimation method to use when estimating parameters.
#' @param localsearch Numeric for use with the Greedy and FV methods; see documentation of those estimation functions for details. Defaults to 0 indicating no local search.
#' 
#' @return List with elements ci (matrix of confidence intervals for Mallows-Binomial parameters; always returned), bootstrap_pi0 (matrix of bootstrap consensus rankings; returned only if all=TRUE), and bootstrap_ptheta (matrix of bootstrap estimates of (p,theta); returned only if all=TRUE).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' ci_mb(Pi=Pi,X=X,M=5,method="ASTAR")
#' ci_mb(Pi=Pi,X=X,M=5,method="FV",interval=0.95,nsamples=200,all=TRUE,localsearch=1)
#'  
#' @export
ci_mb <- function(Pi,X,M,interval=0.90,nsamples=100,all=FALSE,
                  method=c("ASTAR","ASTAR_LP","Greedy","GreedyLocal","FV"),localsearch=0){
  I <- max(nrow(X),nrow(Pi))
  J <- ncol(X)
  if(all){bs_pi0 <- matrix(NA,nrow=nsamples,ncol=J)}
  bs_parameters <- matrix(NA,nrow=nsamples,ncol=J+1)
  for(sample in 1:nsamples){
    bs_sample <- sample(1:I,I,replace=T)
    bs_res <- fit_mb(Pi=Pi[bs_sample,],X=X[bs_sample,],M=M,method=method,localsearch=localsearch)
    if(all){bs_pi0[sample,] <- bs_res$pi0}
    bs_parameters[sample,] <- c(bs_res$p,bs_res$theta)
  }
  ci <- as.data.frame(apply(bs_parameters,2,function(param){quantile(param,probs=c((1-interval)/2,1-(1-interval)/2))}))
  names(ci) <- c(paste0("p",1:J),"theta")
  
  if(all){return(list(ci=ci,bootstrap_pi0=bs_pi0,boostrap_ptheta=bs_parameters))
  }else{return(list(ci=ci))}
}