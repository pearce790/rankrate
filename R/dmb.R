#' Mallows-Binomial density function
#' 
#' This function calculates the density of observation(s) under a Mallows-Binomial distribution.
#' 
#' @import stats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param p Vector of object qualities.
#' @param pi0 Vector specifying the consensus (modal probability) ranking (useful to break ties in p).
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param log Boolean indicating if loglikelihood should be returned.
#' @param Pi_full For rankings based on incomplete object access. See documentation for dmall.
#'  
#' @return (Log) likelihood of rankings under a Mallows-Binomial distribution.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' dmb(Pi=Pi,X=X,p=c(.1,.2,.5,.9),theta=1.1,M=5,log=TRUE)
#' dmb(Pi=Pi,X=X,p=c(.1,.2,.5,.9),theta=2,M=5,log=TRUE)
#'  
#' @export
dmb <- function(Pi,X,p,pi0=NULL,theta,M,log=FALSE,Pi_full=NULL){
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(length(p)!=ncol(X)){stop("p must equal ncol(X)")}
  if(is.null(pi0)){pi0 <- order(p)}
  logd <- dmall(Pi,pi0,theta,log=T,Pi_full)+sum(apply(X,1,function(x){dbinom(x,M,p,log=T)}),na.rm=T)
  if(log){return(logd)
  }else{return(exp(logd))}
}