#' BTL-Binomial density function
#' 
#' This function calculates the density of observation(s) under a Bradley-Terry-Luce-Binomial (BTL-B) distribution.
#' 
#' @import stats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param p Vector of object qualities.
#' @param theta Numeric specifying the consensus strength.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param log Boolean indicating if loglikelihood should be returned.
#' @param Pi_full Matrix of objects considered when formulating each ranking (to account for groupwise comparisons). Should have same number of rows as Pi, and first entries of each row should match Pi. For example, if a judge considers obejcts 1,2,3,4 and provides the top-2 ranking 1<2, then Pi should have a row with entries 1,2 and Pi_full should have a row with entries 1,2,3,4.

#'  
#' @return (Log) likelihood of rankings and ratings under a BTL-B distribution.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,5,2,1,NA,NA,NA),byrow=TRUE,nrow=2)
#' Pi_full <- matrix(c(1,2,3,4,5,2,1,3,4,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,4,1,2,2,5,5),byrow=TRUE,nrow=2)
#' p <- c(.1,.2,.5,.7,.8)
#' theta <- 10
#' M <- 6
#' dbtlb(Pi,X,p,theta,M,log=TRUE)
#' dbtlb(Pi,X,p,theta,M,log=TRUE,Pi_full)
#'  
#' @export
dbtlb <- function(Pi,X,p,theta,M,log=FALSE,Pi_full=NULL){
  if(!is.matrix(X)){stop("X must be a matrix of ratings")}
  if(length(p)!=ncol(X)){stop("length(p) must equal ncol(X)")}
  worth <- exp(-theta*p)
  logd <- dbtl(Pi,worth,log=T,Pi_full=Pi_full)+sum(apply(X,1,function(x){dbinom(x,M,p,log=T)}),na.rm=T)
  if(log){return(logd)
  }else{return(exp(logd))}
}