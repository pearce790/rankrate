#' Mallows density function
#' 
#' This function calculates the density of observation(s) under a Mallows distribution.
#' 
#' @import stats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param pi0 Vector specifying the consensus (modal probability) ranking.
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param log Boolean indicating if loglikelihood should be returned.
#' @param Pi_full Matrix indicating which values were considered when providing a ranking. For each row, first entries should match the corresponding values in Pi. Remaining entries indicate the unranked objects. For example, if a judge only views objects 1,2,4,5 and provides ranking 1<4, then Pi_full should have row 1,4,2,5 or 1,4,5,2.
#'  
#' @return (Log) likelihood of rankings under a Mallows distribution.
#'  
#' @examples 
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' dmall(Pi=Pi,pi0=c(2,1,3,4),theta=1,log=TRUE)
#'  
#' @export
dmall <- function(Pi,pi0,theta,log=FALSE,Pi_full=NULL){
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  if(sum(apply(Pi,1,function(pi){!all(is.na(pi))}))==0){
    if(log){return(0)}else{return(1)}
  }
  
  if(is.null(Pi_full)){
    J <- length(pi0)
    R <- apply(Pi,1,function(pi){length(na.exclude(pi))})
    R[which(R == J-1)] <- J #missing one object implicitly defines a complete ranking
    
    logd <- -theta*sum(apply(Pi,1,function(pi){kendall(pi,pi0)}))-
      sum(unlist(lapply(R,function(r){psi(theta,J,r,log=T)})))
  }else{
    J <- apply(Pi_full,1,function(pi){length(na.exclude(pi))})
    R <- apply(Pi,1,function(pi){length(na.exclude(pi))})
    R[which(R == J-1)] <- R[which(R == J-1)]+1 #missing one object implicitly defines a complete ranking
    sumD <- sum(unlist(lapply(1:nrow(Pi),function(i){
      removed <- setdiff(1:length(pi0),Pi_full[i,])
      kendall(Pi[i,],c(setdiff(pi0,removed),removed))
    })))
    sumlogpsi <- sum(unlist(lapply(1:nrow(Pi),function(i){
      psi(theta,J[i],R[i],log=T)
    })))
    logd <- -theta*sumD-sumlogpsi
  }
  
  if(log){return(logd)
  }else{return(exp(logd))}
}