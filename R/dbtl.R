#' Bradley-Terry-Luce density function
#' 
#' This function calculates the density of observation(s) under a Bradley-Terry-Luce distribution.
#' 
#' @import stats
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param worth Vector of object worth parameters from a Bradley-Terry-Luce type distribution.
#' @param log Boolean indicating if loglikelihood should be returned.
#' @param Pi_full Matrix of objects considered when formulating each ranking (to account for groupwise comparisons). Should have same number of rows as Pi, and first entries of each row should match Pi. For example, if a judge considers obejcts 1,2,3,4 and provides the top-2 ranking 1<2, then Pi should have a row with entries 1,2 and Pi_full should have a row with entries 1,2,3,4.
#'  
#' @return (Log) likelihood of rankings under a Bradley-Terry-Luce distribution.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,5,2,1,NA,NA,NA),byrow=TRUE,nrow=2)
#' Pi_full <- matrix(c(1,2,3,4,5,2,1,3,4,NA),byrow=TRUE,nrow=2)
#' dbtl(Pi=Pi,worth=c(5,4,3,2,1),log=TRUE)
#' dbtl(Pi=Pi,worth=c(5,4,3,2,1),log=TRUE,Pi_full=Pi_full)
#'  
#' @export
dbtl <- function(Pi,worth,log=FALSE,Pi_full=NULL){
  if(!is.matrix(Pi)){stop("Pi must be a matrix of (partial) rankings")}
  if(!is.null(Pi_full)){if(nrow(Pi)!=nrow(Pi_full)){stop("nrow(Pi) must equal nrow(Pi_full).")}}
  if(length(worth)<max(Pi,Pi_full,na.rm=T)){stop("Ensure one worth parameter for every object in Pi and/or Pi_full")}
  
  if(sum(apply(Pi,1,function(pi){!all(is.na(pi))}))==0){ #if all rankings empty, return likelihood of 1 (convenience)
    if(log){return(0)}else{return(1)}
  }
  
  if(!is.null(Pi_full)){
    logd <- sum(unlist(lapply(1:nrow(Pi),function(i){
      pi <- na.exclude(Pi[i,])
      if(length(pi)==0){return(0)
      }else{
        num <- worth[pi]
        denom_terms <- na.exclude(worth[Pi_full[i,]])
        return(sum(log(num))-sum(unlist(lapply(1:length(num),function(j){log(sum(denom_terms[j:length(denom_terms)]))}))))
      }
    })))
  }else{
    logd <- sum(apply(Pi,1,function(pi){
      if(length(na.exclude(pi))==0){return(0)
      }else{
        pi <- na.exclude(pi)
        num <- worth[pi]
        denom_terms <- c(worth[pi],worth[setdiff(1:length(worth),pi)])
        return(sum(log(num)-log(rev(cumsum(rev(denom_terms)))[1:length(pi)])))
      }
    }))
  }
  
  if(log){return(logd)
  }else{return(exp(logd))}
}