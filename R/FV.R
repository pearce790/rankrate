#' Estimate the MLE of a Mallows-Binomial distribution using the FV method
#' 
#' This function estimates the MLE of a Mallows-Binomial distribution using the FV method.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param localsearch Numberic specifying the maximum Kendall distance to the first-estimated consensus ranking of rankings which should be considered during a post-hoc local search; defaults to 0 (indicating no local search). 
#'  
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' FV(Pi=Pi,X=X,M=5)
#' FV(Pi=Pi,X=X,M=5,localsearch=2)
#'  
#' @export
FV <- function(Pi,X,M,localsearch=0){
  J <- ncol(X)
  Pi <- Pi[apply(Pi,1,function(pi){!all(is.na(pi))}),]
  
  central_rank <- order(apply(getQ(Pi,J),2,sum),partial=apply(X,2,function(score){sum(score,na.rm=T)}))
  central_score<- order(apply(X,2,function(score){sum(score,na.rm=T)}),partial=apply(getQ(Pi,J),2,sum))
  
  rankings <- unique(matrix(c(central_rank,central_score),byrow=T,nrow=2))
  if(localsearch>0){
    for(localsearch_iter in 1:localsearch){
      for(row in 1:nrow(rankings)){
        rankings <- rbind(rankings,matrix(unlist(lapply(1:(J-1),function(j){
          new <- rankings[row,]
          new[j:(j+1)] <- rankings[row,(j+1):j]
          new
        })),byrow=T,ncol=J))
      }
      rankings <- unique(rankings)
    }
  }
  
  cost_parameters <- matrix(NA,nrow=nrow(rankings),ncol=2+J)
  for(order_index in 1:nrow(rankings)){
    order <- rankings[order_index,]
    cost_parameters[order_index,1:J] <- phat <- phat_conditional(X,M,order)
    cost_parameters[order_index,J+1] <-thetahat <- theta_conditional(Pi,order)
    cost_parameters[order_index,J+2] <- -dmb(Pi,X,phat,order,thetahat,M,log=T)
  }
  
  cost_index <- which.min(cost_parameters[,J+2])
  cost_index <- which(cost_parameters[,J+2] == cost_parameters[cost_index,J+2])
  
  return(list(pi0=rankings[cost_index,],p=cost_parameters[cost_index,1:J],
              theta=cost_parameters[cost_index,J+1],num_nodes = nrow(rankings)))
}