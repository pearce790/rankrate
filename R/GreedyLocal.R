#' Estimate the MLE of a Mallows-Binomial distribution using the GreedyLocal method
#' 
#' This function estimates the MLE of a Mallows-Binomial distribution using the GreedyLocal method, which is identical to the Greedy method but includes an automatic and targeted post-hoc local search.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' 
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' GreedyLocal(Pi=Pi,X=X,M=5)
#'  
#' @export
GreedyLocal <- function(Pi,X,M){
  J <- ncol(X)
  Pi <- Pi[apply(Pi,1,function(pi){!all(is.na(pi))}),]
  Q <- getQ(Pi,J)
  
  curr_ranking <- c()
  num_nodes <- 0
  while(length(curr_ranking)<(J-1)){
    S <- setdiff(1:J,curr_ranking)
    cost <- rep(NA,length(S))
    for(i in 1:length(S)){
      num_nodes <- num_nodes + 1
      s <- S[i]
      try_ranking <- c(curr_ranking,s)
      cost[i] <- totalcostheuristic_MB(Q,Pi,X,M,try_ranking)
    }
    curr_ranking <- c(curr_ranking,S[which.min(cost)])
  }
  
  rankings <- unique(matrix(c(curr_ranking,setdiff(1:J,curr_ranking)),byrow=T,nrow=1))
  old_rankings <- matrix(NA,nrow=0,ncol=J)
  continue <- TRUE
  
  while(continue){
    rankings <- rbind(rankings,matrix(unlist(lapply(1:(J-1),function(j){
      new <- rankings[1,]
      new[j:(j+1)] <- rankings[1,(j+1):j]
      new
    })),byrow=T,ncol=J))
    rankings <- rankings[apply(rankings,1,function(pi1){!any(apply(old_rankings,1,function(pi2){all(pi1==pi2)}))}),]
    #remove rankings which have already been searched
    cost_parameters <- matrix(NA,nrow=nrow(rankings),ncol=2+J)
    for(order_index in 1:nrow(rankings)){
      order <- rankings[order_index,]
      cost_parameters[order_index,1:J] <- phat <- phat_conditional(X,M,order)
      cost_parameters[order_index,J+1] <-thetahat <- theta_conditional(Pi,order)
      cost_parameters[order_index,J+2] <- -dmb(Pi,X,phat,order,thetahat,M,log=T)
    }
    cost_index <- which.min(cost_parameters[,J+2])
    if(cost_index == 1){continue <- FALSE
    }else{
      old_rankings <- rbind(old_rankings,rankings[-cost_index,])
      rankings <- matrix(rankings[cost_index,],nrow=1)
    }
  }
  
  return(list(pi0=rankings[cost_index,],p=cost_parameters[cost_index,1:J],
              theta=cost_parameters[cost_index,J+1],num_nodes = num_nodes + nrow(rankings) + nrow(old_rankings)))
}