#' Estimate the MLE of a Mallows-Binomial distribution using the Greedy method
#'
#' This function estimates the MLE of a Mallows-Binomial distribution using the Greedy method.
#'
#' @import gtools
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#'
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'
#' @examples
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' ratings <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' Greedy(rankings=rankings,ratings=ratings,M=5)
#'
#' @export
Greedy <- function(rankings,ratings,M){

  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  Q <- getQ(rankings,I,J)

  num_nodes <- 0
  curr_ranking <- matrix(data=NA,nrow=1,ncol=0)
  while(ncol(curr_ranking)<J){
    tmp <- matrix(data=NA,nrow=0,ncol=ncol(curr_ranking)+1)
    for(node_index in 1:nrow(curr_ranking)){
      node <- curr_ranking[node_index,]
      S <- setdiff(1:J,node)
      cost <- rep(NA,length(S))
      for(index in 1:length(S)){
        num_nodes <- num_nodes + 1
        try_ranking <- c(node,S[index])
        total <- totalcostheuristic_MB(Q,rankings,ratings,M,try_ranking)
        cost[index] <- total$totalcostheuristic
      }
      whichmin <- which(cost==min(cost))
      for(which in whichmin){tmp <- rbind(tmp,c(node,S[which]))}
    }
    curr_ranking <- tmp
  }

  if(nrow(curr_ranking)>1){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    pi0 <- curr_ranking
    results <- t(apply(pi0,1,function(curr_node){
      total <- totalcostheuristic_MB(Q,rankings,ratings,M,curr_node)
      c(total$phat,total$thetahat)
    }))
    return(list(pi0=pi0,
                p=results[,1:J],
                theta=results[,J+1,drop=FALSE],
                num_nodes=num_nodes))
  }else{
    pi0 <- c(curr_ranking)
    total <- totalcostheuristic_MB(Q,rankings,ratings,M,pi0)
    return(list(pi0=pi0,
                p=total$phat,
                theta=total$thetahat,
                num_nodes=num_nodes))
  }
}
