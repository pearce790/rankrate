#' Estimate the MLE of a Mallows-Binomial distribution using the GreedyLocal method
#'
#' This function estimates the MLE of a Mallows-Binomial distribution using the GreedyLocal method, which is identical to the Greedy method but includes an automatic and targeted post-hoc local search.
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
#' ratings <- matrix(c(2,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' GreedyLocal(rankings=rankings,ratings=ratings,M=5)
#'
#' @export
GreedyLocal <- function(rankings,ratings,M){

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

  for(row in 1:nrow(curr_ranking)){
    nearby_rankings <-  matrix(unlist(lapply(1:(J-1),function(j){
      new <- curr_ranking[row,]
      new[j:(j+1)] <- curr_ranking[row,(j+1):j]
      new
    })),byrow=T,ncol=J)
    curr_ranking <- rbind(curr_ranking,nearby_rankings)
  }
  curr_ranking <- unique(curr_ranking)
  results <- matrix(data=NA,nrow=nrow(curr_ranking),ncol=J+2)
  for(row in 1:nrow(curr_ranking)){
    curr_node <- curr_ranking[row,]
    results[row,] <- unlist(totalcostheuristic_MB(Q,rankings,ratings,M,curr_node))
  }
  num_nodes <- num_nodes + nrow(curr_ranking)
  which_keep <- which(results[,J+2] == min(results[,J+2]))
  curr_ranking <- curr_ranking[which_keep,,drop=FALSE]

  if(nrow(curr_ranking)>1){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    return(list(pi0=curr_ranking,
                p=results[which_keep,1:J],
                theta=results[which_keep,J+1,drop=FALSE],
                num_nodes=num_nodes))
  }else{
    return(list(pi0=c(curr_ranking),
                p=results[which_keep,1:J],
                theta=results[which_keep,J+1],
                num_nodes=num_nodes))
  }
}
