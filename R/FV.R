#' Estimate the MLE of a Mallows-Binomial distribution using the FV method
#'
#' This function estimates the MLE of a Mallows-Binomial distribution using the FV method.
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
#' FV(rankings=rankings,ratings=ratings,M=5)
#'
#' @export
FV <- function(rankings,ratings,M){

  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  Q <- getQ(rankings,I,J)

  Qsum <- apply(Q,2,sum)
  central_ranking <- matrix(order(Qsum),nrow=1)
  for(j in 2:J){
    if(Qsum[central_ranking[1,j-1]]==Qsum[central_ranking[1,j]]){
      swapped_rank <- central_ranking[1,]
      swapped_rank[c(j-1,j)] <- swapped_rank[c(j,j-1)]
      central_ranking <- rbind(central_ranking,swapped_rank)
    }
  }
  mean_ratings <- apply(ratings,2,function(rating){mean(rating,na.rm=T)})
  central_rating <- matrix(order(mean_ratings),nrow=1)
  for(j in 2:J){
    if(mean_ratings[central_rating[1,j-1]]==mean_ratings[central_rating[1,j]]){
      swapped_rating <- central_rating[1,]
      swapped_rating[c(j-1,j)] <- swapped_rating[c(j,j-1)]
      central_rating <- rbind(central_rating,swapped_rating)
    }
  }
  rankings <- unique(rbind(central_rating,central_ranking))
  rankings_info <- matrix(NA,nrow=nrow(rankings),ncol=J+2)
  for(index in 1:nrow(rankings)){
    node <- rankings[index,]
    total <- totalcostheuristic_MB(Q,rankings,ratings,M,node)
    rankings_info[index,] <- c(total$phat,total$thetahat,total$totalcostheuristic)
  }
  row.names(rankings) <- NULL
  which_min <- which(rankings_info[,J+2] == min(rankings_info[,J+2]))
  if(length(which_min)>1){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    return(list(pi0=rankings[which_min,],
         p=rankings_info[which_min,1:J],
         theta=rankings_info[which_min,J+2,drop=FALSE],
         num_nodes=nrow(rankings)))
  }else{
    return(list(pi0=rankings[which_min,],
                p=rankings_info[which_min,1:J],
                theta=rankings_info[which_min,J+2],
                num_nodes=nrow(rankings)))
  }
}
