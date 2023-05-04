#' Calculate the exact MLE of a Mallows-Binomial distribution using an A* algorithm
#'
#' This function estimates the exact MLE of a Mallows-Binomial distribution using an A* tree search algorithm proposed in Pearce and Erosheva (2022). Algorithm may be very slow when number of objects exceeds 15, but is often still tractable for larger J when consensus is strong.
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param keep_nodes Boolean specifying if function should retain a list of nodes traversed during A* tree search. Defaults to FALSE.
#'
#' @return List with elements pi0 (consensus ranking MLE),  p (object quality parameter MLE), theta (scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity). If keep_nodes=TRUE, a list of traversed nodes is also included.
#'
#' @examples
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' ratings <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' ASTAR(rankings=rankings,ratings=ratings,M=5,keep_nodes=TRUE)
#'
#' @export
ASTAR <- function(rankings,ratings,M,keep_nodes=FALSE){

  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  Q <- getQ(rankings,I,J)

  open <- matrix(NA,nrow=J,ncol=J+1)
  open[,1] <- 1:J
  for(order in 1:J){open[order,J+1] <- totalcostheuristic_MB(Q,rankings,ratings,M,order)$totalcostheuristic}
  closed <- matrix(NA,nrow=0,ncol=J+1)

  curr_node <- na.exclude(open[which.min(open[,J+1]),1:J])
  continue <- TRUE
  while(continue){
    closed <- rbind(closed,open[which.min(open[,J+1]),])
    open <- open[-which.min(open[,J+1]),]

    next_objects <- setdiff(1:J,curr_node)
    open_new <- matrix(NA,nrow=length(next_objects),ncol=J+1)
    for(ind in 1:length(curr_node)){open_new[,ind] <- curr_node[ind]}
    open_new[,length(curr_node)+1] <- next_objects
    for(ind in 1:nrow(open_new)){
      total <- totalcostheuristic_MB(Q,rankings,ratings,M,na.exclude(open_new[ind,1:J]))
      open_new[ind,J+1] <- total$totalcostheuristic
    }
    open <- rbind(open,open_new)
    curr_node <- na.exclude(open[which.min(open[,J+1]),1:J])
    if(length(curr_node)==J){continue <- FALSE}
  }

  if(any(is.na(open[,J+1]) | is.null(open[,J+1]))){stop("ERROR IN ASTAR COST CALCULATION")}

  which_min <- which(open[,J+1]==min(open[,J+1]))
  if(length(which_min)>1){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    pi0 <- open[which_min,1:J,drop=FALSE]
    results <- t(apply(pi0,1,function(curr_node){
      total <- totalcostheuristic_MB(Q,rankings,ratings,M,curr_node)
      c(total$phat,total$thetahat)
    }))
    if(keep_nodes){
      return(list(pi0=pi0,
                  p=results[,1:J],
                  theta=results[,J+1,drop=FALSE],
                  num_nodes=nrow(open)+nrow(closed),
                  nodes = list(open=open,closed=closed)))
    }else{
      return(list(pi0=pi0,
                  p=results[,1:J],
                  theta=results[,J+1,drop=FALSE],
                  num_nodes=nrow(open)+nrow(closed)))
    }
  }else{
    pi0 <- open[which_min,1:J]
    total <- totalcostheuristic_MB(Q,rankings,ratings,M,pi0)
    results <- c(total$phat,total$thetahat)
    if(keep_nodes){
      return(list(pi0=pi0,
                  p=results[1:J],
                  theta=results[J+1],
                  num_nodes=nrow(open)+nrow(closed),
                  nodes = list(open=open,closed=closed)))
    }else{
      return(list(pi0=pi0,
                  p=results[1:J],
                  theta=results[J+1],
                  num_nodes=nrow(open)+nrow(closed)))
    }
  }

}
