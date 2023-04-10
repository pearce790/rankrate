#' Calculate the exact MLE of a Mallows-Binomial distribution using an A* search algorithm
#' 
#' This function estimates the exact MLE of a Mallows-Binomial distribution using an A* tree search algorithm proposed in Pearce and Erosheva (2022). The algorithm may be very slow when number of objects, J, exceeds 15, but is often still tractable for larger J when ranking and rating consensus among judges is strong.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' 
#' @return List with elements pi0 (consensus ranking MLE),  p (object quality parameter MLE), theta (scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' ASTAR(Pi=Pi,X=X,M=5)
#'  
#' @export
ASTAR <- function(Pi,X,M){

  J <- ncol(X)
  Pi <- Pi[apply(Pi,1,function(pi){!all(is.na(pi))}),]
  Q <- getQ(Pi,J)
  
  open <- matrix(NA,nrow=J,ncol=J+1)
  open[,1] <- 1:J
  closed <- matrix(NA,nrow=0,ncol=J+1)
  for(order in 1:J){open[order,J+1] <- totalcostheuristic_MB(Q,Pi,X,M,order)}
  
  curr_node <- na.exclude(open[which.min(open[,J+1]),1:J])
  continue <- TRUE
  while(continue){
    closed <- rbind(closed,open[which.min(open[,J+1]),])
    open <- open[-which.min(open[,J+1]),]
    
    next_objects <- setdiff(1:J,curr_node)
    open_new <- matrix(NA,nrow=length(next_objects),ncol=J+1)
    for(ind in 1:length(curr_node)){open_new[,ind] <- curr_node[ind]}
    open_new[,length(curr_node)+1] <- next_objects
    for(ind in 1:nrow(open_new)){open_new[ind,J+1] <- totalcostheuristic_MB(Q,Pi,X,M,na.exclude(open_new[ind,1:J]))}
    open <- rbind(open,open_new)
    curr_node <- na.exclude(open[which.min(open[,J+1]),1:J])
    if(length(curr_node)==J){continue <- FALSE}
  }
  if(any(is.na(open[,J+1]) | is.null(open[,J+1]))){stop("ERROR IN ASTAR COST CALCULATION")}
  
  phat <- phat_conditional(X,M,curr_node)
  thetahat <- theta_conditional(Pi,curr_node)
  return(list(pi0=curr_node,p=phat,
              theta=thetahat,num_nodes = nrow(open)+nrow(closed)))
}