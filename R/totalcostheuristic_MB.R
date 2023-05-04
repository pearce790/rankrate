#' Calculate the naive total cost heuristic of a Mallows-Binomial model
#'
#' This function calculates the total cost heuristic of a Mallows-Binomial model given Q, rankings, ratings, M, and an order, for use during an A* tree search for the MLE of a Mallows-Binomial model.
#'
#' @import gtools
#'
#' @param Q Matrix of dimension J x J.
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param order Vector specifying a top-r or complete ordering of the desired p vector.
#'
#' @return List of three numeric values: estimated object quality vector (phat), estimated consensus scale parameter (thetahat), and the total cost heuristic (totalcostheuristic)
#'
#' @examples
#' Q <- matrix(c(0,.5,0,0,.5,0,0,0,1,1,0,0,1,1,.5,0),nrow=4,ncol=4)
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' ratings <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' totalcostheuristic_MB(Q=Q,rankings=rankings,ratings=ratings,M=5,order=c(1,2))
#' totalcostheuristic_MB(Q=Q,rankings=rankings,ratings=ratings,M=5,order=c(2,1,4,3))
#'
#' @export
totalcostheuristic_MB <- function(Q,rankings,ratings,M,order){

  J <- ncol(Q)

  # calculate conditional phat
  phat <- phat_conditional(ratings,M,order,J)

  #calculate conditional theta
  S <- setdiff(1:J,order)
  D <- 0
  if(length(S)>=2){D <- D + sum(apply(combinations(length(S),2,S),1,function(uv){
    min(Q[uv[1],uv[2]],Q[uv[2],uv[1]])}))}
  for(i in 1:length(order)){D <- D + sum(Q[setdiff(1:J,order[1:i]),order[i]])}
  thetahat <- theta_conditional(rankings,D,J)

  return(list(phat = phat$phat,thetahat = thetahat$thetahat,
              totalcostheuristic = phat$objective + thetahat$objective))
}
