#' Calculate the LP total cost heuristic of a Mallows-Binomial model
#' 
#' This function calculates the LP total cost heuristic of a Mallows-Binomial model given Q, Pi, X, M, and an order, for use during an A* tree search for the MLE of a Mallows-Binomial model.
#' 
#' @param Q Matrix of dimension J x J.
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param order Vector specifying a top-r or complete ordering of the desired p vector.
#'  
#' @return Numeric specifying the LP total cost heuristic.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' Q <- matrix(c(0,.5,0,0,.5,0,0,0,1,1,0,0,1,1,.5,0),nrow=4,ncol=4)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' totalcostheuristic_MB_LP(Q=Q,Pi=Pi,X=X,M=5,order=c(1,2))
#' totalcostheuristic_MB_LP(Q=Q,Pi=Pi,X=X,M=5,order=c(2,1,4,3))
#' 
#' @export
totalcostheuristic_MB_LP <- function(Q,Pi,X,M,order){
  
  I <- max(nrow(X),nrow(Pi))
  J <- ncol(X)
  
  #calculate conditionally-optimal p
  phat <- phat_conditional(X=X,M=M,order=order)
  
  #calculate conditionally-optimal theta
  D <- suppressWarnings(lp_heuristic(Q=Q,Pi=Pi,I=I,J=J,order=order)$dist)/I
  thetahat <- theta_conditional(Pi=Pi,D=D,J=J)
  
  #calculate lp total cost heuristic
  R <- apply(Pi,1,function(pi){length(na.exclude(pi))})
  R[which(R == J-1)] <- J #missing one object implicitly defines a complete ranking
  return(thetahat*D + mean(unlist(lapply(R,function(r){psi(thetahat,J,r,log=T)})),na.rm=T)+
           -sum(apply(X,2,function(x){sum(x,na.rm=T)})*log(phat)+apply(M-X,2,function(x){sum(x,na.rm=T)})*log(1-phat))/nrow(Pi))
  
}