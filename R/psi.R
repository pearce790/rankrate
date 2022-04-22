#' Psi function
#' 
#' This function calculates the normalizing constant of a Mallows distribution under the Kendall distance
#' 
#' @param theta Non-negative scale parameter.
#' @param J Positive integer indicating total number of objects.
#' @param R Positive integer <=\code{J} indicating size of partial ranking.
#' @param log Boolean indicating if log(Psi) should be returned.
#'  
#' @return Numeric representing normalizing constant of a Mallows distribution.
#'  
#' @examples 
#' psi(theta=1,J=10,R=8)
#' psi(theta=2,J=3,R=3,log=TRUE)
#'  
#' @export
psi <- function(theta,J,R,log=FALSE){
  # calculate the normalizing constant of a Mallows model (J = total number of items, R length of ranking)
  if(R>J){stop("R must be <= J")}
  if(R<=0 | J <= 0){stop("R and J must be positive integers")}
  if(theta<=0){stop("theta must be >=0")}
  
  log_psi <- sum(log((1-exp(-theta*(J-1:R+1))))-log((1-exp(-theta))))
  
  if(log){return(log_psi)
  }else{return(exp(log_psi))}
}