#' Normalizing constant function of a Mallows distribution, psi
#'
#' This function calculates the normalizing constant of a Mallows distribution under the Kendall distance
#'
#' @param theta A numeric entry specifying the Mallows scale parameter.
#' @param J A numeric entry or vector of positive integers indicating total number of objects each judge has access to.
#'   If \code{length(J)>1}, \code{R} must be of same length or a single value.
#' @param R A numeric entry or vector of positive integers indicating the length of the ranking provided by each judge.
#'   If \code{length(R)>1}, \code{J} must be of same length or a single value.
#' @param log A boolean indicating if \code{log(psi)} should be returned.
#'
#' @return A numeric value or vector representing normalizing constant of a Mallows distribution.
#'
#' @examples
#' psi(theta=1,J=10,R=8)
#' psi(theta=2,J=c(4,4,3),R=c(2,2,1),log=TRUE)
#'
#' @export
psi <- function(theta,J,R,log=FALSE){
  if(length(theta)>1){stop("theta must be a single numeric value")}
  if(length(R)>1 & length(J)==1){J <- rep(J,length(R))}
  if(length(J)>1 & length(R)==1){R <- rep(R,length(J))}
  if(length(R)!=length(J)){stop("If length(R)>1 and length(J)>1, length(R) must equal length(J)")}
  if(any(R>J)){stop("R must be <= J")}
  if(any(R<=0) | any(J <= 0)){stop("R and J must be positive integers")}
  if(theta<=0){stop("theta must be >=0")}

  log_psi <- unlist(lapply(1:length(R),function(index){
    Rcurr <- R[index]
    Jcurr <- J[index]
    sum(log((1-exp(-theta*(Jcurr-1:Rcurr+1))))-log((1-exp(-theta))))
  }))

  if(log){return(log_psi)
  }else{return(exp(log_psi))}
}
