#' Random Mallows generation.
#' 
#' This function randomly generates rankings from a Mallows distribution.
#' 
#' @param I Numeric indicating the number of observations to be drawn.
#' @param pi0 Vector specifying the consensus (modal probability) ranking.
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param R Numeric specifying the length of the (partial) rankings to be drawn.
#'  
#' @return Matrix of rankings, one row per ranking.
#'  
#' @examples
#' rmall(I=5,pi0=1:5,theta=1,R=3)
#' rmall(I=10,pi0=1:5,theta=1)
#'  
#' @export
rmall <- function(I,pi0,theta,R=length(pi0)){
  tmp <- function(pi0,theta){
    J <- length(pi0)
    Vj <- pi1 <- c()
    for(j in 1:(J-1)){
      probs <- exp(-theta*(0:(J-j)))/sum(exp(-theta*(0:(J-j))))
      Vj <- c(Vj,sample.int(J-j+1,size=1,prob=probs)-1)
    }
    Vj <- c(Vj,0)
    for(j in 1:J){
      pi1 <- c(pi1,pi0[Vj[j]+1])
      pi0 <- setdiff(pi0,pi0[Vj[j]+1])
    }
    return(pi1)
  }
  t(replicate(I,tmp(pi0,theta)))[,1:R]
}