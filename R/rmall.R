#' Random Mallows generation.
#'
#' This function randomly generates rankings from a Mallows distribution.
#'
#' @param I Numeric indicating the number of observations to be drawn.
#' @param pi0 Vector specifying the consensus (modal probability) ranking.
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param R Numeric specifying the length of the (partial) rankings to be drawn. If NULL, R = length(pi0), i.e., complete rankings are drawn.
#'
#' @return Matrix of rankings, one row per ranking.
#'
#' @examples
#' rmall(I=5,pi0=1:5,theta=1,R=3)
#' rmall(I=5,pi0=1:3,theta=.5,R=c(1,1,1,1,3))
#' rmall(I=5,pi0=1:3,theta=.5)
#'
#' @export
rmall <- function(I,pi0,theta,R=NULL){

  J <- length(pi0)
  if(length(theta)>1){stop("theta must be a single numeric value")}
  if(is.null(R)){R <- rep(J,I)}
  if(length(R)==1){R <- rep(R,I)}
  if(length(R)!=I){stop("R must be of length 1 or I")}
  if(any(R>J)){stop("R cannot be larger than J")}
  if(any(R == (J-1))){
    R[R== (J-1)] <- J
    warning("When R == J-1, automatically set R = J")
  }

  rankings <- matrix(NA,nrow=I,ncol=J)
  for(i in 1:I){
    pi0_curr <- pi0
    Vj <- pi1 <- c()
    for(j in 1:(J-1)){
      probs <- exp(-theta*(0:(J-j)))/sum(exp(-theta*(0:(J-j))))
      Vj <- c(Vj,sample.int(J-j+1,size=1,prob=probs)-1)
    }
    Vj <- c(Vj,0)
    for(j in 1:J){
      pi1 <- c(pi1,pi0_curr[Vj[j]+1])
      pi0_curr <- setdiff(pi0_curr,pi0_curr[Vj[j]+1])
    }
    rankings[i,1:R[i]] <- pi1[1:R[i]]
  }
  return(rankings)
}

