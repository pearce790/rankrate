#' Random Mallows-Binomial generation.
#'
#' This function randomly generates rankings and ratings from a Mallows-Binomial distribution.
#'
#' @param I Numeric indicating the number of observations to be drawn, i.e., the number of judges providing rankings and ratings.
#' @param p Vector specifying the underlying object qualities.
#' @param pi0 Vector specifying the consensus (modal probability) ranking; should be used only for tie-breaking equal values in p.
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param M Numeric specifying the maximum integer rating.
#' @param R Numeric specifying the length of the (partial) rankings to be drawn.
#'
#' @return List containing elements X (I x J matrix of ratings) and Pi (I x R matrix of rankings).
#'
#' @examples
#' rmb(I=5,p=c(.1,.3,.4,.7,.9),theta=1,M=10)
#' rmb(I=10,p=c(.1,.3,.3,.7,.9),pi0=c(1,3,2,4,5),theta=5,M=40,R=3)
#'
#' @export

rmb <- function(I,p,theta,M,pi0=NULL,R=NULL){
  J <- length(p)
  if(is.null(pi0)){
    if(length(unique(p)) == J){pi0 <- order(p)
    }else{stop("Must specify unique p if pi0 is not supplied")}
  }

  ratings <- as.matrix(t(replicate(I,rbinom(n=J,size=M,prob=p))))
  rankings <- rmall(I,pi0,theta,R)

  return(list(ratings=ratings,rankings=rankings,M=M))
}
