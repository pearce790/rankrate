#' Random Mallows-Binomial generation.
#' 
#' This function randomly generates rankings and ratings from a Mallows-Binomial distribution.
#' 
#' @param I Numeric indicating the number of observations to be drawn.
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
rmb <- function(I,p,pi0=NULL,theta,M,R=length(p)){
  if(is.null(pi0)){pi0 <- order(p)}
  J <- length(p)
  
  X <- as.matrix(t(replicate(I,rbinom(n=J,size=M,prob=p))))
  Pi <- matrix(rmall(I=I,pi0=pi0,theta=theta,R=R),nrow=I)
  
  return(list(X=X,Pi=Pi))
}