#' Calculate the density of rankings and ratings under a Mallows-Binomial distribution.
#'
#' This function calculates the density of observation(s) under a Mallows-Binomial distribution.
#'
#' @import stats
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate
#'   reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param p A vector specifying the underlying object qualities. All values between be between 0 and 1, inclusive.
#' @param pi0 A vector specifying the consensus (modal probability) ranking; should be used only for tie-breaking
#'   equal values in \code{p}.
#' @param theta A numeric entry specifying the Mallows scale parameter.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param log A boolean indicating if the log likelihood should be returned.
#'
#' @return A numeric value indicating the (log) likelihood of rankings and ratings under a Mallows distribution.
#'
#' @examples
#' data(ToyData1)
#' dmb(rankings=ToyData1$rankings,ratings=ToyData1$ratings,p=c(.2,.5,.7),theta=1,M=ToyData1$M)
#' dmb(rankings=ToyData1$rankings,ratings=ToyData1$ratings,p=c(.25,.25,.7),theta=1,M=ToyData1$M,
#' pi0=c(1,2,3),log=TRUE)
#'
#' @export
dmb <- function(rankings,ratings,p,theta,M,pi0=NULL,log=FALSE){

  if(is.vector(rankings)){rankings <- matrix(rankings,nrow=1)}
  if(is.vector(ratings)){ratings <- matrix(ratings,nrow=1)}
  if(!is.matrix(rankings)){stop("rankings must be either a vector or matrix")}
  if(!is.matrix(ratings)){stop("ratings must be either a vector or matrix")}

  I <- nrow(rankings)
  J <- length(p)
  if(any(dim(rankings) != c(I,J))){stop("rankings matrix must be dimension I x J")}
  if(any(dim(ratings) != c(I,J))){stop("ratings matrix must be dimension I x J")}
  if(any(ratings<0 | ratings > M | ratings%%1 != 0,na.rm=T)){stop("ratings must be integers between 0 and M")}
  if(is.null(pi0)){
    if(length(unique(p)) == length(p)){pi0 <- order(p)
    }else{stop("Must specify unique p if pi0 is not supplied")}
  }

  log_density <- dmall(rankings,pi0,theta,log=T)+
    sum(apply(ratings,1,function(rating){sum(dbinom(x=rating,size=M,prob=p,log=T),na.rm=T)}))
  if(log){return(log_density)}else{return(exp(log_density))}
}
