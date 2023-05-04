#' Mallows-Binomial density function
#'
#' This function calculates the density of observation(s) under a Mallows-Binomial distribution.
#'
#' @import stats
#'
#' @param rankings A vector or matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One rankings per row.
#' @param ratings Vector or matrix of ratings, with one row per ratings for each judge.
#' @param p Vector of object qualities.
#' @param pi0 Vector specifying the consensus (modal probability) ranking (useful to break ties in p).
#' @param theta Numeric specifying the Mallows scale parameter.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param log Boolean indicating if loglikelihood should be returned.
#'
#' @return (Log) likelihood of rankings under a Mallows-Binomial distribution.
#'
#' @examples
#' rankings <- matrix(c(1,2,3,1,2,3),nrow=2,byrow=TRUE)
#' ratings <- matrix(c(2,4,6,4,2,6),nrow=2,byrow=TRUE)
#' dmb(rankings=c(1,2,3),ratings=c(2,4,6),p=c(.25,.5,.75),theta=2,M=8)
#' dmb(rankings=rankings,ratings=ratings,p=c(0.33,0.33,0.75),theta=6,M=8,pi0=c(1,2,3),log=TRUE)
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
  if(any(ratings<0 | ratings > M | ratings%%1 != 0)){stop("ratings must be integers between 0 and M")}
  if(is.null(pi0)){
    if(length(unique(p)) == length(p)){pi0 <- order(p)
    }else{stop("Must specify unique p if pi0 is not supplied")}
  }

  log_density <- dmall(rankings,pi0,theta,log=T)+
    sum(apply(ratings,1,function(rating){sum(dbinom(x=rating,size=M,prob=p,log=T),na.rm=T)}))
  if(log){return(log_density)}else{return(exp(log_density))}
}
