#' Calculate the density of rankings under a Mallows distribution
#'
#' This function calculates the density of observation(s) under a Mallows distribution.
#'
#' @import stats
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate
#'   reviewer assignments. One ranking per row.
#' @param pi0 A vector specifying the consensus (modal probability) ranking; should be used only for tie-breaking
#'   equal values in \code{p}.
#' @param theta A numeric entry specifying the Mallows scale parameter.
#' @param log A boolean indicating if the log likelihood should be returned.
#'
#' @return A numeric value indicating the (log) likelihood of rankings under a Mallows distribution.
#'
#' @examples
#' rankings1 <- matrix(c(1,2,3,3,1,2),nrow=2,byrow=TRUE)
#' rankings2 <- matrix(c(1,2,3,4,2,3,NA,NA),nrow=2,byrow=TRUE)
#' attr(rankings2,"assignments") <- matrix(c(rep(TRUE,4),FALSE,TRUE,TRUE,TRUE),nrow=2,byrow=TRUE)
#' dmall(rankings=c(1,2,3,NA),pi0=c(1,2,3,4),theta=2)
#' dmall(rankings=rankings1,pi0=c(1,2,3),theta=2)
#' dmall(rankings=rankings2,pi0=c(1,2,3,4),theta=3,log=TRUE)
#'
#' @export
dmall <- function(rankings,pi0,theta,log=FALSE){
  if(is.vector(rankings)){rankings <- matrix(rankings,nrow=1)}
  if(!is.matrix(rankings)){stop("rankings must be either a vector or matrix")}
  I <- nrow(rankings)
  J <- length(pi0)
  pi0 <- as.vector(pi0)
  if(all(is.na(rankings))){return(NA)}
  if(is.null(attr(rankings,"assignments"))){attr(rankings,"assignments") <- matrix(TRUE,nrow=nrow(rankings),ncol=J)}

  Ji <- apply(attr(rankings,"assignments"),1,sum)
  Ri <- apply(rankings,1,function(ranking){sum(!is.na(ranking))})
  which_valid <- Ji > 0 & Ri > 0

  log_density <- -theta*sum(kendall(rankings,pi0),na.rm=T)-
    sum(psi(theta,J=Ji[which_valid],R=Ri[which_valid],log=T))

  if(log){return(log_density)}else{return(exp(log_density))}
}
