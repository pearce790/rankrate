#' Kendall's tau function
#' 
#' This function calculates Kendall's tau distance between a (partial) ranking pi and a complete ranking pi0.
#' 
#' @import stats
#' 
#' @param pi A partial or complete ranking.
#' @param pi0 A complete ranking.
#'  
#' @return Numeric Kendall's tau distance between pi and pi0.
#'  
#' @examples 
#' kendall(pi=c(1,2,3),pi0=c(4,3,2,1))
#' kendall(pi=c(1,2,3),pi0=c(1,3,2))
#'  
#' @export
kendall <- function(pi,pi0){
  # calculate kendall distance between a ranking pi and a central ranking, pi0
  pi <- na.exclude(as.vector(pi))
  pi0 <- as.vector(pi0)
  
  R <- length(pi)
  J <- length(pi0)
  if(length(setdiff(1:J,pi0))!=0){stop("pi0 must be a complete ranking")}
  if(any(pi>J,na.rm=T)){stop("pi cannot contain items not in pi0")}
  if(R>J){stop("R must be <= J")}
  dist <- 0
  for(r in 1:R){
    dist <- dist + (which(pi0 == pi[r]) - 1)
    pi0 <- setdiff(pi0,pi[r])
  }
  return(dist)
}