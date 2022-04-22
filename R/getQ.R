#' Calculate Q Matrix
#' 
#' This function calculates the Q matrix given a collection of (partial) rankings.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param J Total number of objects assessed.
#'  
#' @return Matrix Q of dimensions J x J.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' getQ(Pi=Pi,J=4)
#' 
#' @export
getQ <- function(Pi,J){
  ## Calculate the Q matrix given a collection of (partial) rankings
  ## Inputs: Pi = I x R matrix of partial rankings, J = total number of items
  Q <- matrix(NA,nrow=J,ncol=J)
  for(i in 1:J){for(j in 1:J){
    Q[i,j] <- mean(apply(Pi,1,function(pi){
      if(i %in% pi & j %in% pi){return(which(pi==i)<which(pi==j))}
      if(i %in% pi & !(j %in% pi)){return(TRUE)}
      if(!(i %in% pi) & j %in% pi){return(FALSE)}
      if(!(i %in% pi) & !(j %in% pi)){return(FALSE)}
    }))
  }}
  return(Q)
}