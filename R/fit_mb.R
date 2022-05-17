#' Calculate the exact or approximate MLE of a Mallows-Binomial distribution using various methods
#' 
#' This function calculates the exact or approximate MLE of a Mallows-Binomial distribution using a user-specified method.
#' 
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param method String specifying method, with allowable options "ASTAR","ASTAR_LP","Greedy","GreedyLocal",and "FV". 
#' @param localsearch Numeric for use with the Greedy and FV methods; see documentation of those estimation functions for details. Defaults to 0 indicating no local search.
#' 
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' fit_mb(Pi=Pi,X=X,M=5,method="ASTAR")
#' fit_mb(Pi=Pi,X=X,M=5,method="Greedy",localsearch=0)
#' fit_mb(Pi=Pi,X=X,M=5,method="GreedyLocal")
#' fit_mb(Pi=Pi,X=X,M=5,method="FV",localsearch=3)
#'  
#' @export
fit_mb <- function(Pi,X,M,method=c("ASTAR","ASTAR_LP","Greedy","GreedyLocal","FV"),localsearch=0){
  if(method == "ASTAR"){
    return(ASTAR(Pi,X,M))
  }else if(method == "ASTAR_LP"){
    return(ASTAR_LP(Pi,X,M))
  }else if(method == "Greedy"){
    return(Greedy(Pi,X,M,localsearch))
  }else if(method == "GreedyLocal"){
    return(GreedyLocal(Pi,X,M))
  }else if(method == "FV"){
    return(FV(Pi,X,M,localsearch))
  }else{stop("Need valid estimation method")}
}