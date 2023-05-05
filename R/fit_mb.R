#' Calculate the exact or approximate MLE of a Mallows-Binomial distribution using various methods
#'
#' This function calculates the exact or approximate MLE of a Mallows-Binomial distribution using a user-specified method.
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param method String indicating which estimation method to use when estimating parameters. Allowable options are currently "ASTAR", "Greedy", "GreedyLocal", and "FV". Defaults to exact search, "ASTAR".
#'
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'
#' @examples
#' data("ToyData1")
#' fit_mb(ToyData1$rankings,ToyData1$ratings,ToyData1$M,method="ASTAR")
#' fit_mb(ToyData1$rankings,ToyData1$ratings,ToyData1$M,method="Greedy")
#' fit_mb(ToyData1$rankings,ToyData1$ratings,ToyData1$M,method="GreedyLocal")
#' fit_mb(ToyData1$rankings,ToyData1$ratings,ToyData1$M,method="FV")
#'
#' @export
fit_mb <- function(rankings,ratings,M,method=c("ASTAR","Greedy","GreedyLocal","FV")){
  if(method == "ASTAR"){
    return(ASTAR(rankings,ratings,M))
  }else if(method == "Greedy"){
    return(Greedy(rankings,ratings,M))
  }else if(method == "GreedyLocal"){
    return(GreedyLocal(rankings,ratings,M))
  }else if(method == "FV"){
    return(FV(rankings,ratings,M))
  }else{stop("Need valid estimation method")}
}
