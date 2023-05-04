#' Kendall's tau function
#'
#' This function calculates Kendall's tau distance between ranking(s) and a central permutation, pi0
#'
#' @import stats
#'
#' @param rankings A vector or matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One rankings per row.
#' @param pi0 A vector containing a complete ranking
#'
#' @return Vector of the Kendall's tau distance between each ranking and pi0.
#'
#' @examples
#' ranking1 <- c(2,1,3)
#' ranking2 <- matrix(c(2,1,3,1,2,3),byrow=TRUE,nrow=2)
#' ranking3 <- matrix(c(1,2,3,4,2,4,NA,NA),byrow=TRUE,nrow=2)
#' attr(ranking3,"assignments") <- matrix(c(TRUE,TRUE,TRUE,TRUE,
#'   FALSE,TRUE,FALSE,TRUE),byrow=TRUE,nrow=2)
#' kendall(ranking1,c(1,2,3))
#' kendall(ranking2,c(1,2,3))
#' kendall(ranking3,c(1,2,3,4))
#'
#' @export
kendall <- function(rankings,pi0){

  # check data formats
  if(is.vector(rankings)){rankings <- matrix(rankings,nrow=1)}
  if(!is.matrix(rankings)){stop("rankings must be either a vector or matrix")}
  I <- nrow(rankings)
  J <- ncol(rankings)
  pi0 <- as.vector(pi0)
  if(is.null(attr(rankings,"assignments"))){attr(rankings,"assignments") <- matrix(TRUE,nrow=nrow(rankings),ncol=J)}

  # calculate distance for each ranking individually
  dists <- unlist(lapply(1:I,function(i){
    ranking <- na.exclude(rankings[i,])
    not_in_ballot <- which(attr(rankings,"assignments")[i,] == FALSE)
    R <- length(ranking)
    if(R==0){return(NA)}
    pi0_curr <- setdiff(pi0,not_in_ballot)
    J <- length(pi0_curr)

    if(R>J){stop("R must be <=J for all judges")}
    if(any(!(ranking %in% pi0_curr))){stop("No ranking can contain items not in their ballot")}
    dist <- 0
    for(r in 1:R){
      dist <- dist + (which(pi0_curr == ranking[r])-1)
      pi0_curr <- setdiff(pi0_curr,ranking[r])
    }
    return(dist)
  }))

  return(dists)
}
