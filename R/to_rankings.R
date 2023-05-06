#' Convert ranks into rankings (orderings)
#'
#' This function converts a matrix of ranks into a matrix of rankings (i.e., orderings), potentially including reviewer
#' assignments as an attribute of the ranking matrix. Additionally, it can be used to add an assignments matrix to an existing
#' matrix of rankings.
#'
#' @param ranks A matrix or vector of ranks, such that the (i,j) entry includes the rank given by judge i to proposal j.
#'   \code{NA} is used to indicate that no rank was assigned to a proposal, which may occur for two reasons: (1)
#'   If the \code{assignments} matrix is not specified or the (i,j) entry of assignments is \code{TRUE}, then an \code{NA} indicates that
#'   a proposal was considered worse than all ranked proposals. (2) If the (i,j) entry of \code{assignments} is \code{FALSE}, then \code{NA}
#'   indicates that a proposal was not considered by the judge and no information can be gleaned from the missing rank.
#' @param assignments A matrix of booleans, such that the (i,j) entry is \code{TRUE} if judge i was assigned to review proposal j, and
#'   \code{FALSE} otherwise. If \code{assignments} is \code{NULL}, we assume all judges considered all proposals.
#' @param rankings A matrix or vector of rankings. If a matrix, there should be one ranking per row.
#'
#' @return A matrix of rankings, with one row per ranking. If \code{assignments} argument is specified, then the rankings matrix will have
#'   the attribute "assignments".
#'
#' @examples
#' ranks <- matrix(data=c(4,2,3,1,NA,1,2,3,NA,NA,1,NA),byrow=TRUE,nrow=3)
#' assignments=matrix(TRUE,byrow=TRUE,nrow=3,ncol=4)
#' to_rankings(ranks=ranks)
#' to_rankings(ranks=ranks,assignments=assignments)
#' to_rankings(assignments=matrix(TRUE,nrow=1,ncol=3),rankings=c(3,2,1))
#'
#' @export
to_rankings <- function(ranks,assignments=NULL,rankings=NULL){

  if(!is.null(rankings)){
    if(is.vector(rankings)){rankings <- matrix(rankings,nrow=1)}
    if(!is.matrix(rankings)){stop("rankings must be a vector or matrix")}
    I <- nrow(rankings)
    J <- ncol(rankings)

    if(!is.null(assignments)){
      if(!is.matrix(assignments)){stop("assignments must be a matrix of dimension IxJ")}
      if(nrow(assignments)!=I | ncol(assignments)!=J){stop("assignments must be a matrix of dimensions IxJ")}
      if(any(is.na(assignments) | is.null(assignments))){stop("assignments must only contain boolean values TRUE/FALSE")}
      attr(rankings,"assignments") <- assignments
    }
  }else{
    if(is.vector(ranks)){ranks <- matrix(ranks,nrow=1)}
    I <- nrow(ranks)
    J <- ncol(ranks)
    rankings <- matrix(NA,nrow=I,ncol=J)

    for(i in 1:I){
      curr_ranks <- ranks[i,]
      if(all(is.na(curr_ranks))){ # when no rank data for judge i, rankings is left NA
      }else{
        Ri <- sum(!is.na(curr_ranks))
        rankings[i,1:Ri] <- order(curr_ranks)[1:Ri]
      }
    }

    if(!is.null(assignments)){
      if(!is.matrix(assignments)){
        stop("assignments must be a matrix of dimension IxJ")
      }
      if(nrow(assignments)!=I | ncol(assignments)!=J){
        stop("assignments must be a matrix of dimensions IxJ")
      }
      if(any(is.na(assignments) | is.null(assignments))){
        stop("assignments must only contain boolean values TRUE/FALSE")
      }
      attr(rankings,"assignments") <- assignments
    }
  }

  return(rankings)
}
