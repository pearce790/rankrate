#' Convert ranks into rankings (orderings)
#'
#' This function converts a matrix of ranks into a matrix of orderings, potentially including reviewer assignments as an attribute of the ranking matrix.
#'
#' @param ranks Matrix or vector of ranks, such that the (i,j) entry includes the rank given by judge i to proposal j. NA used to indicate no rank was assigned to a proposal, which may occur for two reasons: (1) If "assignments" is not specified or the (i,j) entry of assignments is TRUE, NA indicates a proposal was considered worse than all ranked proposals. If the (i,j) entry of assignments is FALSE, NA indicates a proposal simply was not considered by the judge and no information can be gleaned from the missing rank.
#' @param assignments Matrix of booleans, such that the (i,j) entry is TRUE if judge i was assigned to review proposal j, and FALSE otherwise. If assignments is NULL, we assume all judges considered all proposals.
#' @param rankings Matrix or vector of rankings, with one row per judge. If rankings is specified, the function is used only to add the assignments attribute to the rankings file; any ranks input is ignored.
#'
#' @return Matrix of rankings, one row per ranking. If "assignments" argument is specified, then the rankings matrix will have the attribute "assignments".
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
