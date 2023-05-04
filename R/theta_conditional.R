#' Estimate theta in a Mallows-Binomial distribution under a constraint.
#'
#' This function calculates the MLE of the consensus scale parameter, theta, in a Mallows-Binomial(p,theta) distribution given a constraint.
#'
#' @import stats
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param D Numeric specifying a minimum total kendall distance between the observed rankings and a supposed central ranking.
#' @param J Numeric specifying the total number of objects to be assessed.
#'
#' @return Numeric specifying the MLE of theta under a constraint.
#'
#' @examples
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' theta_conditional(rankings,D=1,J=4)
#' theta_conditional(rankings,D=2,J=4)
#'
#' @export
theta_conditional <- function(rankings,D,J){

  I <- nrow(rankings)
  if(is.null(attr(rankings,"assignments"))){attr(rankings,"assignments") <- matrix(TRUE,nrow=I,ncol=J)}

  Ji <- apply(attr(rankings,"assignments"),1,sum)
  Ri <- apply(rankings,1,function(ranking){sum(!is.na(ranking))})
  which_valid <- Ji > 0 & Ri > 0

  objective <- function(theta){
    theta*D+sum(psi(theta,J=Ji[which_valid],R=Ri[which_valid],log=T))/I
  }
  theta_hat <- optimize(objective,c(1e-8,1e8))

  return(list(thetahat = theta_hat$minimum,objective = theta_hat$objective))
}
