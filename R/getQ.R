#' Calculate Q Matrix
#'
#' This function calculates the Q matrix given a collection of (partial) rankings.
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate
#'   reviewer assignments. One ranking per row.
#' @param I A numeric entry indicating the total number of judges providing rankings and ratings.
#' @param J A numeric entry or vector of positive integers indicating total number of objects.
#'
#' @return A matrix with dimension \code{J} x \code{J}.
#'
#' @examples
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' getQ(rankings=rankings,I=2,J=4)
#' attr(rankings,"assignments") <- matrix(c(rep(TRUE,7),FALSE),byrow=TRUE,nrow=2,ncol=4)
#' getQ(rankings=rankings,I=2,J=4)
#'
#' @export
getQ <- function(rankings,I,J){
  ## Calculate the Q matrix given a collection of (partial) rankings
  ## Inputs: Pi = I x R matrix of partial rankings, J = total number of items
  Q <- matrix(NA,nrow=J,ncol=J)

  for(i in 1:J){for(j in 1:J){
    Q[i,j] <- mean(unlist(lapply(1:I,function(judge){
      ranking <- rankings[judge,]
      ballot <- attr(rankings,"assignments")[judge,]
      if(any(ballot[c(i,j)] == FALSE)){return(FALSE)
      }else{
        if(i %in% ranking){
          if(j %in% ranking){return(which(ranking == i) < which(ranking == j))
          }else{return(TRUE)}
        }else{
          if(j %in% ranking){return(FALSE)
          }else{return(FALSE)}
        }
      }
    })))
  }}
  return(Q)
}
