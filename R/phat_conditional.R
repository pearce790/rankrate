#' Estimate phat in a Mallows-Binomial given an order constraint
#'
#' This function calculates the MLE of p in a Mallows-Binomial(p,theta) distribution given Order(p).
#'
#' @import stats
#'
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param order Vector specifying a top-r or complete ordering of the desired p vector.
#' @param J Numeric specifying the total number of objects to be assessed.
#'
#' @return Vector of length J that is the MLE of p given order(p) and X.
#'
#' @examples
#' ratings <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' phat_conditional(ratings=ratings,M=5,order=c(1,2,3,4),J=4)
#' phat_conditional(ratings=ratings,M=5,order=c(2,1),J=4)
#'
#' @export
phat_conditional <- function(ratings,M,order,J){

  # Initial calculations
  R <- length(order)
  if(R<J){unordered <- setdiff(1:J,order)}else{unordered <- c()}
  order_and_unordered <- c(order,unordered)
  c1 <- apply(ratings[,order_and_unordered],2,function(ratings){sum(ratings,na.rm=T)})
  c2 <- apply(M-ratings[,order_and_unordered],2,function(ratings){sum(ratings,na.rm=T)})

  # Create Objective and Gradient Functions
  objective <- function(theta){
    return(-sum(c1*log(theta) + c2*log(1-theta)))
  }
  gradient <- function(theta){
    -c1/theta + c2/(1-theta)
  }

  # Create Linear Optimization Constraints
  ui <- rbind(diag(1,nrow=J),diag(-1,nrow=J))
  ci <- rep(c(0,-1),each=J)
  if(R>1){for(place in 2:R){
    ui_addition <- rep(0,J)
    ui_addition[c(place,place-1)] <- c(1,-1)
    ui <- rbind(ui,ui_addition)
    ci <- c(ci,0)
  }}
  if(R<J){for(place in (R+1):J){
    ui_addition <- rep(0,J)
    ui_addition[R] <- -1
    ui_addition[place] <- 1
    ui <- rbind(ui,ui_addition)
    ci <- c(ci,0)
  }}

  # Perform Constrained Optimization
  phat <- constrOptim(theta = seq(1/(J+1),J/(J+1),length=J),
                              f = objective,
                              grad = gradient,
                              outer.iterations = 1000,
                              ui=ui,ci=ci)

  phat$par <- (phat$par)[unlist(lapply(1:J,function(j){which(order_and_unordered==j)}))]
  return(list(phat = phat$par,objective = phat$value))
}
