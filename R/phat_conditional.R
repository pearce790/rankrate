#' Estimate phat in a Mallows-Binomial given an order constraint
#' 
#' This function calculates the MLE of p in a Mallows-Binomial(p,theta) distribution given Order(p).
#' 
#' @import nloptr
#' 
#' @param X Matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param order Vector specifying a top-r or complete ordering of the desired p vector.
#'  
#' @return Vector of length J that is the MLE of p given Order(p) and X.
#'  
#' @examples
#' X <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' phat_conditional(X=X,M=5,order=c(1,2,3,4))
#' phat_conditional(X=X,M=5,order=c(2,1))
#' 
#' @export
phat_conditional <- function(X,M,order){
  J <- ncol(X)
  
  opt_f <- function(p,X,M,order,J){
    # quantity that must be optimized with respect to p
    return(-sum(apply(X,1,function(x){dbinom(x,M,p,log=T)}),na.rm=T))
  }
  opt_g <- function(p,X,M,order,J){
    # contraints, which are met when all quantities returned are nonpositive
    return(c(-diff(p[order]),p[order[length(order)]]-p[setdiff(1:J,order)]))
  }
  res <- nloptr(x0=apply(X,2,function(x){mean(x,na.rm=T)})/M,
                eval_f = opt_f, eval_g_ineq = opt_g,
                lb = rep(0,J), ub = rep(1,J),
                X = X, M = M, order = order, J = J,
                opts = list("algorithm"="NLOPT_LN_COBYLA","xtol_abs"=1e-8,maxeval=1000))
  res$solution[which(res$solution==0)] <- 1e-8
  res$solution[which(res$solution==1)] <- 1-1e-8
  res$solution
}