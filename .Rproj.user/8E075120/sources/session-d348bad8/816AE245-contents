#' Estimate thetahat in a Mallows-Binomial given a constraint
#' 
#' This function calculates the MLE of theta in a Mallows-Binomial(p,theta) distribution given a constraint.
#' 
#' @import stats
#' 
#' @param Pi Matrix of rankings, with one row per ranking
#' @param order Vector specifying a complete ordering of the desired p vector. Either order or D must be specified; if both, use the order.
#' @param D Numeric specifying a minimum total kendall distance between the observed rankings and a supposed central ranking. Either order or D must be specified; if both, use the order.
#' @param J Numeric specifying the total number of objects to be assessed. Not used if order is specified.
#'  
#' @return Numeric specifying the MLE of theta under a constraint.
#'  
#' @examples
#' Pi <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' theta_conditional(Pi=Pi,order=c(1,2,4,3))
#' theta_conditional(Pi=Pi,D=3,J=4)
#' 
#' @export
theta_conditional <- function(Pi,order=NULL,D=NULL,J=NULL){
  # optimize for MLE of scale parameter
  
  if(!is.null(order)){
    return(optim(par=1,fn=function(theta){-dmall(Pi,order,theta,log=T)},
                 lower=1e-5,upper=Inf,method="L-BFGS-B")$par)
  }else if(!is.null(D)){
    if(D<0){stop("D cannot be less than 0")}
    if(is.null(J)){stop("If D specified, then J must be specified as well.")}
    R <- apply(Pi,1,function(pi){length(na.exclude(pi))})
    R[which(R == J-1)] <- J #missing one object implicitly defines a complete ranking
    
    return(optim(1,function(theta){theta*D + mean(unlist(lapply(R,function(r){psi(theta,J,r,log=T)})),na.rm=T)},
                 lower=1e-5,upper=Inf,method="L-BFGS-B")$par)
  }else{
    return(NULL)
  }
}