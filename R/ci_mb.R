#' Bootstrap Confidence Intervals for Mallows-Binomial parameters.
#'
#' This function calculates confidence intervals for parameters in a Mallows-Binomial model
#' using the nonparametric bootstrap.
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate
#'   reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param interval A numeric entry between 0 and 1 specifying the confidence interval (e.g.,
#'   .90 indicates a 90% confidence interval). Defaults to 0.90.
#' @param nsamples A numeric entry indicating desired number of bootstrap samples to be used when
#'   calculating confidence intervals. Defaults to 50.
#' @param all A boolean indicating if estimated parameters from all bootstrap samples should be returned.
#'   Defaults to \code{FALSE}.
#' @param method A character string indicating which estimation method to use when estimating parameters.
#'   Allowable options are currently "ASTAR", "Greedy", "GreedyLocal", and "FV". Defaults to exact search, "ASTAR".
#'
#' @return A list with elements \code{ci}, a matrix of confidence intervals for Mallows-Binomial parameters,
#'   \code{ci_ranks}, a matrix of confidence intervals for object ranks, \code{bootstrap_pi0}, a matrix of
#'   bootstrap consensus rankings (returned only if \code{all==TRUE}), and \code{bootstrap_ptheta}, a
#'   matrix of bootstrap estimates of (p,theta) (returned only if \code{all==TRUE}).
#'
#' @examples
#' data("ToyData1")
#' ci_mb(ToyData1$rankings,ToyData1$ratings,ToyData1$M,method="ASTAR",all=TRUE)
#'
#' @export
ci_mb <- function(rankings,ratings,M,interval=0.90,nsamples=50,all=FALSE,method="ASTAR"){
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  bs_pi0 <- matrix(NA,nrow=nsamples,ncol=J)
  bs_parameters <- matrix(NA,nrow=nsamples,ncol=J+1)

  for(sample in 1:nsamples){
    bs_sample <- sample(1:I,I,replace=T)
    bs_res <- suppressMessages(fit_mb(rankings[bs_sample,],ratings[bs_sample,],M,method))
    if(is.matrix(bs_res$pi0)){
      which_keep <- sample(1:nrow(bs_res$pi0),1)
      bs_res$pi0 <- bs_res$pi0[which_keep,]
      bs_res$p <- bs_res$p[which_keep,]
      bs_res$theta <- bs_res$theta[which_keep,]
    }
    bs_pi0[sample,] <- bs_res$pi0
    bs_parameters[sample,] <- c(bs_res$p,bs_res$theta)
  }
  ci <- as.matrix(apply(bs_parameters,2,function(parameter){quantile(parameter,probs=c((1-interval)/2,1-(1-interval)/2))}))
  colnames(ci) <- colnames(bs_parameters) <- c(paste0("p",1:J),"theta")
  ci_ranks <- matrix(unlist(lapply(1:J,function(j){
    quantile(apply(bs_pi0,1,function(pi0){which(pi0==j)}),probs=c((1-interval)/2,1-(1-interval)/2))
  })),nrow=2)
  rownames(ci_ranks) <- rownames(ci)
  colnames(ci_ranks) <- paste0("Object",1:J)

  if(all){return(list(ci=ci,ci_ranks=ci_ranks,
                      bootstrap_pi0=bs_pi0,
                      bootstrap_ptheta=bs_parameters))
  }else{return(list(ci=ci,ci_ranks=ci_ranks))}
}
