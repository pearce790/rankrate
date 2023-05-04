#' Calculate Confidence Intervals for Mallows-Binomial parameters.
#'
#' This function calculates confidence intervals for parameters in a Mallows-Binomial model using the nonparametric bootstrap.
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#' @param interval Numeric between 0 and 1 specifying the confidence interval (e.g., .90 indicates a 90% confidence interval). Defaults to 0.90.
#' @param nsamples Numeric indicating desired number of bootstrap samples to be used when calculating confidence intervals. Defaults to 100.
#' @param all Boolean indicating if estimated parameters from all bootstrap samples should be returned. Defaults to FALSE.
#' @param method String indicating which estimation method to use when estimating parameters. Allowable options are currently "ASTAR", "Greedy", "GreedyLocal", and "FV". Defaults to exact search, "ASTAR".
#'
#' @return List with elements ci (matrix of confidence intervals for Mallows-Binomial parameters), ci_ranks (matrix of confidence intervals for object ranks), bootstrap_pi0 (matrix of bootstrap consensus rankings; returned only if all=TRUE), and bootstrap_ptheta (matrix of bootstrap estimates of (p,theta); returned only if all=TRUE).
#'
#' @examples
#' rankings <- matrix(c(1,2,3,4,2,1,NA,NA),byrow=TRUE,nrow=2)
#' ratings <- matrix(c(0,1,2,3,1,2,2,5),byrow=TRUE,nrow=2)
#' ci_mb(rankings=rankings,ratings=ratings,M=5,method="FV",all=TRUE)
#'
#' @export
ci_mb <- function(rankings,ratings,M,interval=0.90,nsamples=10,all=FALSE,method="ASTAR"){
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  if(all){bs_pi0 <- matrix(NA,nrow=nsamples,ncol=J)}
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
    if(all){bs_pi0[sample,] <- bs_res$pi0}
    bs_parameters[sample,] <- c(bs_res$p,bs_res$theta)
  }
  ci <- as.matrix(apply(bs_parameters,2,function(parameter){quantile(parameter,probs=c((1-interval)/2,1-(1-interval)/2))}))
  colnames(ci) <- colnames(bs_parameters) <- c(paste0("p",1:J),"theta")
  ci_ranks <- matrix(unlist(lapply(1:J,function(j){
    quantile(apply(bs_pi0,1,function(pi0){which(pi0==j)}),probs=c((1-interval)/2,1-(1-interval)/2))
  })),nrow=2)
  rownames(ci_ranks) <- rownames(ci)
  colnames(ci_ranks) <- paste0("Item",1:J)

  if(all){return(list(ci=ci,ci_ranks=ci_ranks,
                      bootstrap_pi0=bs_pi0,
                      bootstrap_ptheta=bs_parameters))
  }else{return(list(ci=ci,ci_ranks=ci_ranks))}
}
