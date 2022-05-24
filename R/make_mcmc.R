#' Convert Posteriors Samples from fit-btlbmm function into (Split) MCMC Object
#' 
#' This function samples approximately from the posterior of a Bradley-Terry-Luce-Binomial distribution with potentially K mixture components via an adaptive Metropolis-Hastings-within-Gibbs algorithm.
#' 
#' @import coda
#' 
#' @param gibbs Output of the fit_btlbmm function.
#' @param split Number of roughly-equal sized splits to perform of the posterior. Helpful if looking to calculate the Gelman diagnostic on the posterior. Default is no splitting (split=1).
#'  
#' @return MCMC list object.
#'  
#' @export
make_mcmc <- function(gibbs,split=1){
  K <- ncol(gibbs$alpha.out)
  J <- ncol(gibbs$p.out[,,1])
  mcmc <- gibbs$alpha.out
  for(k in 1:K){mcmc <- cbind(mcmc,gibbs$p.out[,,k])}
  mcmc <- cbind(mcmc,gibbs$theta.out)
  mcmc <- as.data.frame(mcmc)
  names(mcmc) <- c(paste0("alpha_",1:K),paste0(rep(paste0("p_",1:J,"_"),K),rep(1:K,each=J)),
                   paste0("theta_",1:K))
  if(split==1){return(as.mcmc(mcmc))
  }else{
    for(ind in 1:split){
      if(ind == 1){mcmc_list <- mcmc.list(as.mcmc(mcmc[ceiling(nrow(mcmc)/split*(ind-1)+1):ceiling(nrow(mcmc)/split*(ind)),]))
      }else{mcmc_list[[ind]] <- as.mcmc(mcmc[ceiling(nrow(mcmc)/split*(ind-1)+1):ceiling(nrow(mcmc)/split*(ind)),])}
    }
    return(mcmc_list)
  }
}