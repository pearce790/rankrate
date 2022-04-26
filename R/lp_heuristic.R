#' Calculate the LP heuristic of a Mallows-Binomial model
#' 
#' This function calculates the LP heuristic of a Mallows-Binomial model, for use during an A* tree search for the MLE of a Mallows-Binomial model.
#' 
#' @import lpSolve
#' 
#' @param Q Matrix of dimension J x J.
#' @param Pi Matrix of partial or complete rankings, one row per ranking.
#' @param I Numeric specifying number of judges
#' @param J Numeric specifying number of objects
#' @param order Vector specifying a top-r or complete ordering of the desired p vector.
#'  
#' @return Numeric specifying the LP heuristic.
lp_heuristic <- function(Q,Pi,I,J,order){
  idx <- function(i,j,J){J*(i-1)+j}
  
  coeffs <- Q*I
  S <- setdiff(1:J,order)
  coeffs <- coeffs[S,S]
  if(length(S)==1){
    pi0 <- c(order,S)
    dist <- sum(apply(Pi,1,function(rank){kendall(pi=rank,pi0=pi0)}))
    return(list(pi0=pi0,dist=dist))
  }
  if(length(S)>=2){
    pairs <- combinations(length(S),2)
    pairs_constraints <- matrix(0,nrow=length(S)*(length(S)-1)/2,ncol=length(S)^2)
    for(row in 1:choose(length(S),2)){
      items <- pairs[row,]
      pairs_constraints[row,c(idx(items[1],items[2],length(S)),
                              idx(items[2],items[1],length(S)))] <- 1
    }
  }else{pairs_constraints <- matrix(0,nrow=0,ncol=length(S)^2)}
  if(length(S)>=3){
    trios <- combinations(length(S),3)
    trios_constraints <- matrix(0,nrow=nrow(trios),ncol=length(S)^2)
    for(row in 1:choose(length(S),3)){
      items <- trios[row,]
      trios_constraints[row,c(idx(items[1],items[2],length(S)),
                              idx(items[2],items[3],length(S)),
                              idx(items[3],items[1],length(S)))] <- 1
    }
  }else{trios_constraints <- matrix(0,nrow=0,ncol=length(S)^2)}
  
  sol <- lp(direction = "min",
            objective = -as.vector(t(coeffs)),
            const.mat = rbind(pairs_constraints,trios_constraints),
            const.dir = c(rep("==",nrow(pairs_constraints)),
                          rep(">=",nrow(trios_constraints))),
            const.rhs = 1)
  vals <- rep(0,J)
  vals[S] <- apply(matrix(sol$solution,byrow=T,ncol=length(S)),2,sum)
  vals[order] <- seq(-length(order),-1,length=length(order))
  vals <- round(vals)
  if(length(vals)-length(unique(vals))>0){
    cond <- TRUE
    Pi0 <- unique(t(apply(matrix(1:50),1,function(i){
      rank <- rank(vals,ties.method="random")
      unlist(lapply(1:J,function(j){which(rank==j)}))
    })))
    while(cond){
      Pi0 <- unique(rbind(Pi0,t(apply(matrix(1:50),1,function(i){
        rank <- rank(vals,ties.method="random")
        unlist(lapply(1:J,function(j){which(rank==j)}))
      }))))
      if(nrow(Pi0)==prod(factorial(table(vals)))){cond<-FALSE}
    }
    dists <- apply(Pi0,1,function(pi0){sum(unlist(apply(Pi,1,function(rank){kendall(pi=rank,pi0=pi0)})))})
    return(list(pi0=Pi0[which.min(dists),],
                dist=dists[which.min(dists)]))
  }else{
    pi0 <- order(vals)
    dist <- sum(unlist(apply(Pi,1,function(rank){kendall(pi=rank,pi0=pi0)})))
    return(list(pi0=pi0,dist=dist))
  }
}