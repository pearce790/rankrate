#' Estimate the MLE of a Mallows-Binomial distribution using the Greedy method
#'
#' This function estimates the MLE of a Mallows-Binomial distribution using the Greedy method.
#'
#' @import gtools
#' @import isotone
#'
#' @param rankings A matrix of rankings, potentially with attribute "assignments" to signify separate reviewer assignments. One ranking per row.
#' @param ratings A matrix of ratings, one row per judge and one column per object.
#' @param M Numeric specifying maximum (=worst quality) integer rating.
#'
#' @return List with elements pi0 (estimated consensus ranking MLE),  p (estimated object quality parameter MLE), theta (estimated scale parameter MLE), and numnodes (number of nodes traversed during algorithm, a measure of computational complexity).
#'
#' @examples
#' data("ToyData1")
#' Greedy(ToyData1$rankings,ToyData1$ratings,ToyData1$M)
#'
#' @export
Greedy <- function(rankings,ratings,M){

  # calculate constants
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings) != c(I,J))){stop("rankings and ratings must have same dimension")}
  if(is.null(attr(rankings,"assignments"))){
    attr(rankings,"assignments") <- matrix(TRUE,nrow=I,ncol=J)
  }

  # calculations for updating theta
  Ji <- apply(attr(rankings,"assignments"),1,sum) #each judge's total number of assignments
  Ri <- apply(rankings,1,function(pi){sum(!is.na(pi))}) #each judge's size of ranking
  i_pi <- which(Ri>0 & Ji>0) #which judges provided rankings
  t1 <- function(theta,D){theta*D/I} #first log probability term
  t2 <- function(theta){sum(psi(theta,Ji[i_pi],Ri[i_pi],log=T))/I} #second log probability term
  get_theta <- function(D){ # function used to optimize for theta conditional on some order
    # D is the minimum summed kendall distance between some central ordering and the observations.
    # As D increases, theta decreases monotonically and the cost increases monotonically.
    opt_theta <- function(theta){t1(theta,D=D)+t2(theta)}
    result <- optimize(f=opt_theta,interval=c(1e-8,10^8))
    return(list(thetahat = result$minimum,objective = result$objective))
  }

  # calculations for updating p
  c1 <- apply(ratings,2,function(x){sum(x,na.rm=T)})/I #set of constants in log binomial density
  c2 <- apply(M-ratings,2,function(x){sum(x,na.rm=T)})/I #set of constants in log binomial density
  t3 <- function(p){sum(c1*log(1/p)+c2*log(1/(1-p)))} #calculation of log binomial density given p
  t3_gradient <- function(p){ -c1/p + c2/(1-p) } # gradient of t3 function, for use in optimization later on
  get_p <- function(order){ # constrained optimization function
    # order is the top-R portion of a central ordering (i.e., if order = c(5,4), then object 5 must be in first place
    # and object 4 must be in second place, the remaining objects must be in places 3:J (any order among those) ).

    # initializing the order
    R <- length(order)
    if(R<J){unordered <- setdiff(1:J,order)}else{unordered <- c()}
    order_and_unordered <- c(order,unordered)

    # calculation of linear constraints
    ui <- rbind(diag(1,nrow=J),diag(-1,nrow=J))
    ci <- rep(c(0,-1),each=J)
    if(R>1){for(place in 2:R){
      ui_addition <- rep(0,J)
      ui_addition[c(order[place],order[place-1])] <- c(1,-1)
      ui <- rbind(ui,ui_addition)
      ci <- c(ci,0)
    }}
    if(R<J){for(place in (R+1):J){
      ui_addition <- rep(0,J)
      ui_addition[c(order_and_unordered[place],order[R])] <- c(1,-1)
      ui <- rbind(ui,ui_addition)
      ci <- c(ci,0)
    }}

    # get intelligent starting values via isotonic regression
    meanratings <- c1/(c1+c2)
    order_meanratings <- c(order,setdiff(order(meanratings),order))
    start <- gpava(1:J,meanratings[order_meanratings])
    start <- start$x[order(order_meanratings)]

    # optimize for p conditional on an ordering of parameters
    result <- constrOptim(theta = start,f = t3,grad = t3_gradient,ui=ui,ci=ci-1e-8,
                          mu=.1,control = list(reltol = 1e-10))
    return(list(phat = result$par,objective = result$value))
  }

  # get D (necessary for updating theta)
  Q <- getQ(rankings,I,J)*I
  get_D <- function(order){
    S <- setdiff(1:J,order)
    D <- 0
    if(length(S)>=2){
      D <- D +
        sum(apply(combinations(length(S),2,S),1,function(uv){
          min(Q[uv[1],uv[2]],Q[uv[2],uv[1]])
        }))
    }
    for(i in 1:length(order)){D <- D + sum(Q[setdiff(1:J,order[1:i]),order[i]])}
    return(D)
  }

  # START ESTIMATION VIA GREEDY ALGORITHM

  num_nodes <- 0
  result <- matrix(NA,nrow=1,ncol=0)
  while(ncol(result)<J){
    tmp <- matrix(NA,nrow=0,ncol=ncol(result)+1)
    for(node_index in 1:nrow(result)){
      current <- result[node_index,]
      neighbors <- setdiff(1:J,current)
      cost <- unlist(lapply(neighbors,function(neighbor){
        order <- c(current,neighbor)
        round(get_theta(D=get_D(order))$objective + get_p(order)$objective,4)
      }))
      num_nodes <- num_nodes + length(neighbors)
      whichmin <- which(cost==min(cost))
      for(which in whichmin){tmp <- rbind(tmp,c(current,neighbors[which]))}
    }
    result <- tmp
  }

  if(nrow(result)>1){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    results <- t(apply(result,1,function(curr_node){
      c(get_p(curr_node)$phat,get_theta(D=get_D(curr_node))$thetahat)
    }))
    return(list(pi0=result,
                p=results[,1:J],
                theta=results[,J+1,drop=FALSE],
                num_nodes=num_nodes))
  }else{
    curr_node <- result[1,]
    results <- c(get_p(curr_node)$phat,get_theta(D=get_D(curr_node))$thetahat)
    return(list(pi0=curr_node,
                p=results[1:J],
                theta=results[J+1],
                num_nodes=num_nodes))
  }
}
