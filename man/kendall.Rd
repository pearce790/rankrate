% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kendall.R
\name{kendall}
\alias{kendall}
\title{Calculate the Kendall's tau between rankings}
\usage{
kendall(rankings, pi0)
}
\arguments{
\item{rankings}{A matrix of rankings, potentially with attribute "assignments" to signify separate
reviewer assignments. One ranking per row.}

\item{pi0}{A vector specifying the consensus (modal probability) ranking.}
}
\value{
A vector of the Kendall's tau distance between each ranking in \code{rankings} and \code{pi0}.
}
\description{
This function calculates Kendall's tau distance between ranking(s) and a central permutation, pi0
}
\examples{
ranking1 <- c(2,1,3)
ranking2 <- matrix(c(2,1,3,1,2,3),byrow=TRUE,nrow=2)
ranking3 <- matrix(c(1,2,3,4,2,4,NA,NA),byrow=TRUE,nrow=2)
attr(ranking3,"assignments") <- matrix(c(TRUE,TRUE,TRUE,TRUE,
  FALSE,TRUE,FALSE,TRUE),byrow=TRUE,nrow=2)
kendall(ranking1,c(1,2,3))
kendall(ranking2,c(1,2,3))
kendall(ranking3,c(1,2,3,4))

}
