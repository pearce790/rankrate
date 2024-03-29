#' Toy data set of rankings and ratings demonstrating tie-breaking
#'
#' This toy data set includes 16 judges and 3 objects, and demonstrates the ability
#' of the Mallows-Binomial model to break ties in ratings via rankings.
#'
#' @format list with three elements: (1) \code{rankings}, a 16 x 3 matrix of rankings with one row per judge;
#'   (2) \code{ratings}, a 16 x 3 matrix of ratings, with one row per judge and one column per object; and
#'   (3) \code{M}, a number indicating the maximum (worst) integer score.
#' @source {Originally analyzed in: Gallo, Stephen A., et al. "A new approach to peer review assessments: Score, then rank" (2023). Research Integrity and Peer Review 8:10 (10). \url{https://researchintegrityjournal.biomedcentral.com/articles/10.1186/s41073-023-00131-7}.}
"ToyData1"
