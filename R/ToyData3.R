#' Toy data set of rankings and ratings when judges express internally inconsistent preferences
#'
#' This toy data set includes 16 judges and 3 objects, and demonstrates the ability
#' of the Mallows-Binomial model to estimate overall object orderings even when judges provide
#' sets of rankings and ratings which are internally inconsistent.
#'
#' @format list with three elements: (1) \code{rankings}, a 16 x 3 matrix of rankings with one row per judge;
#'   (2) \code{ratings}, a 16 x 3 matrix of ratings, with one row per judge and one column per object; and
#'   (3) \code{M}, a number indicating the maximum (worst) integer score.
#' @source {Originally published in: Gallo, Stephen A., et al. "A new approach to peer review assessments: Score, then rank." Research Integrity and Peer Review (2023).}
"ToyData3"
