#' Real peer review data set from the American Institute of Biological Sciences (AIBS)
#'
#' This real data set includes 12 judges (reviewers) and 28 objects (proposals), and demonstrates the ability
#' of the Mallows-Binomial model to combine ratings and rankings for the purpose of demarcating
#' real grant proposals for a funding agency.
#'
#' @format A list with three elements: (1) \code{rankings}, a 12 x 18 matrix of rankings with one row per judge;
#'   (2) \code{ratings}, a 12 x 18 matrix of ratings, with one row per judge and one column per object; and
#'   (3) \code{M}, a number indicating the maximum (worst) integer score.
#' @source {Originally published in: Gallo, Stephen A.. "Grant Peer Review Scoring Data with Criteria Scores" (2023). \url{https://figshare.com/articles/dataset/Grant_Peer_Review_Scoring_Data_with_Criteria_Scores/12728087/1}.}
#' @source {Originally analyzed in: Gallo, Stephen A., et al. "A new approach to peer review assessments: Score, then rank." Research Integrity and Peer Review (2023).}
"AIBS"
