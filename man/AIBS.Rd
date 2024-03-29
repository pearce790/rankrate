% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AIBS.R
\docType{data}
\name{AIBS}
\alias{AIBS}
\title{Real peer review data set from the American Institute of Biological Sciences (AIBS)}
\format{
A list with three elements: (1) \code{rankings}, a 12 x 18 matrix of rankings with one row per judge;
(2) \code{ratings}, a 12 x 18 matrix of ratings, with one row per judge and one column per object; and
(3) \code{M}, a number indicating the maximum (worst) integer score.
}
\source{
{Originally published in: Gallo, Stephen A.. "Grant Peer Review Scoring Data with Criteria Scores" (2023). \url{https://figshare.com/articles/dataset/Grant_Peer_Review_Scoring_Data_with_Criteria_Scores/12728087/1}.}

{Originally analyzed in: Gallo, Stephen A., et al. "A new approach to peer review assessments: Score, then rank" (2023). Research Integrity and Peer Review 8:10 (10). \url{https://researchintegrityjournal.biomedcentral.com/articles/10.1186/s41073-023-00131-7}.}
}
\usage{
AIBS
}
\description{
This real data set includes 12 judges (reviewers) and 28 objects (proposals), and demonstrates the ability
of the Mallows-Binomial model to combine ratings and rankings for the purpose of demarcating
real grant proposals for a funding agency.
}
\keyword{datasets}
