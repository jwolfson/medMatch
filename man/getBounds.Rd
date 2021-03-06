% Generated by roxygen2 (4.0.0): do not edit by hand
\name{getBounds}
\alias{getBounds}
\title{Calculates matching-based bounds on the natural direct effect.}
\usage{
getBounds(DAT, inds, Gamma = 1, Psi = Gamma, effect = "ATT")
}
\arguments{
\item{DAT}{A data frame containing columns for the treatment (\code{Z}), outcome (\code{Y}),
potential mediator (\code{S} and \code{S.cat}), and baseline covariates (\code{W1, W2, ...})}

\item{inds}{A vector of indicators for the bootstrap.}

\item{Gamma}{The sensitivity parameter controlling the odds ratio between the true and estimated odds of a discordant pair favoring the treatment.}

\item{Psi}{The sensitivity parameter controlling the odds ratio between the true and estimated odds of a pair being discordant.}

\item{effect}{The causal effect to be estimated. Can be \code{ATT} or \code{ATC}.}
}
\value{
A vector of three elements containing: \enumerate{
\item The estimated lower bound
\item The estimated value
\item The estimated upper bound
}
}
\description{
Internal function for calculating matching-based bounds for the NDE.
}
\details{
This function is called by the bootstrap procedure in the main \code{\link{DeltaBounds}} function.
}

