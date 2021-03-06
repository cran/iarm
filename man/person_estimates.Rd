% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Person-Fit.R
\name{person_estimates}
\alias{person_estimates}
\title{Person Estimates with MLE and WLE}
\usage{
person_estimates(object, properties = F, allperson = F)
}
\arguments{
\item{object}{An object of class "Rm", a fitted Rasch model or partial
credit model using  the functions RM or PCM in package eRm, or an object of class "raschmodel" or "pcmodel",
 a fitted Rasch model or partial credit model using the functions raschmodel or pcmodel in package psychotools.}

\item{properties}{If TRUE additional properties of the estimates are given (see below).}

\item{allperson}{If TRUE person estimates (MLE and WLE) for all persons in the data set are delivered.}
}
\value{
If properties = False a matrix containing:
\item{Raw score}{raw score}
\item{MLE}{MLE of person parameters}
\item{WLE}{WLE of person parameters}

If properties = TRUE a list with two components, one for MLE and the other for WLE. Each component
contains:
\item{Raw score}{raw score}
\item{MLE or WLE}{person estimates}
\item{SEM}{standard error of measurement}
\item{Bias}{bias}
\item{RMSE}{root mean square error}
\item{Score.SEM}{score sem}
}
\description{
Computes Person estimates with maximum likelihood estimation (MLE) and  weighted likelihood estimation (WLE) for raw scores 0 to m.
}
\examples{
rm.mod <- RM(amts[,4:13])
person_estimates(rm.mod)
}
\references{
Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
\emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 63 - 70.
}
\author{
Marianne Mueller
}
