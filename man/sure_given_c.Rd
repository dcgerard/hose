% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mode_specific_sure.R
\name{sure_given_c}
\alias{sure_given_c}
\title{Calculates SURE given the output of \code{\link{get_c}}.}
\usage{
sure_given_c(obj, func, dfunc, lambda, tau2 = 1)
}
\arguments{
\item{obj}{Output from \code{\link{get_c}}.}

\item{func}{A list of length \code{length(dim(X))} of shrinkage functions.}

\item{dfunc}{A list of length \code{length(dim(X))} of corresponding
derivative functions.}

\item{lambda}{A list of parameter values for shrinking along each mode.}

\item{tau2}{A positive numeric. The variance. Assumed known and defaults to
1.}
}
\description{
Calculates SURE assuming we already did all of the heavy lifting in
calculating the HOSVD and the 'C' matrix from my write-up 'sure_pdf'. Call
\code{\link{get_c}} before using this function.
}
\references{
Gerard, D., & Hoff, P. (2015).
\href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral Estimators}.
\emph{arXiv preprint} arXiv:1505.02114.
}
\seealso{
\code{\link{get_c}}, \code{\link{diverge_given_c}}, \code{\link{sure}}.
}
\author{
David Gerard.
}
