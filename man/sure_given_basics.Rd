% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_sure.R
\name{sure_given_basics}
\alias{sure_given_basics}
\title{Calculates the SURE for estimators that shrink the individual elements of
the core array.}
\usage{
sure_given_basics(c_obj, func, dfunc, params, tau2 = 1)
}
\arguments{
\item{c_obj}{The output from \code{\link{get_basics}}.}

\item{func}{A shrinkage function.}

\item{dfunc}{The derivative function of \code{func}.}

\item{params}{The parameters to be fed into \code{func} and \code{dfunc}.}

\item{tau2}{A positive numeric. The variance. Assumed known and defaults to
1.}
}
\value{
\code{sure_val} The sure value of the higher-order spectral
  estimator.

  \code{est} The mean estimate.
}
\description{
Calculates the SURE for estimators that shrink the individual elements of
the core array.
}
\references{
Gerard, D., & Hoff, P. (2015).
  \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
  Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
}
\seealso{
\code{\link{get_basics}}, \code{\link{sure_S}}.
}
\author{
David Gerard.
}
