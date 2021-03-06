% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block_sure.R
\name{block_sure_given_c_list}
\alias{block_sure_given_c_list}
\title{Calculates the sure of estimators that shrink subtensors using the output of
\code{\link{get_c_list}}.}
\usage{
block_sure_given_c_list(c_obj, func, dfunc, lambda_list, tau2 = 1)
}
\arguments{
\item{c_obj}{The output from \code{\link{get_c_list}}.}

\item{func}{A list of functions to apply to each mode.}

\item{dfunc}{A list of derivatives of the function in \code{func}.}

\item{lambda_list}{A list of parameter values for \code{func} and
\code{dfunc}.}

\item{tau2}{A positive numeric. The variance. Assumed known and defaults to
1.}
}
\value{
\code{sure_final} A numeric. The SURE value.

\code{mean_est} An array of numerics. The mean estimate.
}
\description{
Calculates the sure of estimators that shrink subtensors using the output of
\code{\link{get_c_list}}.
}
\references{
Gerard, D., & Hoff, P. (2015).
  \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
  Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
}
\seealso{
\code{\link{get_c_list}}, \code{\link{block_sure}}.
}
\author{
David Gerard.
}
