% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrinkage_functions.R
\name{df_lasso_mult}
\alias{df_lasso_mult}
\title{Derivative of scaling and soft thresholding shrinkage function.}
\usage{
df_lasso_mult(x, params)
}
\arguments{
\item{x}{A vector of numerics.}

\item{params}{A vector of length 2 of numerics. \code{param[1]} is the
scaling parameter, \code{param[2]} is the thresholding parameter.}
}
\description{
Derivative of scaling and soft thresholding shrinkage function.
}
\seealso{
\code{\link{f_lasso_mult}}.
}
\author{
David Gerard.
}
