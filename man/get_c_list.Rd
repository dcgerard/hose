% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block_sure.R
\name{get_c_list}
\alias{get_c_list}
\title{Calculates the necessary components to calculate the SURE for estimators
that block shrink overlapping blocks.}
\usage{
get_c_list(current_tensor, block_size = dim(current_tensor)[3])
}
\arguments{
\item{current_tensor}{An array of numerics.}

\item{block_size}{The size of the overlapping blocks.}
}
\value{
A list of elements used in \code{\link{block_sure_given_c_list}} to
  calculate the SURE.
}
\description{
Calculates the necessary components to calculate the SURE for estimators
that block shrink overlapping blocks.
}
\references{
Gerard, D., & Hoff, P. (2015).
  \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
  Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
}
\seealso{
\code{\link{block_sure_given_c_list}}, \code{\link{block_sure}}.
}
\author{
David Gerard.
}
