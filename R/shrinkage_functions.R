##' Truncation spectral function.
##'
##' @param x A vector of numerics.
##' @param r A positive integer.
##'
##' @seealso \code{\link{df_truncate}}.
##'
##' @author David Gerard.
##'
##' @export
f_truncate <- function(x, r) {
  n <- length(x)
  y <- rep(0, length = n)
  y[1:r] <- x[1:r]
  return(y)
}

##' Derivative of truncation spectral function.
##'
##' @inheritParams f_truncate
##'
##' @seealso \code{\link{f_truncate}}.
##'
##' @author David Gerard.
##'
##' @export
df_truncate <- function(x, r) {
  n <- length(x)
  y <- rep(0, length = n)
  y[1:r] <- 1
  return(y)
}

##' Soft thresholding shrinkage function. Same as lasso for spectral shrinkage.
##'
##' @param x A vector of numerics.
##' @param lambda A numeric. The threshholding parameter.
##'
##' @seealso \code{\link{df_lasso}}.
##'
##' @author David Gerard.
##'
##' @export
f_lasso <- function(x, lambda) {
  return(pos_part(x - lambda))
}

##' Derivative of soft thresholding shrinkage function.
##'
##' @inheritParams f_lasso
##'
##' @seealso \code{\link{f_lasso}}.
##'
##' @author David Gerard.
##'
##' @export
df_lasso <- function(x, lambda) {
  y <- rep(1, length = length(x))
  y[x < lambda] <- 0
  return(y)
}

##' Scaling and soft thresholding shrinkage function.
##'
##' @param x A vector of numerics.
##' @param params A vector of length 2 of numerics. \code{param[1]} is the
##'   scaling parameter, \code{param[2]} is the threshholding parameter.
##'
##' @seealso \code{\link{df_lasso_mult}}.
##'
##' @author David Gerard.
##'
##' @export
f_lasso_mult <- function(x, params) {
  lambda <- params[1]
  const <- params[2]
  return(const * pos_part(x - lambda))
}

##' Derivative of scaling and soft thresholding shrinkage function.
##'
##'
##' @inheritParams f_lasso_mult
##'
##' @seealso \code{\link{f_lasso_mult}}.
##'
##' @author David Gerard.
##'
##' @export
df_lasso_mult <- function(x, params) {
  lambda <- params[1]
  const <- params[2]
  y <- rep(const, length = length(x))
  y[x < lambda] <- 0
  return(y)
}

##' Soft thresholding a core array.
##'
##' @param S An array of numerics.
##' @param lambda A numeric. The threshholding parameter.
##'
##' @seealso \code{\link{df_S_lasso}}.
##'
##' @author David Gerard.
##'
##' @export
f_S_lasso <- function(S, lambda) {
  ## lasso typo estimator
  S_new <- sign(S) * pos_part_2(abs(S) - lambda)
  return(S_new)
}

##' Derivative of soft thresholding a core array.
##'
##' @inheritParams f_S_lasso
##'
##' @seealso \code{\link{df_S_lasso}}.
##'
##' @author David Gerard.
##'
##' @export
df_S_lasso <- function(S, lambda) {
  ## derivative of lasso type estimator
  diff_S <- array(1, dim = dim(S))
  diff_S[abs(S) < lambda] <- 0
  return(diff_S)
}

##' Scaling and soft thresholding a core array.
##'
##' @param S An array of numerics.
##' @param params A vector of length 2 of numerics. \code{param[1]} is the
##'   scaling parameter, \code{param[2]} is the threshholding parameter.
##'
##' @seealso \code{\link{df_S_lasso_mult}}.
##'
##' @author David Gerard.
##'
##' @export
f_S_lasso_mult <- function(S, params) {
  ## lasso typo estimator
  lambda <- params[1]
  const <- params[2]
  S_new <- const * sign(S) * pos_part_2(abs(S) - lambda)
  return(S_new)
}

##' Derivative of scaling and soft thresholding a core array.
##'
##' @inheritParams f_S_lasso_mult
##'
##' @seealso \code{\link{df_S_lasso_mult}}.
##'
##' @author David Gerard.
##'
##' @export
df_S_lasso_mult <- function(S, params) {
  ## lasso typo estimator
  lambda <- params[1]
  const <- params[2]
  diff_S <- const * array(1, dim = dim(S))
  diff_S[abs(S) < lambda] <- 0
  return(diff_S)
}
