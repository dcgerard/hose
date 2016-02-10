##' SURE for spectral shrinkage functions for complex matrices.
##'
##' @param d A vector of numerics. The singular values.
##' @param lambda A numeric. The shrinkage parameter.
##' @param p A vector of positive integers. The dimension of the data array.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @author David Gerard
##'
##' @references Candes, E., Sing-Long, C. A., & Trzasko, J. D. (2013).
##'   \href{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6545395&tag=1}{Unbiased
##'    risk estimates for singular value thresholding and spectral estimators}.
##'   \emph{Signal Processing, IEEE Transactions on}, 61(19), 4643-4657.
##'
##' @export
sure_matrix_complex <- function(d, lambda, p, tau2 = 1) {

  f_d <- pos_part(d - lambda)
  df_d <- 1 * (d > lambda)

  max_d_lambda <- sapply(d, min, lambda)
  outer_temp <- 1 / outer(d ^ 2, d ^ 2, "-")
  diag(outer_temp) <- 0
  right_mult <- apply(outer_temp, 1, sum)
  div_f <- sum(df_d + (2 * abs(p[1] - p[2]) + 1) * f_d / d) + 4 *
    sum(d * f_d * right_mult)
  -2 * tau2 * prod(p) + sum(max_d_lambda ^ 2) + 2 * tau2 * div_f
}

##' SURE for spectral shrinkage functions for real matrices.
##'
##' Only works for soft thresholding.
##'
##' @param d A vector of numerics. The singular values.
##' @param lambda A numeric. The shrinkage parameter.
##' @param p_dim A vector of positive integers. The dimension of the data array.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @author David Gerard.
##'
##' @references Candes, E., Sing-Long, C. A., & Trzasko, J. D. (2013).
##'   \href{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6545395&tag=1}{Unbiased
##'    risk estimates for singular value thresholding and spectral estimators}.
##'   \emph{Signal Processing, IEEE Transactions on}, 61(19), 4643-4657.
##'
##' @export
sure_matrix <- function(d, lambda, p_dim, tau2 = 1) {

  f_d <- pos_part(d - lambda)
  df_d <- 1 * (d > lambda)

  max_d_lambda <- sapply(d, min, lambda)
  outer_temp <- 1 / outer(d ^ 2, d ^ 2, "-")
  diag(outer_temp) <- 0
  right_mult <- apply(outer_temp, 1, sum)
  ##div_f <- sum(df_d + abs(p_dim[1] - p_dim[2]) *
  ##              pos_part(1 - lambda / d)) + 2 *
  ## sum(d * f_d * right_mult)
  div_f <- sum(df_d + abs(p_dim[1] - p_dim[2]) * f_d / d) +
    2 * sum(d * f_d * right_mult)
  -tau2 * prod(p_dim) + sum(max_d_lambda ^ 2) + 2 * tau2 * div_f
}

##' Min Efron-Morris estimator.
##'
##' @param X An array of numerics. The data tensor.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @author David Gerard
##'
##' @references Efron, B., & Morris, C. (1972).
##'   \href{http://biomet.oxfordjournals.org/content/59/2/335.short}{Empirical
##'   Bayes on vector observations: An extension of Stein's method}.
##'   \emph{Biometrika}, 59(2), 335-347.
##'
##' @export
min_ef <- function(X, tau2 = 1) {
  ## minimizes SURE of efron-moris estimator along each mode
  p <- dim(X)
  n <- length(p)
  ef_est <- list()
  sure_vec <- rep(NA, length = n)
  for (mode_index in 1:n) {
    S_k_inv <- solve(tensr::mat(X, mode_index) %*%
                       t(tensr::mat(X, mode_index)))
    sure_vec[mode_index] <- tau2 * 1 -
      (prod(p[-mode_index]) - p[mode_index] - 1) ^ 2 *
      tau2 ^ 2 / prod(p) * sum(diag(S_k_inv))
    ef_est[[mode_index]] <-
      tensr::amprod(X, (diag(p[mode_index]) -
                          tau2 * (prod(p[-mode_index]) - p[mode_index] - 1) *
                          S_k_inv), mode_index)
  }
  winner <- which.min(sure_vec)
  return(list(ef_est = ef_est[[winner]], sure = prod(p) * sure_vec[winner]))
}

##' Efron-Morris estimator.
##'
##' @inheritParams min_ef
##' @param mode_index The mode upon which to apply the Efron-Morris estimator.
##'   Defaults to the first mode.
##'
##' @author David Gerard.
##'
##' @references Efron, B., & Morris, C. (1972).
##'   \href{http://biomet.oxfordjournals.org/content/59/2/335.short}{Empirical
##'   Bayes on vector observations: An extension of Stein's method}.
##'   \emph{Biometrika}, 59(2), 335-347.
##'
##' @export
efron_morris <- function(X, mode_index = 1, tau2 = 1) {
  p <- dim(X)
  S_k_inv <- solve(tensr::mat(X, mode_index) %*%
                     t(tensr::mat(X, mode_index)))
  ef_est <- tensr::amprod(X, (diag(p[mode_index]) - tau2 *
                                (prod(p[-mode_index]) - p[mode_index] - 1) *
                                S_k_inv), mode_index)
  sure_val <- tau2 * 1 - (prod(p[-mode_index]) - p[mode_index] - 1) ^ 2 *
    tau2 ^ 2 / prod(p) * sum(diag(S_k_inv))
  return(list(ef_est = ef_est, sure = prod(p) * sure_val))
}

##' Stein's estimator.
##'
##' @param X A vector, matrix, or array of numerics. The data.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @author David Gerard.
##'
##' @references James, W., & Stein, C. (1961, June).
##'   \href{http://projecteuclid.org/euclid.bsmsp/1200512173}{Estimation with
##'   quadratic loss}. \emph{In Proceedings of the fourth Berkeley symposium on
##'   mathematical statistics and probability} (Vol. 1, No. 1961, pp. 361-379).
##'
##' @export
stein <- function(X, tau2 = 1) {
  p <- dim(X)
  ## returns Stein's estimator
  sum_x <- sum(X ^ 2)
  est <- (1 - (prod(p) - 2) * tau2 / sum_x) * X
  sure_est <- prod(p) * tau2 - (prod(p) - 2) ^ 2 * tau2 ^ 2 / sum_x
  return(list(est = est, sure_est = sure_est))
}
