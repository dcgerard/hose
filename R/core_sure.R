##' Calculates necessary components to calculate the SURE for estimators that
##' shrink the individual elements of the core array.
##'
##' @param X An array of numerics. The tensor data.
##'
##' @return A list of objects meant to be fed into
##'   \code{\link{sure_given_basics}} to calculate the SURE for higher-order
##'   spectral estimators that shrink the core array directly.
##'
##' @author David Gerard.
##'
##' @seealso \code{\link{sure_given_basics}}, \code{\link{sure_S}}
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
get_basics <- function(X) {
  p <- dim(X)
  hosvd_x <- hosvd_full(X)

  sum_array <- list()
  for (mode_index in 1:length(p)) {
    s_k_mat <- matrix(NA, nrow = p[mode_index], ncol = prod(p[-mode_index]))
    for (i_index in 1:p[mode_index]) {
      if (p[mode_index] > 2) {
        s_k_mat[i_index, ] <-
          colSums((1 / (hosvd_x$D[[mode_index]][i_index] ^ 2 -
                          hosvd_x$D[[mode_index]] ^ 2) *
                     tensr::mat(hosvd_x$S ^ 2, mode_index))[-i_index, ])
      } else {
        s_k_mat[i_index, ] <-
          (1 / (hosvd_x$D[[mode_index]][i_index] ^ 2 -
                  hosvd_x$D[[mode_index]] ^ 2) *
             tensr::mat(hosvd_x$S ^ 2, mode_index))[-i_index, ]
      }
    }
    s_k_temp_array <-
      array(s_k_mat, dim = c(p[mode_index], p[-mode_index]))
    sum_array[[mode_index]] <-
      aperm(s_k_temp_array,
            match(1:length(p), c(mode_index, (1:length(p))[-mode_index])))
  }
  S_mult <- listSum(sum_array)


  sig_list <- list()
  for (mode_index in 1:length(p)) {
    sig_list[[mode_index]] <- rep(NA, length = p[mode_index])
    for (i_index in 1:p[mode_index]) {
      sig_list[[mode_index]][i_index] <-
        sum( (1 / (hosvd_x$D[[mode_index]][i_index] ^ 2 -
                     hosvd_x$D[[mode_index]] ^ 2))[-i_index])
    }
    if (mode_index == 1) {
      sig_leave_one_out <- sig_list[[mode_index]]
    } else {
      sig_leave_one_out <-
        outer(sig_leave_one_out, sig_list[[mode_index]], FUN = "+")
    }
  }
  s_s_array <- sig_leave_one_out * hosvd_x$S

  return(list(hosvd_x = hosvd_x, s_s_array = s_s_array, S_mult = S_mult))
}

##' Calculates the SURE for estimators that shrink the individual elements of
##' the core array.
##'
##' @param c_obj The output from \code{\link{get_basics}}.
##' @param func A shrinkage function.
##' @param dfunc The derivative function of \code{func}.
##' @param params The parameters to be fed into \code{func} and \code{dfunc}.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @return \code{sure_val} The sure value of the higher-order spectral
##'   estimator.
##'
##'   \code{est} The mean estimate.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{get_basics}}, \code{\link{sure_S}}.
##'
##' @export
sure_given_basics <- function(c_obj, func, dfunc, params, tau2 = 1) {
  ## c_obj: object returned from get_basics() tau2: known variance
  hosvd_x <- c_obj$hosvd_x
  s_s_array <- c_obj$s_s_array
  S_mult <- c_obj$S_mult

  p <- dim(s_s_array) ## need to check this.

  s_new <- func(hosvd_x$S, params)
  diff_S <- dfunc(hosvd_x$S, params)
  est <- tensr::atrans(s_new, hosvd_x$U)
  divergence_f <- sum(s_s_array * s_new + diff_S * (1 + S_mult))

  sure_val <- -tau2 * prod(p) + 2 * tau2 * divergence_f +
    sum( (s_new - hosvd_x$S) ^ 2)
  return(list(sure_val = sure_val, est = est))
}


##' Wrapper for \code{\link{get_basics}} and \code{\link{sure_given_basics}}.
##'
##' @inheritParams get_basics
##' @inheritParams sure_given_basics
##'
##' @return \code{sure_val} The sure value of the higher-order spectral
##'   estimator.
##'
##'   \code{est} The mean estimate.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{get_basics}}, \code{\link{sure_given_basics}}.
##'
##' @export
sure_S <- function(X, func, dfunc, params, tau2 = 1) {
  c_obj <- get_basics(X)
  sure_given_basics(c_obj, func, dfunc, params, tau2)
}
