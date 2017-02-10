##' Update of thresholding parameter for scale and soft-thresholding HOSEs.
##'
##' @param c_obj The output from \code{\link{get_c}}.
##' @param lambda_current A vector of numerics of length \eqn{n}. The current
##'   values for the thresholding parameters.
##' @param c_current A positive numeric. The current value of the scaling
##'   parameter.
##' @param k A positive integer. The mode to update.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##' @param epsilon A positive numeric. The min distance of the Newton update.
##'
##' @return \code{lambda_new} A numeric. The update of the threshholding
##'   parameter.
##'
##' @seealso \code{\link{soft_coord}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
update_lambda <- function(c_obj, lambda_current, c_current, k, tau2,
                          epsilon = 10 ^ -4) {
  ## k = current mode to update
  hosvd_x <- c_obj$hosvd_x
  C_array <- c_obj$C_array
  S <- hosvd_x$S
  p <- dim(S)
  n <- length(p)
  sig <- hosvd_x$D

  max_toget <- rep(NA, length = n)
  for (mode_index in 1:n) {
    max_toget[mode_index] <-
      sum(sig[[mode_index]] > lambda_current[mode_index])
  }

  if (any(max_toget == 0)) {
    stop("Regularization is too large")
  }

  indices_toget <- as.matrix(expand.grid(lapply(max_toget, seq, from = 1)))
  C_sub <- array(C_array[indices_toget], dim = max_toget)
  S_sub <- array(S[indices_toget], dim = max_toget)

  f_d <- list()
  sig_sub <- list()
  one_over_f_sig <- list()
  for (mode_index in 1:n) {
    sig_sub[[mode_index]] <- sig[[mode_index]][1:max_toget[mode_index]]
    f_d[[mode_index]] <- sig_sub[[mode_index]] - lambda_current[mode_index]
    one_over_f_sig[[mode_index]] <-
      1 / (f_d[[mode_index]] * sig_sub[[mode_index]])
  }
  one_over_f_sig[[k]] <- rep(0, max_toget[k])



  sig_ik_list <- lapply(max_toget, rep, x = 1)
  sig_ik_list[[k]] <- sig_sub[[k]]

  f_d_ones_added <- f_d
  f_d_ones_added[[k]] <- rep(1, length = max_toget[k])
  f_sig_inv_array <- f_d_ones_added[[1]] / sig_sub[[1]]
  ratio_array <- one_over_f_sig[[1]]
  sig_ik_array <- sig_ik_list[[1]]
  for (mode_index in 2:n) {
    f_sig_inv_array <-
      outer(f_sig_inv_array,
            f_d_ones_added[[mode_index]] / sig_sub[[mode_index]])
    ratio_array <-
      outer(ratio_array, one_over_f_sig[[mode_index]], FUN = "+")
    sig_ik_array <- outer(sig_ik_array, sig_ik_list[[mode_index]])
  }

  b <- sum(c_current * f_sig_inv_array * S_sub ^ 2)
  e <- sum(tau2 * c_current * f_sig_inv_array * C_sub)
  d <- sum(tau2 * c_current * f_sig_inv_array * S_sub ^ 2 * ratio_array)
  a_1 <- sum(c_current ^ 2 * f_sig_inv_array ^ 2 * sig_ik_array * S_sub ^ 2)
  a_2 <- sum(c_current ^ 2 * f_sig_inv_array ^ 2 * S_sub ^ 2)

  lambda_new <- min( (a_1 + e + d - b) / a_2, max(sig[[k]]) - epsilon)
  if (lambda_new < 0) {
    lambda_new <- 0
  }
  return(lambda_new)
}



#' Apply Brent's method to update the kth mode when soft-thresholding.
#'
#' @inheritParams update_lambda
#'
#' @author David Gerard
#'
#' @export
update_lambda_brent <- function(c_obj, lambda_current, c_current, k, tau2) {
    ## k = current mode to update
    hosvd_x <- c_obj$hosvd_x
    C_array <- c_obj$C_array
    S <- hosvd_x$S
    p <- dim(S)
    n <- length(p)
    sig <- hosvd_x$D

    oout <- stats::optim(par = lambda_current[k], fn = sure_given_c_one_mode, method = "Brent",
                         lower = 0, upper = max(sig[[k]]), c_obj = c_obj, 
                         lambda_current = lambda_current, c_current = c_current,
                         current_mode = k, tau2 = tau2)      
    
     return(list(lambda_new = oout$par, sure = oout$value))
}

#' Wrapper for \code{\link{sure_given_c}} when only updating one mode.
#' 
#' @inheritParams update_lambda
#' @param lambda A numeric scalar. The current value of the lambda you are updating.
#' @param current_mode The current mode for which lambda is a member.
#' 
#' @author David Gerard
sure_given_c_one_mode <- function(lambda, c_obj, lambda_current, c_current, current_mode, tau2) {
  n <- length(lambda_current)
  
  ## set up the lasso function -----------------------------------------------
  func_lasso <- list()
  dfunc_lasso <- list()
  func_lasso[[1]] <- f_lasso_mult
  dfunc_lasso[[1]] <- df_lasso_mult
  for (mode_index in 2:n) {
    func_lasso[[mode_index]] <- f_lasso
    dfunc_lasso[[mode_index]] <- df_lasso
  }
  
  lambda_current[current_mode] <- lambda
  
  lambda_list <- list()
  lambda_list[[1]] <- c(lambda_current[1], c_current)
  for (index in 2:n) {
    lambda_list[[index]] <- lambda_current[index]
  }

  sout <- sure_given_c(obj = c_obj, func = func_lasso, dfunc = dfunc_lasso,
                       lambda = lambda_list, tau2 = tau2)
  return(sout$sure_val)
}



##' Update scale parameter in scale and soft-thresholding HOSE's.
##'
##' @param c_obj The output from \code{\link{get_c}}.
##' @param lambda_current A vector of numerics of length \eqn{n}. The initial
##'   starting values for the thresholding parameters.
##' @param c_current A positive numeric. The starting value of the scaling
##'   parameter.
##' @param epsilon A positive numeric. The min distance of the Newton update.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @return \code{c_new} A postive numeric. The update of the scaling parameter.
##'
##' @seealso \code{\link{soft_coord}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
update_c <- function(c_obj, lambda_current, c_current, tau2,
                     epsilon = 10 ^ -4) {
  hosvd_x <- c_obj$hosvd_x
  C_array <- c_obj$C_array
  S <- hosvd_x$S
  p <- dim(S)
  n <- length(p)

  sig <- hosvd_x$D

  max_toget <- rep(NA, length = n)
  for (mode_index in 1:n) {
    max_toget[mode_index] <-
      sum(sig[[mode_index]] > lambda_current[mode_index])
  }

  indices_toget <- as.matrix(expand.grid(lapply(max_toget, seq, from = 1)))
  C_sub <- array(C_array[indices_toget], dim = max_toget)
  S_sub <- array(S[indices_toget], dim = max_toget)

  f_d <- list()
  sig_sub <- list()
  for (mode_index in 1:n) {
    sig_sub[[mode_index]] <- sig[[mode_index]][1:max_toget[mode_index]]
    f_d[[mode_index]] <- sig_sub[[mode_index]] - lambda_current[mode_index]
  }

  f_sig_inv_array <- f_d[[1]] / sig_sub[[1]]
  ratio_array <- 1 / (f_d[[1]] * sig_sub[[1]])
  for (mode_index in 2:n) {
    f_sig_inv_array <-
      outer(f_sig_inv_array, f_d[[mode_index]] / sig_sub[[mode_index]])
    ratio_array <-
      outer(ratio_array,
            1 / (f_d[[mode_index]] * sig_sub[[mode_index]]), FUN = "+")
  }

  a <- sum(f_sig_inv_array ^ 2 * S_sub ^ 2)
  b <- sum(f_sig_inv_array * S_sub ^ 2)
  d <- sum(tau2 * f_sig_inv_array * C_sub)
  e <- sum(tau2 * f_sig_inv_array * S_sub ^ 2 * ratio_array)

  c_new <- max( (b - d - e) / a, epsilon)
  return(c_new)
}

##' Runs an iterative coordinate descent algorithm to minimize the SURE for
##' mode-specific soft thresholding estimators.
##'
##' @param c_obj The output from \code{\link{get_c}}.
##' @param lambda_init A vector of numerics of length \eqn{n}. The initial
##'   starting values for the thresholding parameters.
##' @param c_init A positive numeric. The starting value of the scaling
##'   parameter.
##' @param itermax A positive integer. The maximum number of Newton steps to
##'   iterate through.
##' @param tol A positive numeric. The stopping criterion.
##' @param print_iter A logical. Should we print the updates of the Newton Step?
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##' @param use_sure A logical. Which stopping criterion should we use? The mean
##'   absolute difference in the parameters (\code{FALSE}) or the absolute value
##'   of the deviation of the ratio of adjacent SURE values from 1
##'   (\code{TRUE}).
##'
##' @return \code{c} A numeric. The final value of the scaling parameter.
##'
##'   \code{lambda} A vector of numerics. The final values of the thresholding
##'   parameters.
##'
##'   \code{est} An array of numerics. The final mean estimate.
##'
##' @author David Gerard.
##'
##' @seealso \code{\link{update_c}}, \code{\link{update_lambda}},
##'   \code{\link{get_c}}.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
soft_coord <- function(c_obj, lambda_init, c_init, itermax = 1000,
                       tol = 10 ^ -4, print_iter = TRUE, tau2 = 1,
                       use_sure = TRUE) {
  hosvd_x <- c_obj$hosvd_x
  ##C_array <- c_obj$C_array ## may not be used.
  S <- hosvd_x$S
  p <- dim(S)
  n <- length(p)


  func_lasso <- list()
  dfunc_lasso <- list()
  func_lasso[[1]] <- f_lasso_mult
  dfunc_lasso[[1]] <- df_lasso_mult
  for (mode_index in 2:n) {
    func_lasso[[mode_index]] <- f_lasso
    dfunc_lasso[[mode_index]] <- df_lasso
  }


  lambda_current <- lambda_init
  c_current <- c_init
  current_sure <- prod(p)
  iter_index <- 1
  error_all <- tol + 1
  while (error_all > tol & iter_index < itermax) {
    c_old <- c_current
    lambda_old <- lambda_current
    old_sure <- current_sure
    c_current <- update_c(c_obj, lambda_current, c_current, tau2)
    for (k in 1:n) {
      bout <- update_lambda_brent(c_obj, lambda_current, c_current, k, tau2)
      lambda_current[k] <- bout$lambda_new
      current_sure <- bout$sure
    }

    if (use_sure) {
      cat("     SURE =", round(current_sure, digits = 2), "\n")
      error_all <- abs(1 - current_sure / old_sure)
    } else {
      error_all <- sum(abs(lambda_current - lambda_old)) +
        abs(1 - c_current / c_old)
    }

    if (print_iter) {
      cat("Iteration =", iter_index, "\n")
      cat("        c =", round(c_current, digits = 2), "\n")
      cat("   Lambda =",
          paste("(", paste(round(lambda_current, digits = 2),
                           collapse = ","), ")", sep = ""), "\n\n")
    }
    iter_index <- iter_index + 1
  }

  ## Get final mean estimate.
  param_list <- list()
  param_list[[1]] <- c(lambda_current[1], c_current)
  for (mode_index in 2:n) {
    param_list[[mode_index]] <- lambda_current[mode_index]
  }
  est <- sure_given_c(obj = c_obj, func = func_lasso, dfunc = dfunc_lasso,
                      lambda = param_list, tau2 = tau2)$mean_est

  return(list(c = c_current, lambda = lambda_current, est = est))
}
