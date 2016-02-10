##' Outputs the 'C' array from Gerard and Hoff (2015), along with the HOSVD of
##' the data tensor.
##'
##' This is necessary to calculate the SURE for higher-order spectral
##' estimators.
##'
##' this function only calculates C and the HOSVD so that calculating SURE is
##' faster when doing it over and over again
##'
##' @param X An array of numerics. The data.
##'
##' @return \code{C_array} An array of numerics. The "C" array from Gerard and
##'   Hoff (2015).
##'
##'   \code{hosvd_x} A list containing the higher-order SVD of \code{X}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{diverge_given_c}} for calculating the divergence of
##'   higher-order spectral estimators using the output of \code{get_c}.
##'
##'   \code{\link{sure_given_c}} for calculating the SURE given the output of
##'   \code{get_c}.
##'
##'   \code{\link{sure}} for a wrapper for \code{get_c} and
##'   \code{\link{sure_given_c}}.
##'
##'   \code{\link{soft_coord}} for a coordinate descent algorithm for finding
##'   the optimal sure using the output from \code{get_c}.
##'
##' @export
get_c <- function(X) {

  p <- dim(X)
  n <- length(p)

  hosvd_x <- hosvd_full(X)

  S <- hosvd_x$S

  for (mode_index in 1:n) {
    over_temp <- rep(NA, length = p[mode_index])
    for (individual_index in 1:p[mode_index]) {
      over_temp[individual_index] <-
        sum(1 / (hosvd_x$D[[mode_index]] ^ 2 -
                   hosvd_x$D[[mode_index]][individual_index] ^
                   2)[-individual_index])
    }
    if (mode_index == 1) {
      sig_inv_sq_array <- 1 / hosvd_x$D[[1]] ^ 2
      sig_leave_one_out <- over_temp
    } else {
      sig_inv_sq_array <-
        outer(sig_inv_sq_array,
              1 / hosvd_x$D[[mode_index]] ^ 2, FUN = "+")
      sig_leave_one_out <- outer(sig_leave_one_out, over_temp, FUN = "+")
    }
  }

  D_list <- list()
  for (k_index in 1:n) {
    k_temp_mat <- matrix(NA, nrow = p[k_index], ncol = prod(p[-k_index]))
    for (i_index in 1:p[k_index]) {
      mult_left <- 1 / (hosvd_x$D[[k_index]][i_index] ^ 2 -
                          hosvd_x$D[[k_index]] ^ 2)
      if (p[k_index] > 2) {
        k_temp_mat[i_index, ] <-
          colSums((mult_left * tensr::mat(S ^ 2, k_index))[-i_index, ])
      } else {
        k_temp_mat[i_index, ] <-
          (mult_left * tensr::mat(S ^ 2, k_index))[-i_index, ]
      }
    }
    k_temp_array <- array(k_temp_mat, dim = c(p[k_index], p[-k_index]))
    D_list[[k_index]] <-
      aperm(k_temp_array, match(1:n, c(k_index, (1:n)[-k_index])))
  }

  ## an array of ones
  one_array <- array(1, dim = p)

  C_array <- one_array - S ^ 2 * (sig_inv_sq_array + sig_leave_one_out) +
    listSum(D_list)

  return(list(C_array = C_array, hosvd_x = hosvd_x))
}

##' Calculates divergence of HOSE.
##'
##' Assumes we already did all of the heavy lifting in calculating the HOSVD and
##' the 'C' matrix from my write-up 'sure_pdf'. Call \code{\link{get_c}} before
##' using this function.
##'
##' @inheritParams sure_given_c
##'
##' @return \code{divergence_f} The divergence of the higher-order spectral
##'   estimator.
##'
##'   \code{mean_est} The mean estimate of the higher-order spectral estimator.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{get_c}}, \code{\link{sure_given_c}}.
##'
##' @export
diverge_given_c <- function(obj, func, dfunc, lambda) {
  hosvd_x <- obj$hosvd_x
  C_array <- obj$C_array

  S <- hosvd_x$S

  p <- dim(S)
  n <- length(p)

  D <- list()
  D_inv <- list()
  f_D <- list()
  df_D <- list()
  for (mode_index in 1:n) {
    D[[mode_index]] <- diag(hosvd_x$D[[mode_index]])
    D_inv[[mode_index]] <- diag(1 / hosvd_x$D[[mode_index]])
    f_D[[mode_index]] <-
      diag(do.call(func[[mode_index]],
                   args = list(hosvd_x$D[[mode_index]],
                               lambda[[mode_index]])))
    df_D[[mode_index]] <-
      diag(do.call(dfunc[[mode_index]],
                   list(hosvd_x$D[[mode_index]], lambda[[mode_index]])))
  }

  f_d_inv <- list()
  for (mode_index in 1:n) {
    f_d_inv[[mode_index]] <- f_D[[mode_index]] * D_inv[[mode_index]]
  }


  H_list <- list()
  for (mode_index in 1:n) {
    H_list[[mode_index]] <- f_d_inv
    H_list[[mode_index]][[mode_index]] <-
      df_D[[mode_index]] * D_inv[[mode_index]] ^ 2
  }

  divergence_f <- sum(tensr::atrans(C_array, f_d_inv))
  for (mode_index in 1:n) {
    divergence_f <- divergence_f +
      sum(tensr::atrans(S ^ 2, H_list[[mode_index]]))
  }

  mean_est <- tensr::atrans(tensr::atrans(S, f_d_inv), hosvd_x$U)

  return(list(divergence_f = divergence_f, mean_est = mean_est))
}

##' Calculates SURE given the output of \code{\link{get_c}}.
##'
##' Calculates SURE assuming we already did all of the heavy lifting in
##' calculating the HOSVD and the 'C' matrix from my write-up 'sure_pdf'. Call
##' \code{\link{get_c}} before using this function.
##'
##' @param obj Output from \code{\link{get_c}}.
##' @param func A list of length \code{length(dim(X))} of shrinkage functions.
##' @param dfunc A list of length \code{length(dim(X))} of corresponding
##'   derivative functions.
##' @param lambda A list of parameter values for shinking along each mode.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##' \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral Estimators}.
##' \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{get_c}}, \code{\link{diverge_given_c}}, \code{\link{sure}}.
##'
##' @export
sure_given_c <- function(obj, func, dfunc, lambda, tau2 = 1) {
  hosvd_x <- obj$hosvd_x
  C_array <- obj$C_array

  S <- hosvd_x$S

  p <- dim(S)
  n <- length(p)

  D <- list()
  D_inv <- list()
  f_D <- list()
  df_D <- list()
  for (mode_index in 1:n) {
    D[[mode_index]] <- diag(hosvd_x$D[[mode_index]])
    D_inv[[mode_index]] <- diag(1 / hosvd_x$D[[mode_index]])
    f_D[[mode_index]] <-
      diag(do.call(func[[mode_index]],
                   args = list(hosvd_x$D[[mode_index]],
                               lambda[[mode_index]])))
    df_D[[mode_index]] <-
      diag(do.call(dfunc[[mode_index]], list(hosvd_x$D[[mode_index]],
                                             lambda[[mode_index]])))
  }

  f_d_inv <- list()
  for (mode_index in 1:n) {
    f_d_inv[[mode_index]] <- f_D[[mode_index]] * D_inv[[mode_index]]
  }


  H_list <- list()
  for (mode_index in 1:n) {
    H_list[[mode_index]] <- f_d_inv
    H_list[[mode_index]][[mode_index]] <-
      df_D[[mode_index]] * D_inv[[mode_index]] ^ 2
  }

  divergence_f <- sum(tensr::atrans(C_array, f_d_inv))
  for (mode_index in 1:n) {
    divergence_f <- divergence_f +
      sum(tensr::atrans(S ^ 2, H_list[[mode_index]]))
  }

  rss <- sum((tensr::atrans(S, f_d_inv) - S) ^ 2)
  sure_val <- 2 * tau2 * divergence_f + rss - tau2 * prod(p)
  ## gsure_val <- (rss / prod(p)) / (1 - divergence_f / prod(p)) ^ 2
  gsure_val <- rss / (1 - divergence_f / prod(p)) ^ 2
  mean_est <- tensr::atrans(tensr::atrans(S, f_d_inv), hosvd_x$U)
  return(list(sure_val = sure_val, gsure_val = gsure_val,
              mean_est = mean_est))
}

##' Wrapper for \code{\link{get_c}} and \code{\link{sure_given_c}}.
##'
##' @inheritParams get_c
##' @inheritParams sure_given_c
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @seealso \code{\link{get_c}}, \code{\link{sure_given_c}}.
##'
##' @export
sure <- function(X, func, dfunc, lambda, tau2 = 1) {
  ## wrapper for get_c and sure_given_c
  C_output <- get_c(X)
  return(sure_given_c(obj = C_output, func, dfunc, lambda, tau2 = tau2))
}

##' Iterates through all multilinear ranks less than \code{max_nrank} and
##' returns all SUREs.
##'
##'
##' Iterate through all ranks less than max_nrank and choose rank with smallest
##' sure.
##'
##' @param X An array of numerics. The data.
##' @param max_nrank A vector of positive integers. The maximum rank to inspect for each mode.
##' @param tau2 The variance, assumed known. Defaults to 1.
##'
##' @return \code{min_rank} The multilinear rank that minimizes the SURE.
##'
##' \code{min_sure} The minimum SURE value.
##'
##' \code{sure_vec} A vector of all the SURE values of all the multilinear ranks that were inspected.
##'
##' \code{all_ranks} A matrix of all of the multilinear ranks that were inspected.
##'
##' \code{est} The final mean estimate.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##' \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral Estimators}.
##' \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
sure_rank <- function(X, max_nrank = dim(X), tau2 = 1) {
  p <- dim(X)
  n <- length(p)
  func <- list()
  dfunc <- list()
  for (mode_index in 1:n) {
    func[[mode_index]] <- f_truncate
    dfunc[[mode_index]] <- df_truncate
  }
  C_output <- get_c(X)
  all_ranks <- expand.grid(lapply(max_nrank, seq, from = 1))
  sure_vec <- rep(NA, length = prod(max_nrank))
  for (index in 1:prod(max_nrank)) {
    sure_vec[index] <-
      sure_given_c(C_output, func, dfunc,
                   all_ranks[index, ], tau2 = tau2)$sure_val
  }
  which_rank_min <- which.min(sure_vec)



  min_rank <- c(as.matrix(all_ranks[which_rank_min, ]))
  min_sure <- sure_vec[which_rank_min]

  est <- sure_given_c(C_output, func = func, dfunc = dfunc,
                      lambda = min_rank, tau2 = tau2)$mean_est

  return(list(min_rank = min_rank, min_sure = min_sure,
              sure_vec = sure_vec, all_ranks = all_ranks, est = est))
}

##' Runs a stochastic gradient descent to minimize SURE for higher-order
##' spectral estimators.
##'
##' This function is in beta. Not sure if it works properly or is even
##' worthwhile considering that \code{\link{soft_coord}} works pretty well.
##' Also, we have to call \code{\link{get_c}} anyway, so not sure if this function
##' actually reduces the computational complexity.
##'
##' @param c_obj = Output from \code{\link{get_c}}.
##' @param sgd_lambda Lambda value from sgd for lambda from hose.
##' @param sgd_lambda_c Lambda value from sgd for c from hose.
##' @param c_init Initalization for c from hose.
##' @param lambda_init Initialization for lambda from hose.
##' @param sgd_c c value from sgd.
##' @param itermax Maximum number of iterations for sgd.
##' @param tau2 Known variance.
##' @param print_current Should we print the results at each iteration?
##' @param every_iter If \code{print_current = TRUE}, then every
##'   \code{every_iter} iteration will be printed.
##' @param alpha Burnin for final estimates.
##' @param calc_final Should we calculate the final sure and estimates?
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
sgd_given_c <- function(c_obj, sgd_lambda, sgd_lambda_c, c_init = 1,
                        lambda_init = NULL, sgd_c = 1 / 2 + 0.001,
                        itermax = 10000, tau2 = 1, print_current = TRUE,
                        every_iter = 1000, alpha = 0.2, calc_final = TRUE) {
  hosvd_x <- c_obj$hosvd_x
  c_array <- c_obj$C_array
  p <- dim(c_array)
  n <- length(p)

  c_current <- c_init
  if (is.null(lambda_init)) {
    lambda_current <- rep(0, length = n)
  } else {
    lambda_current <- lambda_init
  }

  lambda_mat <- matrix(NA, nrow = itermax, ncol = n)
  c_vec <- rep(NA, length = itermax)
  for (sgd_index in 1:itermax) {
    max_toget <- rep(NA, length = n)
    for (mode_index in 1:n) {
      max_toget[mode_index] <-
        sum(hosvd_x$D[[mode_index]] > lambda_current[mode_index])
    }


    ## use bernoulli to see if I draw a 0 for gradient
    p_1 <- prod(max_toget) / prod(p)
    is_0 <- rbinom(1, 1, 1 - p_1)

    if (is_0 == 1) {
      lambda_new <- lambda_current
      c_new <- c_current
      lambda_mat[sgd_index, ] <- lambda_new
      c_vec[sgd_index] <- c_new
    } else {
      index_current <- c()
      for (mode_index in 1:n) {
        index_current <-
          c(index_current, sample(1:max_toget[mode_index], size = 1))
      }

      sigma_current <- c()
      f_current <- c()
      for (mode_index in 1:n) {
        sigma_current[mode_index] <-
          hosvd_x$D[[mode_index]][index_current[mode_index]]
        f_current[mode_index] <- f_lasso(sigma_current[mode_index],
                                         lambda_current[mode_index])
      }

      f_times_sigma_inv <- f_current / sigma_current
      f_s_mult <- prod(f_times_sigma_inv)
      S2_current <- hosvd_x$S[matrix(index_current, nrow = 1)] ^ 2
      C_array_current <- c_array[matrix(index_current, nrow = 1)]

      common_to_all <- -c_current * f_s_mult ^ 2 * S2_current + f_s_mult *
        S2_current - tau2 * f_s_mult * C_array_current

      inv_temp <-
        matrix(rep(1 / (f_current * sigma_current), n), nrow = n)
      diag(inv_temp) <- 0

      grad_current_lambda <- 2 * c_current *
        (common_to_all - colSums(inv_temp) *
           tau2 * f_s_mult * S2_current) / f_current

      step_size <- sgd_c / (sgd_lambda * sgd_index)

      lambda_new <- lambda_current - step_size * grad_current_lambda


      grad_c <- 2 * c_current * f_s_mult ^ 2 * S2_current -
        2 * f_s_mult * S2_current +
        2 * tau2 * f_s_mult * C_array_current + 2 * tau2 * f_s_mult *
        S2_current * sum(1 / (f_current * sigma_current))

      step_size_c <- sgd_c / (sgd_lambda_c * sgd_index)
      c_new <- c_current - grad_c * step_size_c


      lambda_mat[sgd_index, ] <- lambda_new
      c_vec[sgd_index] <- c_new

      c_current <- c_new
      lambda_current <- lambda_new
    }

    if (print_current & (sgd_index %% every_iter == 0)) {
      cat("Current Lambda = ", lambda_current, "\n")
      cat("     Current c = ", c_current, "\n")
      cat("     Grad of c = ", grad_c, "\n")
      cat("Grad of lambda = ", grad_current_lambda, "\n\n")
    }
  }

  if (calc_final) {
    ## get final estimates
    c_final <- mean(c_vec[round(alpha * itermax):itermax])
    lambda_final <- colMeans(lambda_mat[round(alpha * itermax):itermax, ])
    lambda_list <- list()
    lambda_list[[1]] <- c(lambda_final[1], c_final)
    func_lasso <- list()
    dfunc_lasso <- list()
    func_lasso[[1]] <- f_lasso_mult
    dfunc_lasso[[1]] <- df_lasso_mult
    for (mode_index in 2:n) {
      lambda_list[[mode_index]] <- lambda_final[mode_index]
      func_lasso[[mode_index]] <- f_lasso
      dfunc_lasso[[mode_index]] <- df_lasso
    }
    sure_est <- sure_given_c(c_obj, func = func_lasso, dfunc = dfunc_lasso,
                             lambda = lambda_list, tau2 = tau2)
    return(list(c_vec = c_vec, lambda_mat = lambda_mat, c_final = c_final,
                lambda_final = lambda_final, sure_final = sure_est))
  } else {
    return(list(c_vec = c_vec, lambda_mat = lambda_mat))
  }
}
