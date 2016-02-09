##' Higher-order SVD using same signed eigenvectors as matix SVD's.
##'
##' Calculates the left singular vectors of each matrix unfolding of an array,
##' then calculates the core array. The resulting output is a Tucker
##' decomposition.
##'
##' \code{Y} is equal to \code{atrans(S, U)}, up to numerical accuracy.
##'
##' This function differs
##' from the \code{hosvd} function in the package \code{tensr} only in  (1) the sign
##' conditions on the core array and (2) it will also return the mode-specific singular values.
##'
##  More details on the HOSVD can be found in
##' \href{http://epubs.siam.org/doi/abs/10.1137/S0895479896305696}{ De Lathauwer
##' et. al. (2000)}.
##'
##' @references De Lathauwer, L., De Moor, B., & Vandewalle, J. (2000).
##'   \href{http://epubs.siam.org/doi/abs/10.1137/S0895479896305696}{A
##'   multilinear singular value decomposition}. \emph{SIAM journal on Matrix
##'   Analysis and Applications}, 21(4), 1253-1278.
##'
##' @return \code{U} A list of matrices with orthonormal columns. Each matrix
##'   contains the mode-specific singular vectors of its mode.
##'
##'   \code{D} A list of vectors of numerics. The \eqn{k}th vector contains the mode-specific singular
##'   values of the \eqn{k}th matricization of \code{Y}
##'
##'   \code{S}An all-orthogonal array. This is the core array from the HOSVD.
##'
##'
##' @author David Gerard.
##'
##' @export
hosvd_full <- function(Y) {
    ## this function is required for sure() below higher order svd
    m <- dim(Y)
    K <- length(m)

    ## get standard tsvd
    U <- list()
    D <- list()
    D_inv <- list()
    for (k in 1:K) {
        X_k_svd <- svd(tensr::mat(Y, k))
        U[[k]] <- X_k_svd$u
        D[[k]] <- X_k_svd$d
        D_inv[[k]] <- 1 / X_k_svd$d
    }
    S <- tensr::atrans(Y, lapply(U, t))
    return(list(U = U, D = D, S = S))
}

##' Sums elements in a list.
##'
##' @param D A list of summable elements.
##'
##' @return \code{list_sum} The sum of the elements in \code{D}.
##'
##' @author David Gerard.
listSum <- function(D) {
    list_sum <- D[[1]]
    for (list_index in 2:length(D)) {
        list_sum <- list_sum + D[[list_index]]
    }
    return(list_sum)
}

## mode specific functions------------------------------------------

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

##' Positive part function.
##'
##' Returns a vector whose elements are the positive parts of the elements of
##' \code{x}.
##'
##' @param x A vector of numerics.
##'
##' @author David Gerard.
pos_part <- function(x) {
    return(sapply(x, max, 0))
}

##' Soft thresholding shrinkage function. Same as lasso for spectral shrinkage.
##'
##' @param A vector of numerics.
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
##' @seealso \code{\link{f_lass}}.
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

############################### These functions shrink S directly

##' Positive part function.
##'
##' Returns an object whose elements are the positive parts of the elements of
##' \code{X}.
##'
##' @param X A vector, matrix, or array of numerics.
##'
##'
##' @author David Gerard.
pos_part_2 <- function(X) {
    X[X < 0] <- 0
    return(X)
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

############### Competitors

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
##' @param p A vector of positive integers. The dimension of the data array.
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
##' @param lamda_init Initialization for lambda from hose.
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



## Block SURE methods this first one overlaps all pixels and circles back

##' Adds a circular padding after matrix.
##'
##' @param X A matrix of numerics.
##' @param padsize A positive integer. The number of columns/rows to add to
##'   \code{X}.
##'
##' @author David Gerard.
##'
##' @export
pad_circ <- function(X, padsize) {
    p <- dim(X)
    Y <- matrix(NA, nrow = p[1] + padsize, ncol = p[2] + padsize)
    Y[1:p[1], 1:p[2]] <- X
    Y[1:p[1], (p[2] + 1):(p[2] + padsize)] <- X[1:p[1], 1:padsize]
    Y[(p[1] + 1):(p[1] + padsize), 1:p[2]] <- X[1:padsize, 1:p[2]]
    Y[(p[1] + 1):(p[1] + padsize), (p[2] + 1):(p[2] + padsize)] <-
      X[1:padsize, 1:padsize]
    return(Y)
}

##' Calculates the necessary components to calculate the SURE for estimators
##' that block shrink overlapping blocks.
##'
##' @param current_tensor An array of numerics.
##' @param block_size The size of the overlapping blocks.
##' @param p The dimension of \code{current_tensor}.
##'
##' @return A list of elements used in \code{\link{block_sure_given_c_list}} to
##'   calculate the SURE.
##'
##' @seealso \code{\link{block_sure_given_c_list}}, \code{\link{block_sure}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
get_c_list <- function(current_tensor, block_size = dim(current_tensor)[3]) {
    p <- dim(current_tensor)
    ##n <- length(p) ## may not be used.

    tensor_expanded <- array(NA, dim = c(p[1] + block_size - 1,
                                         p[2] + block_size - 1, p[3]))
    for (frame_index in 1:p[3]) {
        mat_expanded <- pad_circ(current_tensor[, , frame_index],
                                 padsize = block_size - 1)
        tensor_expanded[, , frame_index] <- mat_expanded
    }

    c_list <- list()
    block_list <- list()
    overall_index <- 1
    for (x_index in 1:p[1]) {
        for (y_index in 1:p[2]) {
            block_list[[overall_index]] <-
              tensor_expanded[x_index:(block_size + x_index - 1),
                              y_index:(block_size + y_index - 1), ]
            c_list[[overall_index]] <- get_c(block_list[[overall_index]])
            overall_index <- overall_index + 1
        }
    }
    return(list(c_list = c_list, block_list = block_list,
                tensor_expanded = tensor_expanded,
                current_tensor = current_tensor, block_size = block_size,
                p = p))
}

##' Calculates the sure of estimators that shrink subtensors using the output of
##' \code{\link{get_c_list}}.
##'
##' @param c_obj The output from \code{\link{get_c_list}}.
##' @param func A list of functions to apply to each mode.
##' @param dfunc A list of derivatives of the function in \code{func}.
##' @param lambda_list A list of parameter values for \code{func} and
##'   \code{dfunc}.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @return \code{sure_final} A numeric. The SURE value.
##'
##' \code{mean_est} An array of numerics. The mean estimate.
##'
##' @seealso \code{\link{get_c_list}}, \code{\link{block_sure}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
block_sure_given_c_list <- function(c_obj, func, dfunc, lambda_list,
                                    tau2 = 1) {
    tensor_expanded <- c_obj$tensor_expanded
    c_list <- c_obj$c_list
    ##block_list <- c_obj$block_list ## may not be used
    block_size <- c_obj$block_size
    current_tensor <- c_obj$current_tensor
    p <- c_obj$p
    ##n <- length(p) ## may not be used.

    p2 <- dim(tensor_expanded)

    overall_index <- 1
    diverge_tot <- 0
    est_tot1 <- array(0, dim = p2)
    for (x_index in 1:p[1]) {
        for (y_index in 1:p[2]) {
            div_current <- diverge_given_c(c_list[[overall_index]],
                                           func = func, dfunc = dfunc,
                                           lambda = lambda_list)
            diverge_tot <- diverge_tot + div_current$divergence_f
            est_tot1[x_index:(block_size + x_index - 1),
                     y_index:(block_size + y_index - 1), ] <-
              est_tot1[x_index:(block_size + x_index - 1),
                       y_index:(block_size + y_index - 1), ] +
              div_current$mean_est
            overall_index <- overall_index + 1
        }
    }

    est_tot <- est_tot1[1:p[1], 1:p[2], ]

    est_tot[1:(block_size - 1), 1:(block_size - 1), ] <-
      est_tot[1:(block_size - 1), 1:(block_size - 1), ] +
      est_tot1[(p[1] + 1):(p[1] + block_size - 1),
               (p[2] + 1):(p[2] + block_size - 1), ]
    est_tot[1:(block_size - 1), 1:p[2], ] <-
      est_tot[1:(block_size - 1), 1:p[2], ] +
      est_tot1[(p[1] + 1):(p[1] + block_size - 1), 1:p[2], ]
    est_tot[1:p[1], 1:(block_size - 1), ] <-
      est_tot[1:p[1], 1:(block_size - 1), ] +
      est_tot1[1:p[1], (p[2] + 1):(p[2] + block_size - 1), ]

    ## because reusing each pixel block_size^2 times.
    est_tot <- est_tot / block_size ^ 2
    diverge_tot <- diverge_tot / block_size ^ 2

    sure_final <- sum( (est_tot - current_tensor) ^ 2) - tau2 * prod(p) +
      2 * tau2 * diverge_tot
    return(list(sure_final = sure_final, mean_est = est_tot))
}

##' Wrapper for \code{\link{get_c_list}} and
##' \code{\link{block_sure_given_c_list}}.
##'
##' @inheritParams get_c_list
##' @inheritParams block_sure_given_c_list
##'
##' @return \code{sure_final} A numeric. The SURE value.
##'
##' \code{mean_est} An array of numerics. The mean estimate.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##' \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral Estimators}.
##' \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
block_sure <- function(current_tensor, func, dfunc, lambda_list, tau2 = 1,
                       block_size = dim(current_tensor)[3]) {
    c_obj <- get_c_list(current_tensor, block_size = block_size)
    block_sure_given_c_list(c_obj, func = func, dfunc = dfunc,
                            lambda_list = lambda_list, tau2)
}


###### now for block SURE with non-overlapping blocks doesn't work if p_k %%
###### block_size = 1 for either k = 1 or 2

##' Calculates the necessary components to calculate the SURE for estimators
##' that shrink non-overlapping sub-tensors.
##'
##' @param current_tensor An array of numerics.
##' @param block_size The size of the non-overlapping blocks.
##'
##' @return A list of elements used in \code{\link{block_sure_given_c_list_non}}
##'   to calculate the SURE.
##'
##' @seealso \code{\link{block_sure_given_c_list_non}},
##'   \code{\link{block_sure_non}}
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
get_c_list_non <- function(current_tensor,
                           block_size = dim(current_tensor)[3]) {
    ## same as get_c_list but with non-overlapping blocks
    p <- dim(current_tensor)
    ##n <- length(p) ## may not be used.

    c_list <- list()
    block_list <- list()
    overall_index <- 1

    x_max <- floor(p[1] / block_size)
    y_max <- floor(p[2] / block_size)
    for (x_index in 1:(x_max + 1)) {
        for (y_index in 1:(y_max + 1)) {
            if (y_index < (y_max + 1) & x_index < (x_max + 1)) {
                block_list[[overall_index]] <-
                  current_tensor[( (x_index - 1) * block_size + 1):(x_index * block_size),
                                 ( (y_index - 1) * block_size + 1):(y_index * block_size), ]
                c_list[[overall_index]] <- get_c(block_list[[overall_index]])
                overall_index <- overall_index + 1
            } else if (y_index < (y_max + 1) & x_index == (x_max + 1)) {
                block_list[[overall_index]] <-
                  current_tensor[( (x_index - 1) * block_size + 1):p[1],
                                 ( (y_index - 1) * block_size + 1):(y_index * block_size), ]
                c_list[[overall_index]] <- get_c(block_list[[overall_index]])
                overall_index <- overall_index + 1
            } else if (y_index == (y_max + 1) & x_index < (x_max + 1)) {
                block_list[[overall_index]] <-
                  current_tensor[( (x_index - 1) * block_size + 1):(x_index * block_size),
                                 ( (y_index - 1) * block_size + 1):p[2], ]
                c_list[[overall_index]] <- get_c(block_list[[overall_index]])
                overall_index <- overall_index + 1
            } else {
                block_list[[overall_index]] <-
                  current_tensor[( (x_index - 1) * block_size + 1):p[1],
                                 ( (y_index - 1) * block_size + 1):p[2], ]
                c_list[[overall_index]] <- get_c(block_list[[overall_index]])
                overall_index <- overall_index + 1
            }
        }
    }
    return(list(c_list = c_list, block_list = block_list,
                current_tensor = current_tensor, block_size = block_size))
}

##' Calculates the SURE of estimators that shrink non-overlapping sub-tensors
##' using the output of \code{\link{get_c_list_non}}.
##'
##' @param c_obj The output from \code{\link{get_c_list}}.
##' @param func A list of functions to apply to each mode.
##' @param dfunc A list of derivatives of the function in \code{func}.
##' @param lambda_list A list of parameter values for \code{func} and
##'   \code{dfunc}.
##' @param tau2 A positive numeric. The variance. Assumed known and defaults to
##'   1.
##'
##' @return \code{sure_final} A numeric. The SURE value.
##'
##'   \code{mean_est} An array of numerics. The mean estimate.
##'
##' @seealso \code{\link{get_c_list_non}}, \code{\link{block_sure_non}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
block_sure_given_c_list_non <- function(c_obj, func, dfunc, lambda_list,
                                        tau2 = 1) {
    current_tensor <- c_obj$current_tensor
    c_list <- c_obj$c_list
    ##block_list <- c_obj$block_list ## may not be used.
    block_size <- c_obj$block_size
    p <- dim(current_tensor)
    ##n <- length(p) ## may not be used.

    overall_index <- 1
    diverge_tot <- 0
    est_tot <- array(0, dim = p)
    x_max <- floor(p[1] / block_size)
    y_max <- floor(p[2] / block_size)
    for (x_index in 1:(x_max + 1)) {
        for (y_index in 1:(y_max + 1)) {
            if (y_index < (y_max + 1) & x_index < (x_max + 1)) {
                div_current <- diverge_given_c(c_list[[overall_index]],
                                               func = func, dfunc = dfunc,
                                               lambda = lambda_list)
                diverge_tot <- diverge_tot + div_current$divergence_f
                est_tot[( (x_index - 1) * block_size + 1):(x_index * block_size),
                        ( (y_index - 1) * block_size + 1):(y_index * block_size), ] <-
                  div_current$mean_est
                overall_index <- overall_index + 1
            } else if (y_index < (y_max + 1) & x_index == (x_max + 1)) {
                div_current <- diverge_given_c(c_list[[overall_index]],
                                               func = func, dfunc = dfunc,
                                               lambda = lambda_list)
                diverge_tot <- diverge_tot + div_current$divergence_f
                est_tot[( (x_index - 1) * block_size + 1):p[1],
                        ( (y_index - 1) * block_size + 1):(y_index * block_size), ] <-
                  div_current$mean_est
                overall_index <- overall_index + 1
            } else if (y_index == (y_max + 1) & x_index < (x_max + 1)) {
                div_current <- diverge_given_c(c_list[[overall_index]],
                                               func = func, dfunc = dfunc,
                                               lambda = lambda_list)
                diverge_tot <- diverge_tot + div_current$divergence_f
                est_tot[( (x_index - 1) * block_size + 1):(x_index * block_size),
                        ( (y_index - 1) * block_size + 1):p[2], ] <-
                  div_current$mean_est
                overall_index <- overall_index + 1
            } else {
                div_current <-
                  diverge_given_c(c_list[[overall_index]],
                                  func = func, dfunc = dfunc,
                                  lambda = lambda_list)
                diverge_tot <- diverge_tot + div_current$divergence_f
                est_tot[( (x_index - 1) * block_size + 1):p[1],
                        ( (y_index - 1) * block_size + 1):p[2], ] <-
                  div_current$mean_est
                overall_index <- overall_index + 1
            }
        }
    }

    sure_final <- sum( (est_tot - current_tensor) ^ 2) - tau2 * prod(p) +
      2 * tau2 * diverge_tot
    return(list(sure_final = sure_final, mean_est = est_tot))
}

##' Wrapper for \code{\link{get_c_list_non}} and
##' \code{\link{block_sure_given_c_list_non}}.
##'
##' @inheritParams get_c_list_non
##' @inheritParams block_sure_given_c_list_non
##'
##' @return \code{sure_final} A numeric. The SURE value.
##'
##'   \code{mean_est} An array of numerics. The mean estimate.
##'
##' @seealso \code{\link{get_c_list_non}},
##'   \code{\link{block_sure_given_c_list_non}}.
##'
##' @author David Gerard.
##'
##' @references Gerard, D., & Hoff, P. (2015).
##'   \href{http://arxiv.org/abs/1505.02114}{Adaptive Higher-order Spectral
##'   Estimators}. \emph{arXiv preprint} arXiv:1505.02114.
##'
##' @export
block_sure_non <- function(current_tensor, func, dfunc, lambda_list, tau2 = 1,
                           block_size = dim(current_tensor)[3]) {
    c_obj <- get_c_list_non(current_tensor, block_size = block_size)
    block_sure_given_c_list_non(c_obj, func = func, dfunc = dfunc,
                                lambda_list = lambda_list, tau2)
}


## Coordinate Descent for Soft Thresholding ---------------------------

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
        for (k in 1:n) {
            lambda_current[k] <-
              update_lambda(c_obj, lambda_current, c_current, k, tau2)
        }
        c_current <- update_c(c_obj, lambda_current, c_current, tau2)

        if (use_sure) {
            old_sure <- current_sure
            lambda_list <- list()
            lambda_list[[1]] <- c(lambda_current[1], c_current)
            for (mode_index in 2:n) {
                lambda_list[[mode_index]] <- lambda_current[mode_index]
            }
            current_sure <-
              sure_given_c(c_obj, func_lasso, dfunc_lasso,
                           lambda = lambda_list, tau2 = tau2)$sure_val
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
