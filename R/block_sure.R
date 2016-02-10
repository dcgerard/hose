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
