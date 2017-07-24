## Functions to implement MEET from Yokota, Lee, and Cichocki (2017) for comparisons.
## https://doi.org/10.1109/TSP.2016.2620965

#' Implements the Modified eivgenvalues estimator for Tucker rank determination (MEET) algorithm
#' from Yokota, Lee, and Cichocki (2017).
#'
#' The model is \deqn{Y = A + E,} where A is a low-rank matrix and E is
#' the noise matrix is i.i.d. N(0,s) entries. s is not assumed to be known.
#'
#' @param Y A numeric array.
#' @param tau2 The known variance. Will be estimated if missing. Not supported yet.
#'
#' @return A list with two elements: \code{rank} the estimated multilinear rank of the tensor and
#'     \code{sigma} a vector of the estimated variance when using each mode.
#'
#'
#' @author David Gerard
#'
#'
#' @references Yokota, Tatsuya, Namgil Lee, and Andrzej Cichocki. "Robust multilinear tensor rank estimation using higher order singular value decomposition and information criteria." IEEE Transactions on Signal Processing 65.5 (2017): 1196-1206. DOI: 10.1109/TSP.2016.2620965
#'
meet <- function(Y, tau2 = NULL) {
  warning("meet is not fully supported yet")
  if (!is.null(tau2)) {
    stop("Not supported yet: tau2 must be NULL.")
  }
  p <- dim(Y)
  y_hosvd <- hosvd_full(Y)

  rank_vec <- rep(NA, length = length(p))
  sigma_vec <- rep(NA, length = length(p))
  for (index in 1:length(p)) {
    I_nbar <- prod(p[-index])
    lambda_vec <- (y_hosvd$D[[index]] ^ 2) / I_nbar


    if (is.null(tau2)) { # estimate variance
      mu_vec <- colSums(tensr::mat(y_hosvd$S, index) ^ 2)
      nu_vec <- sort(mu_vec, decreasing = FALSE)
      ## Get modified eigenvalues
      cum_vec <- cumsum(nu_vec)
      min_ind <- sum(lambda_vec[length(lambda_vec)] * prod(p) > cum_vec)
      rho_est <- (I_nbar- min_ind) / I_nbar
      sigma_est <- 1 / (prod(p) * (1 - rho_est)) * cum_vec[min_ind]
    } else {
      sigma_est <- tau2
      rho_est   <- 1 - lambda_vec[length(lambda_vec)] / sigma_est
    }

    mod_lambda_vec <- lambda_vec / rho_est - (1 - rho_est) / rho_est * sigma_est

    ## MDL using modified eigenvalues
    rout <- mdl_eigen(lambda_vec = mod_lambda_vec, rho = rho_est, I_nbar = I_nbar)
    rank_est <- rout$rank_est
    rank_vec[index] <- rank_est
    sigma_vec[index] <- sigma_est
  }

  ## Deal with special case of estimated 0 rank
  if (any(rank_vec == 0)) {
    rank_vec <- rep(0, length = length(p))
  }

  return(list(rank = rank_vec, sigma = sigma_vec))
}

#' SCORE method from Yokota, Lee, and Cichoki (2017).
#'
#' Implements the "Sparse CORE" method from Yokota, Lee, and Cichocki (2017).
#' The idea is that one uses only a subset of the core array from the HOSVD
#' to calculate the "modified" singular values. Then one uses the minimum
#' description length (see \code{\link{mdl_eigen}}) criteria to choose the rank
#' of the array. The porportion of the core should be small. We have one percent
#' as the default as suggested in their paper, but my guess is that this should
#' actually depend on the rank of the mean tensor.
#'
#' @param Y An array of numerics.
#' @param rho The proportion of the core array to use to calcuate the
#'     mode-specific modified singular values.
#' @param sampmin The minimum number of columns to of the core to use
#'     in calculating the mode-specific modified singular values.
#' @param return_est A logical. Should we return the truncated HOSVD
#'     (\code{TRUE}) or not (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @references Yokota, Tatsuya, Namgil Lee, and Andrzej Cichocki. "Robust multilinear tensor rank estimation using higher order singular value decomposition and information criteria." IEEE Transactions on Signal Processing 65.5 (2017): 1196-1206. DOI: 10.1109/TSP.2016.2620965
#'
#' @export
score_ylc <- function(Y, rho = 0.01, sampmin = min(c(dim(Y), 5)),
                      return_est = TRUE) {
  ## Check input -----------------------------------------------------
  assertthat::are_equal(is.logical(return_est), TRUE)
  assertthat::assert_that(rho < 1, rho > 0)
  assertthat::assert_that(sampmin > 0, sampmin <= min(dim(Y)))

  p <- dim(Y)
  y_hosvd <- hosvd_full(Y)
  rank_vec <- rep(NA, length = length(p))
  for (index in 1:length(p)) {
    S_k2 <- tensr::mat(y_hosvd$S, index) ^ 2
    mu_vec <- colSums(S_k2)
    numcols <- round(rho * ncol(S_k2))
    if (numcols < sampmin) {
      numcols <- sampmin
      mod_rho <- sampmin / ncol(S_k2)
    } else {
      mod_rho <- rho
    }
    which_cols <- order(mu_vec, decreasing = TRUE)[1:numcols]
    mod_lambda <- sort(rowMeans(S_k2[, which_cols]), decreasing = TRUE)
    mdlout <- mdl_eigen(lambda_vec = mod_lambda, rho = mod_rho, I_nbar = ncol(S_k2))
    rank_vec[index] <- mdlout$rank
  }
  return_list <- list()
  return_list$rank = rank_vec
  if (return_est) {
    est <- array(y_hosvd$S[as.matrix(expand.grid(lapply(rank_vec, FUN = seq, from = 1))), drop = FALSE] , dim = rank_vec)
    for (index in 1:length(p)) {
      est <- tensr::amprod(est, y_hosvd$U[[index]][, 1:rank_vec[index], drop = FALSE], index)
    }
  }
  return(return_list)
}

#' Returns the minimum description length given the eigenvalues.
#'
#' @param lambda_vec A numeric vector of eigenvalues.
#' @param rho The estimated rho parameter from Yokota, Lee, and Cichocki (2017)
#' @param I_nbar The product of the other dimensions.
#'
#' @author David Gerard
#'
#' @references Yokota, Tatsuya, Namgil Lee, and Andrzej Cichocki. "Robust multilinear tensor rank estimation using higher order singular value decomposition and information criteria." IEEE Transactions on Signal Processing 65.5 (2017): 1196-1206. DOI: 10.1109/TSP.2016.2620965
#'
#' @seealso \code{\link{meet}}
#'
#' @return A list with two elements: \code{rank_est} the estimated rank; and \code{dl} a data-frame
#'     with the description lengths for each rank from \code{0} to \code{length(lambda_vec) - 1}.
#'
mdl_eigen <- function(lambda_vec, rho, I_nbar) {

  assertthat::assert_that(length(lambda_vec) > 0)
  assertthat::are_equal(order(lambda_vec), length(lambda_vec):1)
  lambda_vec <- lambda_vec[length(lambda_vec):1] ## reverse for easier handeling
  n_lambda <- length(lambda_vec)

  r_vec <- (n_lambda - 1):0
  numer_vec <- cumprod(lambda_vec) ^ (1 / (n_lambda - r_vec))
  denom_vec <- cumsum(lambda_vec) / (n_lambda - r_vec)
  mdl_vec <- -2 * log(numer_vec / denom_vec) * (rho * I_nbar * (n_lambda - r_vec)) +
    r_vec * (2 * n_lambda - r_vec) * log(rho * I_nbar)

  ## full rank description length?
  ## mdl_vec <- c(mdl_vec, n_lambda ^ 2 * log(rho * I_nbar))

  rank_est <- r_vec[which.min(mdl_vec)]
  return(list(rank_est = rank_est, dl = data.frame(rank = r_vec, dl = mdl_vec)))
}

mdl_simple <- function(lambda_vec, rho, I_nbar) {
  I_n <- length(lambda_vec)
  dl_vec <- rep(NA, length = I_n)
  for (r in 0:(I_n - 1)) {
    dl_vec[r + 1] <- double_dl(lambda_vec = lambda_vec,
                               r = r,
                               rho = rho,
                               I_nbar = I_nbar)
  }
  rank_est <- which.min(dl_vec) - 1
  return(list(rank_est = rank_est, dl = data.frame(rank = 0:(I_n - 1), dl = dl_vec)))
}

double_dl <- function(lambda_vec, r, rho, I_nbar) {
  I_n <- length(lambda_vec)
  sublambda <- lambda_vec[(r + 1):I_n]
  numer <- prod(sublambda ^ (1 / (I_n - r)))
  denom <- mean(sublambda)
  dl <- -2 * log(numer / denom) * (rho * I_nbar * (I_n - r)) +
    r * (2 * I_n - r) * log(rho * I_nbar)
  return(dl)
}
