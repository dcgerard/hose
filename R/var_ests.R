#' Estimate of the variance from Gavish and Donoho (2014) based on the
#' Marcenko-Pasture distribution.
#'
#' Returns a normalized squared median singular value that is a
#' consistent estimate of the variance under a specific asymptotic
#' regime.
#'
#' Under the asymptotic regime where the rank of the mean is fixed and
#' the signal to noise ratio is fixed, the squared median singular
#' value divided by the variance and the larger dimension size will
#' converge to the median of the Marcenko-Pasture distribution. Hence,
#' the squared median singular value divided by the larger dimension size
#' and the median of the Marcenko-Pastur distribution will converge to
#' the variance of the data matrix. This estimator works really well
#' when the rank is much less than the dimensions of the matrix and
#' really poorly otherwise.
#'
#'
#' @param dmed A positive numeric. The median singular value of the
#'     data matrix
#' @param N A positive integer. The row dimension.
#' @param p A positive integer. The column dimension.
#'
#' @return \code{sig2_est} A positive numeric. The estimate of the
#'     variance of the data matrix.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Gavish, M., & Donoho,
#'     D. L. (2014). \href{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6846297&tag=1}{The
#'     optimal hard threshold for singular values is
#'     4/sqrt(3)}. Information Theory, IEEE Transactions on, 60(8),
#'     5040-5053.
sig_mp <- function(dmed, N, p) {
    mu <- RMTstat::qmp(0.5, svr = max(N, p) / min(N, p))

    sig2_est <- dmed ^ 2 / (max(N, p) * mu)

    return(sig2_est)
}


#' Simplest estimator for variance.
#'
#' An estimator for the variance of a matrix of independent random
#' variables when assuming the mean matrix is low rank with the rank
#' known.
#'
#' Let \eqn{Y} be a matrix with row dimension \eqn{n} and column
#' dimension \eqn{p} where \eqn{n \ge p}. The model is \eqn{Y =
#' \Theta + \sigma E} where \eqn{\Theta} is low rank and \eqn{E}
#' contains independent elements with mean zero and variance one. The
#' MLE of \eqn{\sigma^2} is the sum of squares of the last \eqn{p - r}
#' singular values divided by \eqn{pn}. This estimator has negative
#' bias so Choi et al (2014) suggested to instead divide by \eqn{n(p -
#' r)}. This is the estimator implemented here. It seems to have
#' reasonable performance when the rank is small.  Though the biggest
#' drawback here is that you have to know the rank ahead of time.
#'
#' @param d A vector of positive numerics. The singular values of the
#'     data matrix.
#' @param N A positive integer. The row dimension.
#' @param p A positive integer. The column dimension.
#' @param r A positive integer. The known rank of the matrix. Must be
#'     less than \code{min(N, p)}
#'
#' @return \code{sig2_est} A positive numeric. The estimate of the
#'     variance.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Choi, Yunjin, Jonathan Taylor, and Robert
#'     Tibshirani. \href{http://arxiv.org/abs/1410.8260}{"Selecting the number of principal components: Estimation of the true rank of a noisy matrix."}
#'     arXiv preprint arXiv:1410.8260 (2014).
#'
sig_naive <- function(d, N, p, r) {
    if (r >= min(N, p)) {
        stop("r needs to be less than min(N, p)")
    }

    sig2_est <- sum(d[(r+1):length(d)] ^ 2) / (max(N, p) * (min(N, p) - r))
    return(sig2_est)
}

#' Soft-impute Cross-validation as described in Choi et al (2014).
#'
#' Hold out some data from a matrix and use \code{softImpute} to complete the
#' matrix. The tuning parameter with the smallest prediction error is selected.
#'
#' @param Y The data matrix.
#' @param k A positive integer. The fold for the soft-impute cross validation.
#'   Default is 10.
#' @param lambda_grid A vector of positive numerics. The values of lambda to
#'   compute. The default is 20 values from the minimum to the maximum singular
#'   value of \code{Y}.
#' @param print_update A logical. Should we print to the screen the status of
#'   the cross-validation-ish procedure at each iteration (TRUE) or not (FALSE)?
#'
#' @return \code{lambda_min} A positive numeric. The lambda that minimizes the
#'   prediction error.
#'
#'   \code{lambda_grid} A vector of positive numerics. The putative lambdas.
#'
#'   \code{pred_err_vec} A vector of positive numerics. The prediction errors
#'   for the lambdas in \code{lambda_grid}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Choi, Yunjin, Jonathan Taylor, and Robert Tibshirani.
#'   \href{http://arxiv.org/abs/1410.8260}{"Selecting the number of principal
#'   components: Estimation of the true rank of a noisy matrix."} arXiv preprint
#'   arXiv:1410.8260 (2014).
#'
#' @seealso \code{\link{sig_soft}} that calls \code{soft_cv} to find the optimal
#'   lambda.
#'
#'
soft_cv <- function(Y, k = 10, lambda_grid = NULL, print_update = FALSE) {

    svdY <- svd(Y)
    d <- svdY$d
    if(is.null(lambda_grid)) {
        lambda_grid <- seq(min(d), max(d), length = 20)
    }

    num_miss <- prod(dim(Y)) / k

    miss_seq <- sample(1:prod(dim(Y)))

    pred_err_vec <- rep(0, length = length(lambda_grid))
    for(index in 1:k) {
        Ymiss <- Y
        which_miss <- miss_seq[((index - 1) * num_miss + 1):(index * num_miss)]
        Ymiss[which_miss] <- NA

        if(print_update) {
            cat("Fold Number:", index, "\n")
        }
        for(lambda_index in 1:length(lambda_grid)) {
            fit1 <- softImpute::softImpute(Ymiss, rank.max = min(dim(Y)) - 1,
                                           lambda = lambda_grid[lambda_index])
            Yest <- softImpute::complete(Ymiss, fit1)
            pred_err_vec[lambda_index] <- pred_err_vec[lambda_index] * (index - 1) / index +
                sum((Y - Yest)^2) / index

            if(print_update) {
                cat("Lambda =", round(lambda_grid[lambda_index], digits = 1), ":",
                    round(pred_err_vec[lambda_index] / max(c(1, pred_err_vec)), digits = 2), "\n")
            }
        }
    }

    lambda_min <- lambda_grid[which.min(pred_err_vec)]

    return(list(lambda_min = lambda_min, lambda_grid = lambda_grid,
                pred_err_vec = pred_err_vec, svdY = svdY))
}

#' Variance estimators from Choi et al (2014).
#'
#' Returns a degrees of freedom corrected residual mean-squared error
#' of a soft-thresholding estimator where the tuning parameter is
#' chosen by a cross-validation-ish procedure.
#'
#' See Choi et al (2014) for details. This seems to be the best
#' estimation procedure so far, but also takes the longest.
#'
#' You can try out different values of \code{c_val} with the outputs
#' of \code{sse} and \code{dfLambda}
#'
#'
#' @param c_val The ad-hoc adjustment to the degrees of freedom. Choi et
#'     al (2014) found that 2/3 worked well in simulations.
#' @inheritParams soft_cv
#'
#' @return \code{sig2_est} A positive numeric. The estimate of the
#'     variance.
#'
#'     \code{sse} A positive numeric. The sum of squared errors for
#'     estimated Y.
#'
#'     \code{dfLambda} A positive integer. The estimated df.
#'
#'
#' @author David Gerard
#'
#' @export
#'
#' @references Choi, Yunjin, Jonathan Taylor, and Robert
#'     Tibshirani. \href{http://arxiv.org/abs/1410.8260}{"Selecting the number of principal components: Estimation of the true rank of a noisy matrix."}
#'     arXiv preprint arXiv:1410.8260 (2014).
#'
#' @seealso \code{\link{soft_cv}} for the cross-validation-ish procedure.
#'
sig_soft <- function(Y, c_val = 2/3, k = 10, lambda_grid = NULL, print_update = FALSE) {
    soft_out <- soft_cv(Y = Y, k = k, lambda_grid = lambda_grid, print_update = print_update)

    svdY <- soft_out$svdY

    dNew <- svdY$d - soft_out$lambda_min
    dNew[dNew < 0] <- 0
    dfLambda <- sum(dNew > 0)

    Yest <- svdY$u %*% diag(dNew) %*% t(svdY$v)

    sse <- sum((Y - Yest) ^ 2)
    sig2_est <-  sse / (max(dim(Y)) * (min(dim(Y)) - c_val * dfLambda))
    return(list(sig2_est = sig2_est, sse = sse, dfLambda = dfLambda))
}

#' Non-parametric variance estimation for tensor datasets.
#'
#' Implements one of three variance estimation procedures for matrix data sets
#' on each matricization of a data tensor. These estimates are then averaged to
#' come up with a final variance estimate.
#'
#' If a data tensor has low multilinear rank mean, then this function will use
#' matrix-specific variance estimation procedures along each matricization of the
#' data tensor to estimate the variance. These estimates are then averaged to
#' provide the final estimate of the variance.
#'
#' @param Y An array of numerics. The data tensor. Assumed to have a mean that
#'   is low multilinear rank.
#' @param sig_method A string. The variance estimation procedure to use on each
#'   matricization of the data tensor. Either "soft" for calling
#'   \code{\link{sig_soft}}, "mp" for calling \code{\link{sig_mp}}, or "naive"
#'   for calling \code{\link{sig_naive}}. See those functions for descriptions
#'   of these procedures.
#' @param r A vector of positive integers. The known multilinear rank of the
#'   mean tensor. Only needed if \code{sig_method == "naive"}.
#'
#' @return \code{sig2_est} A positive numeric. The final variance estimate.
#'
#'   \code{sig2_vec} A vector of positive numerics. The estimates of the
#'   variance when using each matricization of the data tensor. If these numbers
#'   are wildly different, then you should be skeptical of this variance
#'   estimation approach.
#'
#' @references Choi, Yunjin, Jonathan Taylor, and Robert Tibshirani.
#'   \href{http://arxiv.org/abs/1410.8260}{"Selecting the number of principal
#'   components: Estimation of the true rank of a noisy matrix."} arXiv preprint
#'   arXiv:1410.8260 (2014).
#'
#'   Gavish, M., & Donoho, D. L. (2014).
#'   \href{http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=6846297&tag=1}{The
#'    optimal hard threshold for singular values is 4/sqrt(3)}. Information
#'   Theory, IEEE Transactions on, 60(8), 5040-5053.
#'
#' @author David Gerard
#'
#' @export
#'
#' @seealso \code{\link{sig_soft}} for the soft-thresholding procedure,
#'   \code{\link{sig_mp}} for a procedure based on a specific asymptotic
#'   framework, and \code{\link{sig_naive}} for a basic estimator when the
#'   multilinear rank is known.
#'
tensor_var_est <- function(Y, sig_method = c("soft", "mp", "naive"), r = NULL) {
    sig_method <- match.arg(sig_method, c("soft", "mp", "naive"))
    p <- dim(Y)
    n <- length(p)

    sig_vec <- rep(NA, length = n)

    if (sig_method == "soft") {
        for(mode_index in 1:n) {
            Ymat <- tensr::mat(Y, mode_index)
            sig_vec[mode_index] <- sig_soft(Ymat)$sig2_est
            cat("Mode", mode_index, ", Sigma Hat =", sig_vec[mode_index], "\n")
        }
    } else if (sig_method == "mp") {
        for(mode_index in 1:n) {
            Ymat <- tensr::mat(Y, mode_index)
            dmed <- stats::median(svd(Ymat, nu = 0, nv = 0)$d)
            sig_vec[mode_index] <- sig_mp(dmed = dmed, N = dim(Ymat)[1], p = dim(Ymat)[2])
        }
    } else if (sig_method == "naive") {
        if (is.null(r)) {
            stop('r cannot be NULL if sig_method = "naive"')
        }
        for(mode_index in 1:n) {
            Ymat <- tensr::mat(Y, mode_index)
            Ysv <- svd(Ymat, nu = 0, nv = 0)$d
            sig_vec[mode_index] <- sig_naive(d = Ysv, N = dim(Ymat)[1], p = dim(Ymat)[2],
                                             r = r[mode_index])
        }

    } else {
        stop("Invalid sig_method")
    }
    sig2_est <- mean(sig_vec)
    return(list(sig2_est = sig2_est, sig2_vec = sig_vec))
}

