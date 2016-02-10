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
##' @param Y An array of numerics.
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
##'   \code{S} An all-orthogonal array. This is the core array from the HOSVD.
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
