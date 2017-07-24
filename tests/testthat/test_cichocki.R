context("cichocki")

test_that("meet works OK", {
  set.seed(9)
  p <- c(10, 10, 10)
  m <- c(5, 5, 5)
  E <- array(stats::rnorm(prod(p), sd = 1), dim = p)
  S <- array(stats::rnorm(prod(m)), dim = m)
  U <- list(matrix(stats::rnorm(m[1] * p[1]), nrow = p[1]),
            matrix(stats::rnorm(m[2] * p[2]), nrow = p[2]),
            matrix(stats::rnorm(m[3] * p[3]), nrow = p[3]))
  theta <- tensr::atrans(S, U)
  theta <- theta * sqrt(sum(E ^ 2) / sum(theta ^ 2))
  Y <- theta + E
  sout <- score_ylc(Y, return_est = TRUE)
  sout$rank
  dout <- trunc_hosvd(Y = Y, method = "mdl")
  dout$rank
  svaout <- trunc_hosvd(Y = Y, method = "pa")
  svaout$rank
  bcvout <- svaout <- trunc_hosvd(Y = Y, method = "bcv")
  bcvout$rank

  ## test two implementations of MDL
  d <- svd(tensr::mat(Y, 2))$d ^ 2
  t1 <- mdl_eigen(d, rho = 0.0103, I_nbar = 100)
  double_dl(d, r = 2, rho = 0.0103, I_nbar = 100)
  t2 <- mdl_simple(lambda_vec = d, rho = 0.0103, I_nbar = 100)

  expect_equal(t1$dl$dl[10:1], t2$dl$dl)

}
)
