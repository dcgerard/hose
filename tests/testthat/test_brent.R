context("Brent's Method Updates")

test_that("update_lambda_brent works", {
  set.seed(3)

  theta <- readRDS(file = "low_rank_array.Rds")
  p <- dim(theta)
  Y <- theta + array(rnorm(prod(p)), dim = p)

  tau2 <- 1
  c_obj <- get_c(Y)
  lambda_current <- c(3, 3, 3)
  c_current <- 1

  for (index in 1:10) {
    lambda_current[1] <- update_lambda_brent(c_obj = c_obj, lambda_current = lambda_current,
                                             c_current = c_current, k = 1, tau2 = tau2)[[1]]
    lambda_current[2] <- update_lambda_brent(c_obj = c_obj, lambda_current = lambda_current,
                                             c_current = c_current, k = 2, tau2 = tau2)[[1]]
    lambda_current[3] <- update_lambda_brent(c_obj = c_obj, lambda_current = lambda_current,
                                             c_current = c_current, k = 3, tau2 = tau2)[[1]]
  }

  lambda_current <- c(1, 1, 1)
  sout <- soft_coord(c_obj = c_obj, lambda_init = lambda_current, c_init = 1)
}
)
