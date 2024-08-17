## code to prepare `posterior_res` dataset goes here

library(fastverse)

set.seed(123)

n_iter = 50 # 101 to make sure the test hits the padding fn
n_chain = 4

test_ddf = rnorm(n_iter*n_chain*3) |>
  matrix(ncol = 3) |>
  qDT() |>
  mtt(`.chain` = rep(1:4, each = n_iter),
      `.iteration`= rep(1:n_iter, times = n_chain),
      `.draw` = 1:(n_iter*n_chain)) |>
  posterior::as_draws_df()

test_res = test_ddf |>
  posterior::summarise_draws()

usethis::use_data(test_ddf, overwrite = TRUE)
usethis::use_data(test_res, overwrite = TRUE)
