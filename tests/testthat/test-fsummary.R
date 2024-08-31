test_that("same answer as posterior::summarise_draws()", {

  fsumm_res = fsummary(test_ddf)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

test_that("same answer on big result", {

  mirai::daemons(1)

  n = 4000
  d = 1000

  set.seed(123)

  datM = runif(n*d) |>
    matrix(ncol = d)

  expect_lt(datM[1,1] - 0.2875775201246142,
            1e-8)

  big_ddf = datM |>
    qDT() |>
    mtt(`.chain` = rep(1:4, each = 1000),
        `.iteration` = rep(1:1000, times = 4),
        `.draw` = 1:4000)

  class(big_ddf) = c("draws_df",   class(big_ddf))

  big_fsumm_res = fsummary(big_ddf)

  expect_equal(big_fsumm_res,
               big_test_res,
               tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("same answer as posterior::summarise_draws(), parallel", {
  mirai::daemons(6)
  fsumm_res = fsummary(test_ddf)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

test_that("same answer on big result, parallel", {

  mirai::daemons(6)

  n = 4000
  d = 1000

  set.seed(123)

  datM = runif(n*d) |>
    matrix(ncol = d)

  expect_lt(datM[1,1] - 0.2875775201246142,
            1e-8)

  big_ddf = datM |>
    qDT() |>
    mtt(`.chain` = rep(1:4, each = 1000),
        `.iteration` = rep(1:1000, times = 4),
        `.draw` = 1:4000)

  class(big_ddf) = c("draws_df",   class(big_ddf))

  big_fsumm_res = fsummary(big_ddf)

  expect_equal(big_fsumm_res,
               big_test_res,
               tolerance = 1e-6,
               ignore_attr = TRUE)
})
