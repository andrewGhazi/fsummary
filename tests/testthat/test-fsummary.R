test_that("same answer as posterior::summarise_draws()", {

  fsumm_res = fsummary(test_ddf)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

test_that("same answer as posterior::summarise_draws() with fft", {

  fsumm_res = fsummary(test_ddf, fft_acov = TRUE)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

test_that("same answer on big result", {

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
  mirai::daemons(0)
  mirai::daemons(6)
  fsumm_res = fsummary(test_ddf)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

test_that("same answer on big result, parallel", {

  mirai::daemons(0)
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

test_that("poorly mixed chains", {

  mirai::daemons(0)

  bad_conv_fsum = bad_conv_ddf |>
    fsummary::fsummary()

  expect_equal(bad_conv_post, bad_conv_fsum,
               tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("zero variance case", {
  fsumm_res = fsummary(na_ddf)
  expect_equal(fsumm_res,
               na_res,
               ignore_attr = TRUE,
               tolerance = 1e-6)
})

