test_that("data sim works", {


  n = 4000
  d = 1000

  set.seed(123)

  datM = runif(n*d) |>
    matrix(ncol = d)

  expect_lt(datM[1,1] - 0.2875775201246142,
            1e-8)
})
