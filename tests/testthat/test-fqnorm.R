test_that("fqnorm works", {

  x1 = runif(1e3)
  x2 = runif(1e3, .999)
  x3 = runif(1e3, max = .001)
  x4 = runif(1e3, min = .45, max = .55)

  p = c(x1, x2, x3, x4)

  expect_equal(fsummary:::fqnorm(p)[,1], qnorm(p),
               tolerance = 1e-10)
})
