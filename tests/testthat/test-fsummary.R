test_that("same answer as posterior::summarise_draws()", {
  fsumm_res = fsummary(test_ddf)
  expect_equal(fsumm_res,
               test_res,
               ignore_attr = TRUE, tolerance = 1e-6)
})
