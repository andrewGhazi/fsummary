#' Test draws data frame
#'
#' A test draws data frame.
#' @format ## `test_ddf`
#' A draws_df of three parameters, V1, V2, and V3, all numeric with 100 iterations and four chains.
"test_ddf"

#' Test summary result
#'
#' A test summary result to check that fsummary gets the same numbers.
#' @format ## `test_res`
#' A tibble from posterior::summarise_draws(test_ddf)
"test_res"

#' Test summary result
#'
#' A big test summary result to check that fsummary gets the same numbers. No accompanying
#' ddf, see the generating script in data-raw/
#' @format ## `big_test_res` A tibble from posterior::summarise_draws(test_ddf)
"big_test_res"

#' zero variance draws data frame
#'
#' A test draws data frame where some parameters have 0 variance, which should produce some NAs.
#' @format ## `na_ddf`
#' A draws_df of three parameters, V1, V2, and V3, all numeric with 100 iterations and four chains.
"na_ddf"

#' zero variance summary result
#'
#' A test summary result where some parameters have 0 variance, which should produce some NAs.
#' @format ## `na_res`
#' A tibble from posterior::summarise_draws(test_ddf)
"na_res"
