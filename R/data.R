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
