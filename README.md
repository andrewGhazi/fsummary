
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fsummary

<!-- badges: start -->
<!-- badges: end -->

The goal of fsummary is to compute posterior summaries from `cmdstanr`
FAST!

## Installation

You can install the development version of fsummary from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("andrewGhazi/fsummary")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fsummary)
library(collapse)
library(posterior)

options(digits = 3)

ddf = rnorm(4e3*100) |> # fake 1k draws of 100 example variables from each of 4 chains
  matrix(ncol = 100) |> 
  qDT() |> 
  mtt(`.draw` = 1:4e3,
      `.iteration` = rep(1:1000, times = 4),
      `.chain` = rep(1:4, each = 1000)) |> 
  posterior::as_draws_df()

fsummary(ddf) |> # instead of posterior::summarise_draws(ddf)
  head() 
#>    variable     mean    median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>      <char>    <num>     <num> <num> <num> <num> <num> <num>    <num>    <num>
#> 1:       V1  0.00376  0.011541 0.988 0.997 -1.62  1.60     1     4041     3921
#> 2:       V2  0.02068  0.030167 1.014 1.008 -1.67  1.65     1     3698     3948
#> 3:       V3  0.00683 -0.008664 1.019 1.018 -1.68  1.68     1     3830     3739
#> 4:       V4  0.00374  0.000502 1.008 1.013 -1.61  1.65     1     3364     3774
#> 5:       V5 -0.01915 -0.021491 1.009 1.003 -1.69  1.68     1     3962     3892
#> 6:       V6  0.02870  0.000575 0.994 1.011 -1.58  1.68     1     3737     3989
```

On one core, it’s roughly twice as fast as `summarise_draws()` when
computing convergence metrics and ~5x faster without. Example:

``` r
microbenchmark::microbenchmark({summarise_draws(ddf)},
                               {fsummary(ddf, .cores = 1)},
                               times = 3)
```

    Unit: milliseconds
                                  expr min  lq mean median  uq max neval
          {     summarise_draws(ddf) } 666 671  691    677 703 729     3
     {     fsummary(ddf, .cores = 1) } 278 278  283    278 285 292     3

``` r
microbenchmark::microbenchmark({summarise_draws(ddf, default_summary_measures())},
                               {fsummary(ddf, .cores = 1, conv_metrics = FALSE)},
                               times = 3)
```

    Unit: milliseconds
                                                    expr   min    lq  mean median    uq   max neval
    { summarise_draws(ddf, default_summary_measures()) } 114.4 114.5 130.4  114.6 138.3 162.1     3
     { fsummary(ddf, .cores = 1, conv_metrics = FALSE) }  18.6  18.7  19.3   18.9  19.6  20.3     3
