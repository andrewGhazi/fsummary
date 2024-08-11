
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

ddf = rnorm(4e3*1000) |> # fake 1k draws of 100 example variables from each of 4 chains
  matrix(ncol = 1000) |> 
  qDT() |> 
  mtt(`.draw` = 1:4e3,
      `.iteration` = rep(1:1000, times = 4),
      `.chain` = rep(1:4, each = 1000)) |> 
  posterior::as_draws_df()

fsummary(ddf) |> # instead of posterior::summarise_draws(ddf)
  head() 
#>    variable     mean   median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>      <char>    <num>    <num> <num> <num> <num> <num> <num>    <num>    <num>
#> 1:       V1 -0.04636 -0.06079 1.013 1.009 -1.73  1.62     1     3871     4000
#> 2:       V2 -0.00991 -0.01376 0.995 0.990 -1.67  1.62     1     3973     3867
#> 3:       V3  0.01591  0.00210 1.007 1.000 -1.65  1.67     1     4138     3972
#> 4:       V4 -0.00985 -0.01078 0.992 0.994 -1.63  1.60     1     4103     3881
#> 5:       V5 -0.00786 -0.00513 0.998 1.021 -1.68  1.58     1     4079     4013
#> 6:       V6  0.00435 -0.00394 1.001 1.013 -1.60  1.69     1     3934     3852
```

On one core, it’s about 2-3 times faster than `summarise_draws()` when
computing convergence metrics and 5-6 times faster without. It also uses
less memory in both cases. A quick test on my machine:

``` r

check_fun = function(x,y) {
  waldo_res = waldo::compare(x, y, tolerance = 1e-6)

  # fsummary returns a data.table instead of a tibble
  non_classattr_diffs = grep(value = TRUE,
                             "class\\(old|attr\\(old",
                             waldo_res,
                             invert = TRUE)

  values_all_near = all(mapply(\(z, w) all((z - w) < 1e-6),
                               x |> num_vars(),
                               y |> num_vars()))

  (length(non_classattr_diffs) == 0) & values_all_near
}

bench::mark(fsummary = {fsummary(ddf, .cores = 1)},
            posterior = {summarise_draws(ddf)},
            iterations = 10,
            check = check_fun,
            filter_gc = FALSE)
```

      expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time 
      <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> 
    1 fsummary      2.24s    2.55s     0.388    1.66GB     5.12    10   132      25.8s 
    2 posterior     6.08s     6.3s     0.159    5.13GB     2.62    10   165      1.05m 

``` r
bench::mark(fsummary = {fsummary(ddf, .cores = 1,
                                 conv_metrics = FALSE)},
            posterior = {summarise_draws(ddf, 
                                         default_summary_measures())},
            iterations = 10,
            check = check_fun,
            filter_gc = FALSE)
```

      expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time 
      <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> 
    1 fsummary    182.7ms 183.54ms     5.33     92.1MB     2.13    10     4      1.87s 
    2 posterior      1.1s    1.11s     0.872     759MB     2.88    10    33     11.47s 

![](man/figures/comparison.png)

# TODO

- parallelization
- faster convergence metrics with better ranking / qnorm (I’ve gotten
  the inverse normal transformation going twice as fast in Julia at
  least…)
