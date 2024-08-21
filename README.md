
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

set.seed(123)
ddf = rnorm(4e3*1000) |> # fake 1k draws of 100 example variables from each of 4 chains
  matrix(ncol = 1000) |> 
  qDT() |> 
  mtt(`.draw` = 1:4e3,
      `.iteration` = rep(1:1000, times = 4),
      `.chain` = rep(1:4, each = 1000)) |> 
  posterior::as_draws_df()

fsummary(ddf) |> # instead of posterior::summarise_draws(ddf)
  head() 
#>    variable     mean    median    sd   mad    q5   q95  rhat ess_bulk ess_tail
#>      <char>    <num>     <num> <num> <num> <num> <num> <num>    <num>    <num>
#> 1:       V1  0.00733 -0.004966 0.993 0.975 -1.65  1.64     1     4161     3930
#> 2:       V2 -0.00628 -0.007577 0.996 1.007 -1.63  1.60     1     4099     3978
#> 3:       V3  0.00578 -0.010096 1.006 1.021 -1.65  1.66     1     4005     3807
#> 4:       V4 -0.00163 -0.006376 1.001 0.983 -1.64  1.68     1     3851     3729
#> 5:       V5 -0.03390 -0.041806 1.004 1.008 -1.66  1.63     1     3878     3891
#> 6:       V6 -0.00366 -0.000109 1.000 1.013 -1.66  1.66     1     3547     3824
```

On one core, it’s about 3 times faster than `summarise_draws()` when
computing convergence metrics and 5-6 times faster without. It also uses
less memory in both cases. A quick test on my machine:

``` r

check_fun = function(x,y) {
  check_res = waldo::compare(x, y, tolerance = 1e-6, ignore_attr = TRUE)
  length(check_res) == 0
}

bench::mark(fsummary = {fsummary(ddf, .cores = 1)},
            posterior = {summarise_draws(ddf)},
            iterations = 10,
            check = check_fun,
            filter_gc = FALSE)
```


      expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
      <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
    1 fsummary       1.9s  1.98s     0.485  996.81MB     4.51    10    93     20.61s
    2 posterior     6.09s  6.21s     0.161    5.13GB     4.29    10   267      1.04m

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
    1 fsummary   192.84ms 199.21ms     4.97     92.1MB     2.48    10     5      2.01s
    2 posterior     1.11s    1.16s     0.842     759MB     4.80    10    57     11.88s

![](man/figures/comparison.png)

# TODO

- parallelization with `mirai`
- faster convergence metrics with better ranking / qnorm (I’ve gotten
  the inverse normal transformation going twice as fast in Julia at
  least…)
