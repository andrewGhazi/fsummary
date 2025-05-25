
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fsummary

<!-- badges: start -->

<!-- badges: end -->

The goal of fsummary is to compute posterior summaries of `draws_df`
objects - FAST!

## Installation

You can install the development version of fsummary from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("andrewGhazi/fsummary")
```

## Example

This package’s main function, `fsummary()`, computes the same thing as
`posterior::summarise_draws()`. Here we run it on on a simulated draws
data frame `ddf`:

``` r
library(fsummary)
library(collapse)
library(posterior)

options(digits = 3)

set.seed(123)

n_chain = 4
n_iter = 1000
n_var = 1000

ddf = rnorm(n_chain * n_iter * n_var) |> 
  matrix(ncol = 1000) |> 
  qDT() |> 
  mtt(`.draw` = 1:(n_iter*n_chain),
      `.iteration` = rep(1:n_iter, times = n_chain),
      `.chain` = rep(1:n_chain, each = n_iter)) |> 
  posterior::as_draws_df()

fsummary(ddf) |> head()
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
computing convergence metrics and 5-6 times faster without. A couple
quick tests on my machine:

``` r

check_fun = function(x,y) {
  check_res = waldo::compare(x, y, 
                             tolerance = 1e-6, ignore_attr = TRUE)
  length(check_res) == 0
}

bench::mark(fsummary = {fsummary(ddf)},
            posterior = {summarise_draws(ddf)},
            iterations = 10,
            check = check_fun,
            filter_gc = FALSE)
```

      expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
      <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
    1 fsummary      1.64s    1.67s     0.581    1.05GB     5.11    10    88     17.23s
    2 posterior     6.24s    6.44s     0.156    5.13GB     3.36    10   215      1.07m

``` r
bench::mark(fsummary = {fsummary(ddf,
                                 conv_metrics = FALSE)},
            posterior = {summarise_draws(ddf, 
                                         default_summary_measures())},
            iterations = 10,
            check = check_fun,
            filter_gc = FALSE)
```

      expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
      <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
    1 fsummary   188.89ms 194.65ms     4.50     92.1MB     2.70    10     6      2.22s
    2 posterior     1.08s    1.11s     0.871     759MB     4.10    10    47     11.47s

![](man/figures/comparison.png)

You can set up daemons with
[`mirai`](https://shikokuchuo.net/mirai/index.html) for parallelization:

``` r
mirai::daemons(4)
system.time({fsummary(ddf)})
```

       user  system elapsed 
      0.018   0.012   1.227 

(This helps more once you have summaries that take longer than a
second!)

# TODO

- ✔ ~~parallelization with `mirai`~~
  - Something broke. Fix it.
- faster convergence metrics with better ✔ ~~ranking~~ / ~~qnorm~~
- ✔ ~~Option for FFT autocovariance if user knows they have badly
  converged parameters~~
- Submission to CRAN. If you use this package and it would be helpful to
  you for it to be available through CRAN, submit an issue saying so.
  I’ve basically got it ready to submit, it’s just not clear to me
  whether it would be worth doing the paperwork.
