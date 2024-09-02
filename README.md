
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

fsummary(ddf) 
#>       variable     mean   median    sd   mad    q5   q95  rhat ess_bulk
#>         <char>    <num>    <num> <num> <num> <num> <num> <num>    <num>
#>    1:       V1  0.00733 -0.00497 0.993 0.975 -1.65  1.64     1     4161
#>    2:       V2 -0.00628 -0.00758 0.996 1.007 -1.63  1.60     1     4099
#>    3:       V3  0.00578 -0.01010 1.006 1.021 -1.65  1.66     1     4005
#>    4:       V4 -0.00163 -0.00638 1.001 0.983 -1.64  1.68     1     3851
#>    5:       V5 -0.03390 -0.04181 1.004 1.008 -1.66  1.63     1     3878
#>   ---                                                                  
#>  996:     V996 -0.02662 -0.03425 0.999 1.007 -1.65  1.58     1     3993
#>  997:     V997 -0.01837 -0.02869 0.993 1.011 -1.60  1.60     1     3944
#>  998:     V998  0.00392  0.00214 1.009 1.047 -1.61  1.67     1     3879
#>  999:     V999  0.00600  0.01776 1.006 0.984 -1.67  1.68     1     3991
#> 1000:    V1000  0.04056  0.04775 0.988 0.999 -1.60  1.63     1     4175
#>       ess_tail
#>          <num>
#>    1:     3930
#>    2:     3978
#>    3:     3807
#>    4:     3729
#>    5:     3891
#>   ---         
#>  996:     3801
#>  997:     4016
#>  998:     3851
#>  999:     4043
#> 1000:     3892
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

      expression      min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
      <bch:expr> <bch:tm> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
    1 fsummary      1.68s  1.73s     0.549     1.2GB     6.10    10   111     18.21s
    2 posterior     6.03s  6.27s     0.158    5.13GB     4.43    10   281      1.06m

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
    1 fsummary   190.43ms 198.94ms     4.46     92.1MB     2.67    10     6      2.24s
    2 posterior     1.09s    1.12s     0.822     759MB     4.85    10    59     12.17s

![](man/figures/comparison.png)

# TODO

- ✔ ~~parallelization with `mirai`~~
- faster convergence metrics with better ✔ ~~ranking~~ / qnorm (I’ve
  gotten the inverse normal transformation going twice as fast in Julia
  at least…)
- Option for FFT autocovariance if user knows they have badly converged
  parameters
