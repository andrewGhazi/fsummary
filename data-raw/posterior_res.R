## code to prepare `posterior_res` dataset goes here

library(fastverse)

set.seed(123)

n_iter = 51 # 101 to make sure the test hits the padding fn and uneven split chains
n_chain = 4
n_param = 3

test_ddf = rnorm(n_iter*n_chain*n_param) |>
  matrix(ncol = n_param) |>
  qDT() |>
  mtt(`.chain` = rep(1:4, each = n_iter),
      `.iteration`= rep(1:n_iter, times = n_chain),
      `.draw` = 1:(n_iter*n_chain)) |>
  posterior::as_draws_df()

test_res = test_ddf |>
  posterior::summarise_draws()

usethis::use_data(test_ddf, overwrite = TRUE)
# usethis::use_data(test_res, overwrite = TRUE, internal = TRUE)

# Large example - don't include the (big) ddf in the package for this one:

n = 4000
d = 1000

set.seed(123)

datM = runif(n*d) |>
  matrix(ncol = d)

if ((datM[1,1] - 0.2875775201246142) > 1e-8) cli::cli_abort("Something with the RNG went haywire, you didn't get the right simulated draws")

big_ddf = datM |>
  qDT() |>
  mtt(`.chain` = rep(1:4, each = 1000),
      `.iteration` = rep(1:1000, times = 4),
      `.draw` = 1:4000) |>
  posterior::as_draws_df()

big_test_res = big_ddf |>
  posterior::summarise_draws()

# usethis::use_data(big_test_res, overwrite = TRUE, internal = TRUE)

#### Generate example with NAs / 0 variance ----

library(cmdstanr)
tf = write_stan_file("
data {
  int<lower=1> N;
  array[N] real x;
  vector[N] y;
}
transformed data {
  vector[N] mu = rep_vector(0, N);
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[N, N] L_K;

  matrix[N, N] K = add_diag(gp_exp_quad_cov(x, alpha, rho), sigma^2);

  L_K = cholesky_decompose(K);
}
model {

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();

  y ~ multi_normal_cholesky(mu, L_K);
}
")

set.seed(123)
m = cmdstan_model(tf)
l = list(N = 5, x = rnorm(5), y = rnorm(5))
na_ddf = m$sample(l, refresh = 0,
                  show_messages = FALSE, show_exceptions = FALSE,
                  iter_sampling = 200,
                  seed = 123) |>
  posterior::as_draws_df()

na_res = na_ddf |>
  posterior::summarise_draws()

## Badly converged example ----

n_chain = 4
n_iter = 1000
n_param = 5

set.seed(123)

bad_conv_ddf = matrix(rnorm(n_chain*n_iter*n_param), ncol = n_param) |>
  qDT() |>
  mtt(`.draw` = 1:(n_chain*n_iter),
      `.iteration` = rep(1:n_iter, times = n_chain),
      `.chain` = rep(1:n_chain, each = n_iter)) |>
  mtt(V1 = V1 + 1*c(`.chain` == 2),
      V2 = V2 + 1*(`.chain` %in% c(3,4))) |> # shift the mean up by 1 on one or two chains
  posterior::as_draws_df()

bad_conv_post = bad_conv_ddf |>
  posterior::summarise_draws()

usethis::use_data(test_ddf,
                  test_res,
                  bad_conv_ddf,
                  bad_conv_post,
                  big_test_res,
                  na_ddf,
                  na_res,
                  overwrite = TRUE, internal = TRUE)
# usethis::use_data(na_res, overwrite = TRUE, internal = TRUE)

