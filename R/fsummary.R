frhat_grp_df = function(gdf, n_iter, variables) {

  grp_df = gdf |> gby(`.chain`)

  chain_mean = grp_df |> smr(across(variables, fmean), keep.group_vars = FALSE)
  chain_var =  grp_df |> smr(across(variables, fvar), keep.group_vars = FALSE)

  var_bw = n_iter * fvar(chain_mean)
  var_wi = fmean(chain_var)

  sqrt((var_bw / var_wi + n_iter - 1) / n_iter)
}


frhat = function(z_scaled, z_scaled_folded, n_iter, n_chain, variables) {
  # too_clever() only works when there's no ties for continuous parameter values
  # sampled by stan. It breaks for parameters with small ranges that exhibit
  # recurring values.

  half_iter = floor(n_iter/2)

  frhat_bulk = frhat_grp_df(z_scaled, half_iter, variables)
  frhat_tail = frhat_grp_df(z_scaled_folded, half_iter, variables)

  pmax(frhat_bulk, frhat_tail) |>
    qDT("variable") |>
    setNames(c("variable", "rhat"))
}

fess = function(ddff, n_iter, n_chain, variables) {

  # posterior::autocovariance is actually amazingly fast

  # This is by far the slowest part, 76% of time spent up through the first while loop is here. 2.23s for dyingforacup
  # acov = ddff |>
  #   gby(`.chain`) |>
  #   mtt(across(variables,
  #              \(x) facov(x - fmean(x), # This is somehow faster than fwithin()
  #                         fvar(x))))

  # The Rcpp version is about 80% faster. Pretty good. I think it could be even faster if
  # the C++ function took the whole array at once, along with the indices corresponding to
  # chains, but then it would be harder to use the simple fft() function from armadillo.
  # Plus I think it would be harder to parallelize.

  .pad_X = function(X, k) {
    if (nrow(X) != 2*k) {
      X = rbind(X, matrix(0, nrow = 2*k - nrow(X), ncol = ncol(X)))
    }

    return(X)
  }

  .check_for_nans = function(res, vx, call = rlang::env_parent()) {
    # fftm() returns NaN for parameters with 0 variance. This happens sometimes during
    # ess_tail() calculations with short chains where an entire chain can sometimes fail
    # to have even a single entry that falls below the 5th percentile (or above the 95th).

    nans = is.nan(res[1,])

    if (any(nans)) {
      if (any(nans & (vx != 0))) {
        cli::cli_abort("Got NaN autocovariance on a chain with non-zero variance. WTF?! File a GitHub issue: {.url https://github.com/andrewGhazi/fsummary/issues}", call = call)
      } else {
        res[,vx == 0] = 0
      }
    }

    res
  }

  .fftm_chain = function(chain_dt) {
    nr = nrow(chain_dt)

    k = nextn(nr)

    vx = fvar(chain_dt)

    X = qM(chain_dt) |> .pad_X(k)

    fX = fftm(X)*(2*k)

    res = TRA(fX[1:nr,,drop=FALSE], vx / fX[1,] * (nr-1)/nr, FUN = "*") |>
      .check_for_nans(vx)

    res |> qDT()
  }

  # RcppArmadillo version about 1.75x faster
  acov = ddff |>
    gby(`.chain`) |>
    mtt(across(variables,
               fwithin)) |>
    fungroup() |>
    slt(".chain", variables) |>
    split(by = ".chain", keep.by = FALSE) |>
    lapply(.fftm_chain) |> # replace with future_map here
    rbindlist() |>
    setNames(variables) |>
    add_vars(ddff |> slt(".chain", ".iteration", ".draw"))

  acov_means = acov |>
    gby(`.iteration`) |>
    smr(across(variables,
               fmean))

  #
  mean_var_df = acov_means |>
    sbt(`.iteration` == 1) |>
    pivot(ids = ".iteration") |>
    mtt(mv = value * n_iter / (n_iter - 1),
        var_p = value) # it seems like posterior:::.ess has a redundant multiply

  if (n_chain > 1) {
    # should be possible to do this block and the command above without pivots,
    # which might be expensive for many variables. But they're small pivots (1xD), so
    # probably high-hanging fruit.
    chain_means = ddff |>
      gby(`.chain`) |>
      smr(across(variables,
                 fmean)) |>
      smr(across(variables,
                 fvar)) |>
      pivot()

    mean_var_df[,var_p := var_p +  chain_means$value]
  }

  # while loop 1 prep ----
  tacov_mean_mat = acov_means[,..variables] |> qM() |> t() |> unname()

  rh_m = matrix(0, nrow = length(variables), ncol = n_iter)
  t = 0
  rhe = rep(1, length(variables))
  rhe_final = rep(0, length(variables))
  rh_m[,1] = 1
  rho = 1 - (mean_var_df$mv - tacov_mean_mat[,2])/mean_var_df$var_p
  rh_m[,2] = rho
  epo = rhe+rho
  track = seq_along(variables)
  max_t = rep(2, length(variables))

  # while loop 1 ----
  while (t < ((n_iter) - 5) && any(!is.nan(epo[track])) && any(epo[track] > 0)) {
    max_t[track] = t
    t = t + 2

    rhe[track] = 1 - (mean_var_df$mv[track] - tacov_mean_mat[track,t+1]) / mean_var_df$var_p[track]
    rho[track] = 1 - (mean_var_df$mv[track] - tacov_mean_mat[track,t+2]) / mean_var_df$var_p[track]

    epo[track] = rhe[track] + rho[track]

    ei = epo[track] > 0
    drp = !ei
    rhe_final[track[drp]] = rhe[track[drp]] #  & rhe[track[drp]] > 0

    # if (t == 30) browser()

    if (any(ei)) {
      rh_m[track[ei],t + 1] = rhe[track][ei]
      rh_m[track[ei],t + 2] = rho[track][ei]
    }

    track = track[ei]
  }

  max_t = max_t + 2
  rhe_pos = rhe_final > 0
  to_add = which(rhe_pos)
  rh_m[cbind(to_add, max_t[rhe_pos]+1)] = rhe_final[rhe_pos]

  mt_df = data.table(i = seq_along(variables),
                     mt = max_t)

  # while loop 2 ----

  max_max = max(max_t)

  t = 0

  while (t <= max_max - 4) {
    t = t+2
    relevant = mt_df |> sbt((mt - 2) >= t)
    var_i = relevant$i

    two_ahd =  t +  1:2
    two_bed =  t + -1:0

    two_ahd_sum = rh_m[var_i, two_ahd, drop = FALSE] |> matrixStats::rowSums2()
    two_bhd_sum = rh_m[var_i, two_bed, drop = FALSE] |> matrixStats::rowSums2()

    update_i = which(two_ahd_sum > two_bhd_sum)

    if (!rlang::is_empty(update_i)) {
      res = two_bhd_sum[update_i]/2
      rh_m[var_i[update_i], t + 1] = res
      rh_m[var_i[update_i], t + 2] = res
    }
  }

  ess = n_chain * n_iter
  tau_bound = 1/log10(ess)

  rh_m_t = rh_m[,seq_len(min(max_max+4, n_iter))] |> t()

  tau_hat = rep(0, length(max_t))

  for (i in seq_along(variables)) {
    imax = max_t[i]
    tau_hat[i] = -1 + 2 * fsum(rh_m_t[1:imax,i]) + rh_m_t[imax+1,i]
  }

  if (any(tau_hat < tau_bound)) cli::cli_warn("ESS capped to avoid unstable estimates for {variables[which(tau_hat < tau_bound)]}")

  ess / tau_hat
}

fess_tail = function(ddff, q_df, half_iter, two_chain, variables) {

  q5_I = ddff |>
    get_vars(variables) |>
    qM() |>
    TRA(q_df$q5, FUN = "-") |>
    magrittr::is_weakly_less_than(0) |>
    cbind(ddff |> get_vars(c(".chain", ".iteration", ".draw")))

  q95_I = ddff |>
    get_vars(variables) |>
    qM() |>
    TRA(q_df$q95, FUN = "-") |>
    magrittr::is_weakly_less_than(0) |>
    cbind(ddff |> get_vars(c(".chain", ".iteration", ".draw")))

  ess_tail5  = fess( q5_I, half_iter, two_chain, variables)
  ess_tail95 = fess(q95_I, half_iter, two_chain, variables)

  pmin(ess_tail5, ess_tail95)
}

#' Fast summary function for cmdstanr draws
#'
#' @param ddf a draws_df
#' @param conv_metrics logical indicating whether to compute convergence metrics (slower)
#' @param .cores cores to use when computing mean/median
#' @details Currently .cores is only used by fmean and fmedian, so it likely won't help
#' much unless you have loads of parameters and/or extremely long chains.
#' @returns a data.table of summary metrics
#' @export
fsummary = function(ddf,
                    conv_metrics = TRUE,
                    .cores = getOption("mc.cores", 1)) {

  if (!inherits(ddf, "draws_df")) cli::cli_abort("Input {.var ddf} must be a {.cls draws_df}")

  ddf = ddf |> qDT()

  variables = names(ddf) |>
    grep(pattern = "\\.chain|\\.iteration|\\.draw", value = TRUE, invert = TRUE)

  n_iter = max(ddf$`.iteration`)
  n_chain = max(ddf$`.chain`)
  half_iter = floor(n_iter/2)
  two_chain = floor(n_chain*2)

  no_dots = ddf |> get_vars(variables)

  q_df = no_dots |>
    dapply(fquantile, probs = c(.05, .95)) |>
    data.table::transpose() |>
    setNames(c("q5", "q95"))

  res = data.table(variable =        variables,
                   mean     =   fmean(no_dots, nthreads = .cores),
                   median   = fmedian(no_dots, nthreads = .cores),
                   sd       =     fsd(no_dots),
                   mad      = 1.4826 * fmedian(abs(TRA(no_dots, fmedian(no_dots))))) |>
    add_vars(q_df)

  # The hard/slow part: rhat & ess ----

  if (conv_metrics) {
    rank_fun = ifelse(n_chain*n_iter > 3e4,
                      data.table::frank, # more overhead, but faster for very large vectors
                      base::rank)

    # almost all time spent on the INT. A faster rank function would help a lot.
    # mnorm::qnormFast is a tiny bit better (~10% faster for the full thing) but adds a
    # dependency and some approximation error. Combining both ranking and qnorm in
    # a Cpp function would probably be substantially better. Especially if it took the
    # variables in grouped long format (would require a pivot). Above my pay grade.

    #draw data frame with fold
    ddff = ddf |>
      qDT() |>
      mtt(fold = `.iteration` > ceiling(n_iter/2)) |>
      mtt(`.chain` = data.table::fifelse(fold,
                                         `.chain` + n_chain,
                                         `.chain`),
          `.iteration` = data.table::fifelse(fold,
                                             `.iteration` - ceiling(n_iter/2),
                                             `.iteration`))

    fold_meds = ddff |>
      slt(variables) |>
      fmedian()

    ddff = ddff |>
      sbt(`.iteration` <= half_iter) |>
      get_vars(c(variables, ".chain", ".iteration", ".draw"))

    # nrow(ddff) != nrow(ddf) always. Odd number of iterations -> uneven folded chain
    # lengths -> posterior:::.split_chains drops some values.

    z_scaled = ddff |>
      mtt(across(variables,
                 \(x) qnorm((rank_fun(x) - 3/8) / (nrow(ddff) - 2 * 3/8 + 1))))

    ess_tail_df = data.table(variable = variables,
                             ess_tail = fess_tail(ddff, q_df, half_iter, two_chain, variables))

    settransformv(ddff, variables, TRA, STATS = fold_meds, apply = FALSE)
    settransformv(ddff, variables, abs, apply = FALSE)

    z_scaled_folded = ddff |>
      mtt(across(variables,
                 \(x) qnorm((rank_fun(x) - 3/8) / (nrow(ddff) - 2 * 3/8 + 1))))

    rh_df = frhat(z_scaled, z_scaled_folded, n_iter, n_chain, variables)

    ess_bulk_df = data.table(variable = variables,
                             ess_bulk = fess(z_scaled, half_iter, two_chain, variables))



    res = join(res, rh_df,
               on = "variable",
               validate = "1:1", verbose = FALSE) |>
      join(ess_bulk_df,
           on = "variable",
           validate = "1:1",
           verbose = FALSE) |>
      join(ess_tail_df,
           on = "variable",
           validate = "1:1",
           verbose = FALSE)
  }

  res[]
}

