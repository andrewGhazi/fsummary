# library(furrr)
#
# plan(multisession, workers = 4)

frhat2 = function(z_scaled, z_scaled_folded, n_iter) {

  half_iter = floor(n_iter/2)

  variables = names(z_scaled) |>
    grep(pattern = "\\.chain|\\.iteration|\\.draw", value = TRUE, invert = TRUE)

  frhat_bulk = frhat_grp_df(z_scaled, half_iter, variables)
  frhat_tail = frhat_grp_df(z_scaled_folded, half_iter, variables)

  pmax(frhat_bulk, frhat_tail) |>
    qDT("variable") |>
    setNames(c("variable", "rhat"))
}

fess2 = function(ddff, n_iter, n_chain) {

  variables = names(ddff) |>
    grep(pattern = "\\.chain|\\.iteration|\\.draw", value = TRUE, invert = TRUE)

  # posterior::autocovariance is actually amazingly fast

  acov = ddff |> # This is by far the slowest part, 76% of time spent up through the first while loop is here
    gby(`.chain`) |>
    mtt(across(variables,
               posterior::autocovariance))

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

  data.table(variable = variables,
             ess_bulk = ess / tau_hat)
}

fess_tail2 = function(ddff, q_df, half_iter, two_chain) {

  variables =  names(ddff) |>
    grep(pattern = "\\.chain|\\.iteration|\\.draw", value = TRUE, invert = TRUE)

  q_df = q_df |>
    sbt(variable %in% variables) |>
    mtt(variable = factor(variable,
                          levels = variables)) |>
    roworder(variable)

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

  ess_tail5  = fess2( q5_I, half_iter, two_chain)$ess_bulk
  ess_tail95 = fess2(q95_I, half_iter, two_chain)$ess_bulk

  data.table(variable = variables,
             ess_tail = pmin(ess_tail5, ess_tail95))
}

fsummary_par = function(ddf, conv_metrics = TRUE, .cores = 4) {

  k = 1
  my_opts = furrr_options(globals = FALSE,
                          packages = c( "collapse"),
                          seed = 123,
                          scheduling = FALSE)

  ddf = ddf |> qDT()

  variables = names(ddf) |>
    grep(pattern = "\\.chain|\\.iteration|\\.draw", value = TRUE, invert = TRUE)

  n_iter = max(ddf$`.iteration`)
  n_chain = max(ddf$`.chain`)
  half_iter = floor(n_iter/2)
  two_chain = floor(n_chain*2)

  no_dots = ddf |> get_vars(variables)

  dot_i = names(ddf) |> startsWith(".") |> which()

  chunk_per_core = 1 # TODO play around with this.

  core_i = seq_along(variables) %% floor(.cores * chunk_per_core) + 1

  df_list = lapply(1:floor(.cores * chunk_per_core),
                   \(i) ddf |> get_vars(c(whichv(core_i, i),
                                          dot_i)))

  res = data.table(variable =        variables,
                   mean     =   fmean(no_dots, nthreads = .cores),
                   median   = fmedian(no_dots, nthreads = .cores))

  # sd parallel ----
  # sd_df = df_list |>
  #   parallel::mclapply(\(x) x |> slt(-`.iteration`, -`.draw`, -`.chain`) |> fsd() |> qDT("variable"),
  #                      mc.cores = .cores) |>
  #   rbindlist() |>
  #   setnames(2, "sd")

  sd_df = df_list |>
    future_map(\(x) x |> slt(-`.iteration`, -`.draw`, -`.chain`) |> fsd() |> qDT("variable"),
               .options = my_opts) |>
    rbindlist() |>
    setnames(2, "sd")

  res = join(res, sd_df,
             on = "variable",
             validate = "1:1",
             verbose = FALSE)

  # mad parallel ----
  .get_mad = \(x) {
    nd = x |> slt(-`.iteration`, -`.draw`, -`.chain`)
    (1.4826 * fmedian(abs(TRA(nd, fmedian(nd))))) |>
      qDT('variable')
  }

  # mad_df = df_list |>
  #   parallel::mclapply(.get_mad,
  #                      mc.cores = .cores) |>
  #   rbindlist() |>
  #   setnames(2, "mad")
  mad_df = df_list |>
    future_map(.get_mad,
               .options = my_opts) |>
    rbindlist() |>
    setnames(2, "mad")

  res = join(res, mad_df,
             on = "variable",
             validate = "1:1",
             verbose = FALSE)

  # quantiles parallel ----
  .get_qs = \(x) {
    x |> slt(-`.iteration`, -`.draw`, -`.chain`) |>
      dapply(fquantile, probs = c(.05, .95)) |>
      mtt(quant_name = c("q5", "q95")) |>
      data.table::transpose(keep.names = "variable",
                            make.names = "quant_name")
  }

  # q_df = df_list |>
  #   parallel::mclapply(.get_qs,
  #                      mc.cores = .cores) |>
  #   rbindlist()
  q_df = df_list |>
    future_map(.get_qs,
               .options = my_opts) |>
    rbindlist()

  res = join(res, q_df,
             on = "variable",
             validate = "1:1",
             verbose = FALSE)

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

    #draw data frame with fold. Should really rename to this to split instead of fold...
    .split = \(x, half_iter, n_chain) {
      x |>
      mtt(split = `.iteration` > half_iter) |>
        mtt(`.chain` = data.table::fifelse(split,
                                           `.chain` + n_chain,
                                           `.chain`),
            `.iteration` = data.table::fifelse(split,
                                               `.iteration` - half_iter,
                                               `.iteration`)) |>
        slt(-split)

    }

    # ddfs_list = df_list |>
    #   parallel::mclapply(.split, mc.cores = .cores)
    ddfs_list = df_list |>
      future_map(.split,
                 .options = my_opts,
                 half_iter = half_iter,
                 n_chain = n_chain)

    .z_scale = \(x, rank_fun, n_iter, n_chain) {
      x |>
        mtt(across(!(names(x) %in% c(".chain", ".iteration", ".draw")),
                   \(x) qnorm((rank_fun(x) - 3/8) / (n_iter*n_chain - 2 * 3/8 + 1))))
    }

    .fold_z_scale = \(x, rank_fun, n_iter, n_chain) {
      x |>
        mtt(across(!(names(x) %in% c(".chain", ".iteration", ".draw")),
                   \(x) abs(x - fmedian(x)))) |> # This is posterior:::fold_draws()
        mtt(across(!(names(x) %in% c(".chain", ".iteration", ".draw")),
                   \(x) qnorm((rank_fun(x) - 3/8) / (n_iter*n_chain - 2 * 3/8 + 1))))
    }

    # z_scaled = ddfs_list |>
    #   parallel::mclapply(.z_scale,
    #                      mc.cores = .cores)
    z_scaled = ddfs_list |>
      future_map(.z_scale,
                 .options = my_opts,
                 rank_fun = rank_fun, n_iter = n_iter, n_chain = n_chain)

    # z_scaled_folded = ddfs_list |>
    #   parallel::mclapply(.fold_z_scale,
    #                      mc.cores = .cores)
    z_scaled_folded = ddfs_list |>
      future_map(.fold_z_scale,
                 .options = my_opts,
                 rank_fun = rank_fun, n_iter = n_iter, n_chain = n_chain)

    # rh_df = parallel::mcmapply(frhat2,
    #                    z_scaled, z_scaled_folded,
    #                    MoreArgs = list(n_iter = n_iter),
    #                    mc.cores = .cores, SIMPLIFY = FALSE) |>
    #   rbindlist()

    rh_df = future_map2(z_scaled, z_scaled_folded,
                        frhat2,
                        n_iter = n_iter,
                        .options = furrr_options(globals = "frhat_grp_df",
                                                 packages = c("collapse"),
                                                 seed = 123,
                                                 scheduling = FALSE)) |>
      rbindlist()

    bulk_df = z_scaled |>
      future_map(fess2,
                 n_iter = half_iter,
                 n_chain = two_chain,
                 .options = furrr_options(globals = "frhat_grp_df",
                                          packages = c("collapse", "data.table"),
                                          seed = 123,
                                          scheduling = FALSE)) |>
      rbindlist()

    tail_df = ddfs_list |>
      future_map(fess_tail2,
                 q_df = q_df,
                 half_iter = half_iter,
                 two_chain = two_chain,
                 .options = furrr_options(globals = "fess2",
                                          packages = c("collapse", "data.table"),
                                          seed = 123,
                                          scheduling = FALSE)) |>
      rbindlist()

    res = join(res, rh_df,
               on = "variable",
               validate = "1:1", verbose = FALSE) |>
      join(bulk_df,
           on = "variable",
           validate = "1:1",
           verbose = FALSE) |>
      join(tail_df,
           on = "variable",
           validate = "1:1",
           verbose = FALSE)
  }

  res[]
}

fsummary_par(ddf, .cores = 4)
