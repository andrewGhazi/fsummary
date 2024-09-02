multiple_daemons = function() {
  is.matrix(mirai::status()$daemons) && nrow(mirai::status()$daemons) > 1
}

frhat_grp_df = function(gdf, n_iter, variables) {

  grp_df = gdf |> gby(`.chain`)

  chain_mean = grp_df |> smr(across(variables, fmean), keep.group_vars = FALSE)
  chain_var =  grp_df |> smr(across(variables, fvar), keep.group_vars = FALSE)

  var_bw = n_iter * fvar(chain_mean)
  var_wi = fmean(chain_var)

  sqrt((var_bw / var_wi + n_iter - 1) / n_iter)
}

frhat = function(z_scaled, z_scaled_folded, n_iter, n_chain, variables) {

  half_iter = floor(n_iter/2)

  frhat_bulk = frhat_grp_df(z_scaled       , half_iter, variables)

  frhat_tail = frhat_grp_df(z_scaled_folded, half_iter, variables)

  pmax(frhat_bulk, frhat_tail) |>
    qDT("variable") |>
    setNames(c("variable", "rhat"))
}


.ch1_chain_dt = function(chain_M) {
  fsummary:::cov_head(chain_M, n=1, offset=0)[1,]
}

get_ch1 = function(split_chains) {
  lapply(split_chains, .ch1_chain_dt)
}

get_chain_info = function(ddff, n_cov, offset) {
  ddff |>
    slt(".chain", ".iteration", ".draw") |>
    sbt(`.iteration` <= (n_cov + offset) &
          `.iteration` >= (offset + 1))
}

.check_for_nans = function(res, vx, call = rlang::env_parent()) {
  # fftm() returns NaN for parameters with 0 variance. This happens sometimes during
  # ess_tail() calculations with short chains where an entire chain can sometimes fail
  # to have even a single entry that falls below the 5th percentile (or above the 95th).
  # TODO: check how cov_head handles 0 variance chains.

  nans = is.nan(res[1,])

  if (any(nans)) {
    if (any(nans & (vx != 0))) {
      cli::cli_abort("Got NaN autocovariance on a chain with non-zero variance. WTF?! File a GitHub issue: {.url https://github.com/andrewGhazi/fsummary/issues}",
                     call = call)
    } else {
      res[,vx == 0] = 0
    }
  }

  res
}

.cov_head = function(chain_Ml, ch1, n_cov=n_cov, offset=offset) {

  M = chain_Ml

  nr = nrow(M)

  vx = fvar(M)

  chX = fsummary:::cov_head(M, n=n_cov, offset=offset)

  res = TRA(chX, vx / ch1 * (nr-1)/nr, FUN = "*") |>
    fsummary:::.check_for_nans(vx)

  res |> qDT()
}

get_acov_means = function(split_chains, ch1_by_chain, variables, n_cov, offset, chain_info) {

  n_iter = nrow(split_chains[[1]])

  minn = min(n_cov, n_iter)
  mino = min(offset, n_iter)

  acov = mapply(.cov_head,
                split_chains,
                ch1_by_chain,
                MoreArgs = list(n_cov = minn, offset = mino),
                SIMPLIFY = FALSE)

  acov |>
    rowbind(use.names = FALSE) |>
    setColnames(variables) |>
    add_vars(chain_info |> slt(".iteration")) |>
    gby(".iteration") |>
    fmean(use.g.names = FALSE)
}

fess = function(ddff, n_iter, n_chain, variables) {

  # This is by far the slowest part, 76% of time spent up through the first while loop is here. 2.23s for dyingforacup
  # acov = ddff |>
  #   gby(`.chain`) |>
  #   mtt(across(variables,
  #              \(x) facov(x - fmean(x), # This is somehow faster than fwithin()
  #                         fvar(x))))

  n_cov = 16
  offset = 0

  # RcppArmadillo version about 1.75x faster
  # acov = ddff |>
  #   gby(`.chain`) |>
  #   mtt(across(variables,
  #              fwithin)) |>
  #   fungroup() |>
  #   slt(".chain", variables) |>
  #   split(by = ".chain", keep.by = FALSE) |>
  #   lapply(.fftm_chain) |> # replace with future_map here
  #   rbindlist() |>
  #   setNames(variables) |>
  #   add_vars(ddff |> slt(".chain", ".iteration", ".draw") )

  split_chains = fsummary:::center_split_df(ddff |> slt(variables),
                                            ddff$.chain - 1,
                                            n_chain,
                                            n_iter)

  for (i in 1:n_chain) {
    split_chains[[i]] = setColnames(split_chains[[i]], variables)
  }

  ch1_by_chain = fsummary:::get_ch1(split_chains)

  if (length(variables) == 1) ch1_by_chain = lapply(ch1_by_chain, as.matrix)
  # ^ this gets re-used in the while loop(s) if additional covariance terms are needed

  chain_info = fsummary:::get_chain_info(ddff, n_cov, offset)

  acov_means = fsummary:::get_acov_means(split_chains, ch1_by_chain, variables, n_cov, offset, chain_info)

  # collapse::pivot() bugs out if it's only one variable
  mean_var_df = acov_means |>
    sbt(`.iteration` == 1) |>
    data.table::melt.data.table(id.vars = ".iteration") |> # pivot(ids = ".iteration") |>
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
      data.table::melt.data.table(measure.vars = variables) |>
      slt(cmv = "value")

    # Adding collapse:: namespaces here^ on across() totally fucks it up. This might
    # change but that leaves the options as:
    # 1. Rewrite the handful of commands that use across() to horrifying data.table commands or
    # 2. require everywhere(library(collapse)), which is not that bad I think.

    # chain_means = slt(ddff, ".chain", variables)[,lapply(.SD, mean),by=".chain",.SDcols = variables]
    # chain_mean_var = chain_means[,lapply(.SD, var),.SDcols = variables] |>
    #   collapse::pivot()
    mean_var_df = mean_var_df |>
      add_vars(chain_means)
    settransform(mean_var_df, var_p = var_p + cmv)

  }

  # while loop 1 prep ----
  tacov_mean_mat = acov_means |>
    slt(variables) |>
    qM() |>
    t() |>
    unname()

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

  n_cov_terms = n_cov # number of autocovariance terms we currently have
  offset = n_cov # where to start if we need to add more
  offset_inc = 8 # how many more we'll collect if needed

  # while loop 1 ----
  while (t < ((n_iter) - 5) && any(!is.nan(epo[track])) && any(epo[track] > 0)) {

    max_t[track] = t
    t = t + 2

    if (t > (n_cov_terms - 2)) {

      # cli::cli_warn("Entered addl acov loop with {length(track)} variable{?s}.")
      # oops, didn't collect enough acov terms, go get some more. This will be way slower
      # than the fft approach if you need to do it more than once for many variables.
      # TODO: test how often this happens with poorly mixed chains.

      n_cov_terms = n_cov_terms + offset_inc

      zm = matrix(0, nrow = length(variables), ncol = offset_inc)

      addl_chain_info = fsummary:::get_chain_info(ddff, offset_inc, offset)

      if (offset + offset_inc > n_iter) offset_inc = n_iter - offset

      addl_acov_means = fsummary:::get_acov_means(split_chains |> lapply(\(x) x[,variables[track],drop=FALSE]),
                                                  ch1_by_chain |> lapply(\(x) x[track]),
                                                  variables[track],
                                                  offset_inc, offset,
                                                  addl_chain_info)

      addl_tacov = addl_acov_means |>
        slt(-`.iteration`) |>
        qM() |>
        t() |>
        unname()

      zm[track,] = addl_tacov

      tacov_mean_mat = cbind(tacov_mean_mat, zm)
      offset = offset + offset_inc
    }

    rhe[track] = 1 - (mean_var_df$mv[track] - tacov_mean_mat[track,t+1,drop=FALSE]) / mean_var_df$var_p[track]
    rho[track] = 1 - (mean_var_df$mv[track] - tacov_mean_mat[track,t+2,drop=FALSE]) / mean_var_df$var_p[track]

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

  mt_df = data.table::data.table(i = seq_along(variables),
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

  rh_m_t = rh_m[,seq_len(min(max_max+4, n_iter)), drop = FALSE] |> t()

  tau_hat = rep(0, length(max_t))

  for (i in seq_along(variables)) {
    imax = max_t[i]
    tau_hat[i] = -1 + 2 * fsum(rh_m_t[1:imax,i]) + rh_m_t[imax+1,i]
  }

  if (any(tau_hat < tau_bound)) cli::cli_warn("ESS capped to avoid unstable estimates for {variables[which(tau_hat < tau_bound)]}")

  ess / tau_hat
}

get_q_df = function(ddf, variables, chunks_list) {

  res = ddf |>
    slt(variables) |>
    dapply(fquantile, probs = c(.05, .95)) |>
    data.table::transpose() |>
    setNames(c("q5", "q95"))

  return(res)
}

ltet = function(x,y) {
  x <= y
}

get_quantile_ind_df = function(ddff, q_df, variables, q) {
  ddff |>
    get_vars(variables) |>
    qM() |>
    TRA(get_elem(q_df, q), FUN = "-") |>
    fsummary:::ltet(0) |>
    qDT() |>
    add_vars(ddff |> get_vars(c(".chain", ".iteration", ".draw")))
}

fess_tail = function(ddff, q_df, half_iter, two_chain, variables) {
  q5_I = ddff |>
    get_vars(variables) |>
    qM() |>
    TRA(q_df$q5, FUN = "-") |>
    fsummary:::ltet(0) |>
    cbind(ddff |> get_vars(c(".chain", ".iteration", ".draw")))

  q95_I = ddff |>
    get_vars(variables) |>
    qM() |>
    TRA(q_df$q95, FUN = "-") |>
    fsummary:::ltet(0) |>
    cbind(ddff |> get_vars(c(".chain", ".iteration", ".draw")))

  ess_tail5  = fess( q5_I, half_iter, two_chain, variables)
  ess_tail95 = fess(q95_I, half_iter, two_chain, variables)

  data.table(variable = variables,
             ess_tail = pmin(ess_tail5, ess_tail95))
}

get_chunks = function(variables) {

  status_res = mirai::status()

  n_job = 2 * nrow(status_res$daemons)

  chunks = data.table::data.table(i = rep(seq_len(n_job),
                                          length.out = length(variables)),
                                  v = variables,
                                  vi = seq_along(variables)) |>
    roworder(i) |>
    split(by = "i")
}

fess_bulk = function(ddff, half_iter, two_chain, variables) {
  data.table(variable = variables,
             ess_bulk = fess(ddff, half_iter, two_chain, variables))
}

get_stats_df = function(ddf, variables) {
  no_dots = ddf |> get_vars(variables)

  res = data.table(variable =        variables,
                   mean     =   fmean(no_dots),
                   median   = fmedian(no_dots),
                   sd       =     fsd(no_dots),
                   mad      = 1.4826 * fmedian(abs(TRA(no_dots, fmedian(no_dots)))))

  return(res)
}

z_scale = function(x, n) {
  qnorm((fsummary:::myrank(x)[,1] - 3/8) / (n - 2 * 3/8 + 1))
}

z_scale_df = function(ddff, variables) {
  ddff |> # 3.7s
    mtt(across(variables,
               z_scale,
               n = nrow(ddff)))
}

get_folded_with_meds = function(ddf, variables, n_iter, n_chain, half_iter) {
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

  return(list(ddff, fold_meds))
}

.fsummary = function(ddf,
                     conv_metrics = TRUE,
                     verbose = TRUE,
                     chunks_list = NULL) {

  # if chunks are provided, recursively call this function on each chunk
  if (!is.null(chunks_list)) {
    # Man I hate this 2D mapping behavior of mirai_map

    # m_res = mirai::mirai_map(chunks_list,
    #                          .fsummary,
    #                          .args = list(chunks_list = NULL))[]

    .args = list(chunks_list = NULL,
                 conv_metrics = conv_metrics,
                 verbose = verbose)

    m_res <- vector(mode = "list", length = length(chunks_list))
    names(m_res) = names(chunks_list)

    for (i in seq_along(m_res)) {
      m_res[[i]] = mirai::mirai({.f(chunk_i,
                                    conv_metrics = cm,
                                    verbose = v,
                                    chunks_list = NULL)},
                                .args = list(.f = .fsummary,
                                             cm = conv_metrics,
                                             v = verbose,
                                             chunk_i = chunks_list[[i]]),
                                .compute = "default")
    }

    class(m_res) = "mirai_map"

    res = rowbind(m_res[])
  } else {

    strt = Sys.time()

    if (!inherits(ddf, "draws_df")) cli::cli_abort("Input {.var ddf} must be a {.cls draws_df}")

    variables = names(ddf) |> head(-3)
    ddf = ddf |> qDT()

    n_iter = max(ddf$`.iteration`)
    n_chain = max(ddf$`.chain`)

    half_iter = floor(n_iter/2)
    two_chain = floor(n_chain*2)

    strt = Sys.time()

    q_df = fsummary:::get_q_df(ddf, variables)

    if (verbose) cli::cli_alert("quantiles  {round(digits = 2, Sys.time() - strt)}")

    strt = Sys.time()

    stats_df = fsummary:::get_stats_df(ddf, variables)

    if (verbose) cli::cli_alert("stats {round(digits = 2, Sys.time() - strt)}")

    res = stats_df |>
      add_vars(q_df)

    # The hard/slow part: rhat & ess ----
    if (conv_metrics) {
      # almost all time spent on the INT. A faster rank function would help a lot.
      # mnorm::qnormFast is a tiny bit better (~10% faster for the full thing) but adds a
      # dependency and some approximation error. Combining both ranking and qnorm in
      # a Cpp function would probably be substantially better. Especially if it took the
      # variables in grouped long format (would require a pivot). Above my pay grade.

      strt = Sys.time()

      #draw data frame with fold
      folded_with_meds = fsummary:::get_folded_with_meds(ddf, variables, n_iter, n_chain, half_iter)
      ddff = folded_with_meds[[1]]
      fold_meds = folded_with_meds[[2]]

      if (verbose) cli::cli_alert("ddff setup {round(digits = 2, Sys.time() - strt)}")

      # nrow(ddff) != nrow(ddf) always. Odd number of iterations -> uneven folded chain
      # lengths -> posterior:::.split_chains drops some values.

      strt = Sys.time()

      z_scaled = fsummary:::z_scale_df(ddff, variables)

      if (verbose) cli::cli_alert("z_scale {round(digits = 2, Sys.time() - strt)}")

      strt = Sys.time()

      ess_tail_df = fsummary:::fess_tail(ddff, q_df, half_iter, two_chain, variables)

      if (verbose) cli::cli_alert("fess_tail {round(digits = 2, Sys.time() - strt)}")

      strt = Sys.time()

      demedian_abs = \(x, y=fold_meds) abs(TRA(x, STATS = y))

      settransformv(ddff, variables, demedian_abs, apply = FALSE)

      z_scaled_folded = fsummary:::z_scale_df(ddff, variables)

      if (verbose) cli::cli_alert("z_scale_fold {round(digits = 2, Sys.time() - strt)}")

      strt = Sys.time()

      rh_df = fsummary:::frhat(z_scaled, z_scaled_folded, n_iter, n_chain, variables)

      if (verbose) cli::cli_alert("rhat {round(digits = 2, Sys.time() - strt)}")

      strt = Sys.time()

      ess_bulk_df = fsummary:::fess_bulk(z_scaled, half_iter, two_chain, variables)

      if (verbose) cli::cli_alert("fess_bulk {round(digits = 2, Sys.time() - strt)}")

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

  }

  res[]
}

#' Fast summary function for cmdstanr draws
#'
#' @description
#' This function quickly computes the same thing as \link[posterior]{summarise_draws} with two minor differences:
#' \itemize{
#'   \item{it produces a \link[data.table]{data.table} instead of a \link[tibble]{tibble}}
#'   \item{there can be negligible floating point error way out past the decimal point}
#' }
#'
#' @param ddf a draws_df
#' @param conv_metrics logical indicating whether to compute convergence metrics (slower)
#' @details
#' \itemize{
#' \item{The input \code{ddf} can ONLY be a \link[posterior]{draws_df}. No arrays, lists, etc. Use \link[posterior]{as_draws_df} if needed.}
#' \item{Activate parallelization by calling \link[mirai]{daemons} with your desired
#'   number of cores before calling \code{fsummary}.}
#' }
#' @examples
#' fsummary(test_ddf)
#'
#' @returns a data.table of summary metrics
#' @export
fsummary = function(ddf,
                    conv_metrics = TRUE,
                    verbose = FALSE) {

  variables = names(ddf) |> head(-3)

  if (multiple_daemons()) {

    if (verbose) cli::cli_alert_info("Loading `collapse` on daemons with {.code mirai::everywhere(library(collapse))}")

    mirai::everywhere(library(collapse))

    chunks = get_chunks(variables)

    chunks_list = lapply(chunks, \(x) slt(ddf, c(x$v, ".chain", ".iteration", ".draw")));

    res = .fsummary(ddf = NULL,
                    conv_metrics = conv_metrics, verbose = verbose,
                    chunks_list = chunks_list) |>
      add_vars(rowbind(chunks)) |>
      roworder(vi) |>
      slt(-i, -v, -vi)
  } else {
    res = .fsummary(ddf,
                    conv_metrics = conv_metrics, verbose = verbose)
  }

  res
}

