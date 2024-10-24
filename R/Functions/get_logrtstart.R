# Use EpiEstim to initialize
get_logrtstart <- function(data,
                           window = 1,
                           gi_mean = 11.5 / 7) {
  window <- window
  gi_mean <- gi_mean
  gi_var <- 2 * (gi_mean / 2)^2
  ts <- data$time
  ts <- ts[ts > 1 & ts <= (max(ts) - window + 1)]
  te <- ts + (window - 1)
  ee_outs <- EpiEstim::estimate_R(
    incid = data$total_cases,
    method = "uncertain_si",
    config = EpiEstim::make_config(
      list(
        mean_si = gi_mean,
        min_mean_si = 1,
        max_mean_si = gi_mean + 1,
        std_mean_si = 1.5,
        std_std_si = 1.5,
        std_si = sqrt(gi_var),
        min_std_si = sqrt(gi_var) * .8,
        max_std_si = sqrt(gi_var) * 1.2,
        n1 = 50,
        n2 = 100,
        t_start = ts,
        t_end = te
      )
    )
  )
  ee_quantile <- ee_outs[["R"]] |>
    dplyr::select(t_start,
                  rt_mean = `Mean(R)`,
                  rt_median = `Median(R)`,
                  rt_CI95l = `Quantile.0.025(R)`,
                  rt_CI95u = `Quantile.0.975(R)`) |>
    dplyr::mutate(time  = t_start)
  log_ee_median <- log(ee_quantile |> dplyr::pull(rt_median))
  first_one <- log_ee_median[1]
  rt_start <- c(first_one, log_ee_median)
  return(rt_start)
}