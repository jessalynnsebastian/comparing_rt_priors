choose_kappa_params <- function(spline_posterior) {
  posterior_pars <- summary(spline_posterior)
  start_mean <- posterior_pars[["spec_pars"]][[1]]
  start_sd <- posterior_pars[["spec_pars"]][[2]]
  true_lb <- posterior_pars[["spec_pars"]][[3]]
  true_ub <- posterior_pars[["spec_pars"]][[4]]
  start_params <- c(start_mean, start_sd)
  true_quantiles <- c(true_lb, true_ub)
  optim_params <- optim(par = start_params,
                        fn = compare_kappa_quantiles,
                        true_quantiles = true_quantiles)
  return(optim_params)
}

compare_kappa_quantiles <- function(candidate_params, true_quantiles) {
  candidate_quantiles <- quantile(truncnorm::rtruncnorm(10000,
                                                        a = 0,
                                                        mean = candidate_params[1],
                                                        sd = candidate_params[2]),
                                                        c(0.025, 0.975))
  loss <- (true_quantiles[1] - candidate_quantiles[1])^2 +
    (true_quantiles[2] - candidate_quantiles[2])^2
  return(loss)
}