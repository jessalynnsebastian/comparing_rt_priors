library(rstan)
source("R/Functions/discretize.R")
fit_estimgamma_model <- function(data,
                                 gen_params,
                                 delay_params,
                                 prev_vals,
                                 log_nu_mean = -2,
                                 log_nu_sd = 0.7,
                                 log_sigma_mean = -0.6,
                                 log_sigma_sd = 0.6,
                                 log_rho_mean,
                                 log_rho_sd,
                                 log_r0_mean = log(1),
                                 log_r0_sd = 0.75,
                                 kappa_mean,
                                 kappa_sd,
                                 init_func,
                                 iterations = 2000,
                                 thin = 2,
                                 adapt_delta = 0.99,
                                 treedepth = 12,
                                 seed = 45,
                                 chain = 4,
                                 gen_dist = "hypo-exp",
                                 delay_dist = "gamma",
                                 prior = "random walk",
                                 log_delta_mean = -0.6,
                                 log_delta_sd = 0.6,
                                 theta_mean = 0.3,
                                 M = 10,
                                 param_v = 3) {
  data_length <- dim(data)[1]

  if (gen_dist == "hypo-exp") {
    gen_weights <- epidemia_hypoexp(data_length, gen_params)
  }

  if (gen_dist == "log-normal") {
    gen_weights <- epidemia_lognormal(data_length, gen_params)
  }

  if (gen_dist == "weibull") {
    gen_weights <- epidemia_weibull(data_length, gen_params)
  }

  if (delay_dist == "gamma") {
    delay_weights <- zero_epidemia_gamma(data_length,
                                         delay_params[1],
                                         delay_params[2])
  }

  if (prior == "random walk" || prior == "integrated brownian motion") {
    model_object <- list(n = data_length,
                         d = data_length,
                         w = gen_weights,
                         delay_weights = delay_weights,
                         obs = data$total_cases,
                         test = data$total_tests,
                         prev_vals = 4,
                         log_incid_rate_mean = log_nu_mean,
                         log_incid_rate_sd = log_nu_sd,
                         log_sigma_mu = log_sigma_mean,
                         log_sigma_sd = log_sigma_sd,
                         log_rho_mu = log_rho_mean,
                         log_rho_sd = log_rho_sd,
                         log_r0_mu = log_r0_mean,
                         log_r0_sd = log_r0_sd,
                         kappa_mu = kappa_mean,
                         kappa_sd = kappa_sd)
    stanfile <- ifelse(prior == "random walk",
                       "Stan/rt_estim_gamma_RW.stan",
                       "Stan/rt_estim_gamma_IBM.stan")
  } else if (prior == "ornstein-uhlenbeck") {
    model_object <- list(n = data_length,
                         d = data_length,
                         w = gen_weights,
                         delay_weights = delay_weights,
                         obs = data$total_cases,
                         test = data$total_tests,
                         prev_vals = 4,
                         log_incid_rate_mean = log_nu_mean,
                         log_incid_rate_sd = log_nu_sd,
                         log_delta_mu = log_delta_mean,
                         log_delta_sd = log_delta_sd,
                         theta_mu = theta_mean,
                         log_rho_mu = log_rho_mean,
                         log_rho_sd = log_rho_sd,
                         log_r0_mu = log_r0_mean,
                         log_r0_sd = log_r0_sd,
                         kappa_mu = kappa_mean,
                         kappa_sd = kappa_sd)
    stanfile <- "Stan/rt_estim_gamma_OU.stan"
  } else if (prior == "hsgp") {
    model_object <- list(n = data_length,
                         d = data_length,
                         w = gen_weights,
                         delay_weights = delay_weights,
                         obs = data$total_cases,
                         test = data$total_tests,
                         prev_vals = 4,
                         log_incid_rate_mean = log_nu_mean,
                         log_incid_rate_sd = log_nu_sd,
                         log_rho_mu = log_rho_mean,
                         log_rho_sd = log_rho_sd,
                         log_r0_mu = log_r0_mean,
                         log_r0_sd = log_r0_sd,
                         kappa_mu = kappa_mean,
                         kappa_sd = kappa_sd,
                         M = M,
                         L = data_length + 1,
                         param_v = param_v)
    stanfile <- "Stan/rt_estim_gamma_HSGP.stan"
  }

  control_list <- list(adapt_delta = adapt_delta,
                       max_treedepth = treedepth)


  model_fit <- stan(file = stanfile,
                    data = model_object,
                    seed = seed,
                    iter = iterations,
                    thin = thin,
                    chain = chain,
                    init = init_func,
                    control = control_list)
  return(model_fit)
}