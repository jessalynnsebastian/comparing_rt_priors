# Fit models to example simulated data
# Source functions
setwd(here::here())
list.files("R/Functions", full.names = TRUE) |>
  purrr::walk(source)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load example sim data
#sim_data <- read.csv("Data/S3_weeklydata.csv")
#data_length <- nrow(sim_data)
# Now testing new sim data
sim_data <- readRDS("Data/test_sim.rds")
sim_data <- sim_data$observed
# only run on 26 weeks
sim_data <- sim_data[18:44, ]
data_length <- nrow(sim_data)

# Init function
init_func <- function() {
  list(mu0_raw = 1,
       log_incid_rate_raw = 0,
       log_rt0_raw = 0,
       rho = 9E-5,
       kappa = 5,
       seed_incid_one_raw = 1,
       i_raw = rep(1, length.out = data_length))
}

# GRW prior
grw_fit <- fit_estimgamma_model(sim_data,
                                gen_params = c(7 / 4, 7 / 7.5),
                                delay_params = c(1, 7 / 4),
                                prev_vals = 8,
                                log_nu_mean = -2,
                                log_nu_sd = 0.7,
                                log_sigma_mean = -0.6,
                                log_sigma_sd = 0.6,
                                log_rho_mean = log(0.016 / 1023),
                                log_rho_sd = 0.3,
                                log_r0_mean = log(1),
                                log_r0_sd = 0.75,
                                kappa_mean = 70,
                                kappa_sd = 80,
                                init_func = init_func,
                                iterations = 4000,
                                thin = 2,
                                adapt_delta = 0.99,
                                treedepth = 12,
                                seed = 45,
                                chain = 4,
                                gen_dist = "hypo-exp",
                                delay_dist = "gamma",
                                prior = "random walk")
saveRDS(grw_fit, file = "Results/SimData/grw_fit_new_sim.rds")

# IBM prior
ibm_fit <- fit_estimgamma_model(sim_data,
                                gen_params = c(7 / 4, 7 / 7.5),
                                delay_params = c(1, 7 / 4),
                                prev_vals = 8,
                                log_nu_mean = -2,
                                log_nu_sd = 0.7,
                                log_sigma_mean = -2,
                                log_sigma_sd = 0.6,
                                log_rho_mean = log(0.016 / 1023),
                                log_rho_sd = 0.3,
                                log_r0_mean = log(1),
                                log_r0_sd = 0.75,
                                kappa_mean = 70,
                                kappa_sd = 80,
                                init_func = init_func,
                                iterations = 4000,
                                thin = 2,
                                adapt_delta = 0.99,
                                treedepth = 12,
                                seed = 45,
                                chain = 4,
                                gen_dist = "hypo-exp",
                                delay_dist = "gamma",
                                prior = "integrated brownian motion")
saveRDS(ibm_fit, file = "Results/SimData/ibm_fit_new_sim.rds")

# OU prior
ou_fit  <- fit_estimgamma_model(sim_data,
                                gen_params = c(7 / 4, 7 / 7.5),
                                delay_params = c(1, 7 / 4),
                                prev_vals = 8,
                                log_nu_mean = -2,
                                log_nu_sd = 0.7,
                                log_delta_mean = -0.6,
                                log_delta_sd = 0.6,
                                theta_mean = 0.2,
                                log_rho_mean = log(0.016 / 1023),
                                log_rho_sd = 0.3,
                                log_r0_mean = log(1),
                                log_r0_sd = 0.75,
                                kappa_mean = 70,
                                kappa_sd = 80,
                                init_func = init_func,
                                iterations = 4000,
                                thin = 2,
                                adapt_delta = 0.99,
                                treedepth = 12,
                                seed = 45,
                                chain = 4,
                                gen_dist = "hypo-exp",
                                delay_dist = "gamma",
                                prior = "ornstein-uhlenbeck")
saveRDS(ou_fit, file = "Results/SimData/ou_fit_new_sim.rds")

# Hilbert Space Gaussian Process Approx, Matern kernel
gp_fit  <- fit_estimgamma_model(sim_data,
                                gen_params = c(7 / 4, 7 / 7.5),
                                delay_params = c(1, 7 / 4),
                                prev_vals = 8,
                                log_nu_mean = -2,
                                log_nu_sd = 0.7,
                                log_rho_mean = log(0.016 / 1023),
                                log_rho_sd = 0.3,
                                log_r0_mean = log(1),
                                log_r0_sd = 0.75,
                                kappa_mean = 70,
                                kappa_sd = 80,
                                init_func = init_func,
                                iterations = 4000,
                                thin = 2,
                                adapt_delta = 0.99,
                                treedepth = 12,
                                seed = 45,
                                chain = 4,
                                gen_dist = "hypo-exp",
                                delay_dist = "gamma",
                                M = 10,
                                param_v = 3,
                                prior = "hsgp")
saveRDS(gp_fit, file = "Results/SimData/gp _fit_new_sim.rds")
