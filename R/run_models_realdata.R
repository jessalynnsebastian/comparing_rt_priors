# Run models on real data using Feretti generation time
# Source functions
setwd(here::here())
list.files("R/Functions", full.names = TRUE) |>
  purrr::walk(source)

## CHANGE LATER:
indic <- 3 # this is for OC data only

# Load CA county data
ca_data <- read_csv("covid19cases_test.csv")
county_list <- c("Alameda", "Los Angeles", "Orange", "San Diego",
                 "San Francisco", "Santa Clara", "Riverside",
                 "San Bernardino", "Contra Costa", "Sacramento",
                 "Fresno", "Merced", "Monterey", "Stanislaus",
                 "Tulare")
county_name <- county_list[indic]
# Create weekly data, currently using data 08/02/2020 - 01/09/2022
county_data <- ca_data |>
  filter(area == county) |>
  mutate(week = epiweek(date),
         year = epiyear(date)) |>
  group_by(year, week) |>
  summarise(total_cases = sum(cases),
            total_tests = sum(total_tests),
            min_date = min(date)) |>
  ungroup() |>
  filter(min_date >= start_sunday & min_date <= end_sunday) |>
  mutate(time = row_number(),
         epidemia_time = time + 1)

# Run spline for kappa priors
county_spline <- brms::brm(brms::bf(total_cases ~ s(time)),
                           data = county_data, family = negbinomial(),
                           cores = 4, seed = 17, iter = 4000,
                           warmup = 1000, thin = 10, refresh = 0,
                           control = list(adapt_delta = 0.99))

# Calculate kappa priors
county_kappa <- choose_kappa_params(county_spline)

# Calculate quantile for tests
test_quantile <- quantile(county_data$total_tests)

# Fit the model

# Choose starting points
# Choose Rt starting point using EpiEstim
logrt_start <- get_logrtstart(county_data)

# Choose incidence starting point
incid_start <- 1 / 0.066 * county_data$total_cases

# Init function
init_func <- function() {
  list(log_incid_rate_raw = 0,
       log_rt0_raw = 0,
       rho = 0.066 / test_quantile[2],
       kappa = county_kappa$par[1],
       seed_incid_one_raw = 1,
       incid = incid_start,
       log_rt = logrt_start)
}

# Parameters:
# gen_params are the generation time parameters:
## first one is the latent period rate
## second is the infectious period rate
# delay_params are for a gamma dist'n, also used for latent.
# rho is the "case detection prior" rho * tests * cases
## is the mean of case observation model
# if you're fitting CA data, the rho prior should work well as is

# Gaussian Random Walk prior
grw_county_posterior <- fit_estimgamma_model(county_data,
                                             gen_params = c(3.2862,
                                                            6.1244 * 1 / 7),
                                             delay_params =
                                               c(4.05, 2 * 7 * 0.74),
                                             prev_vals = 4,
                                             log_rho_mean =
                                               log(0.066 /
                                                     test_quantile[2]),
                                             log_rho_sd = 0.3,
                                             log_r0_mean = log(1),
                                             log_r0_sd = 0.75,
                                             log_sigma_mean = -0.6,
                                             log_sigma_sd = 0.6,
                                             kappa_mean = county_kappa$par[1],
                                             kappa_sd = county_kappa$par[2],
                                             iterations = 6000,
                                             init_func = init_func,
                                             gen_dist = "weibull",
                                             seed = 225,
                                             thin = 3,
                                             prior = "random walk")
saveRDS(grw_county_posterior,
        file = "Results/RealData/grw_ca_fit.rds")

# Integrated Brownian Motion prior
ibm_county_posterior <- fit_estimgamma_model(county_data,
                                             gen_params = c(3.2862,
                                                            6.1244 * 1 / 7),
                                             delay_params =
                                               c(4.05, 2 * 7 * 0.74),
                                             prev_vals = 4,
                                             log_rho_mean =
                                               log(0.066 /
                                                     test_quantile[2]),
                                             log_rho_sd = 0.3,
                                             log_r0_mean = log(1),
                                             log_r0_sd = 0.75,
                                             log_sigma_mean = -2,
                                             log_sigma_sd = 0.6,
                                             kappa_mean = county_kappa$par[1],
                                             kappa_sd = county_kappa$par[2],
                                             iterations = 6000,
                                             init_func = init_func,
                                             gen_dist = "weibull",
                                             seed = 225,
                                             thin = 3,
                                             prior =
                                               "integrated brownian motion")
saveRDS(ibm_county_posterior,
        file = "Results/RealData/ibm_ca_fit.rds")

# Ornstein-Uhlenbeck prior
ou_county_posterior <- fit_estimgamma_model(county_data,
                                            gen_params = c(3.2862,
                                                           6.1244 * 1 / 7),
                                            delay_params =
                                              c(4.05, 2 * 7 * 0.74),
                                            prev_vals = 4,
                                            log_rho_mean =
                                              log(0.066 /
                                                    test_quantile[2]),
                                            log_rho_sd = 0.3,
                                            log_r0_mean = log(1),
                                            log_r0_sd = 0.75,
                                            log_delta_mean = -0.6,
                                            log_delta_sd = 0.6,
                                            theta_mu = 0.3,
                                            kappa_mean = county_kappa$par[1],
                                            kappa_sd = county_kappa$par[2],
                                            iterations = 6000,
                                            init_func = init_func,
                                            gen_dist = "weibull",
                                            seed = 225,
                                            thin = 3,
                                            prior = "ornstein-uhlenbeck")
saveRDS(ou_county_posterior,
        file = "Results/RealData/ou_ca_fit.rds")