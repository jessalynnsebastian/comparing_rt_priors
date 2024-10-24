source("R/Functions/summarize_Rt.R")

# Read in sim data
sim_data <- readRDS("Data/test_sim.rds")$truth
sim_data <- sim_data[18:44, ]

grw_filepath <- "Results/SimData/grw_fit_new_sim.rds"
ou_filepath <- "Results/SimData/ou_fit_new_sim.rds"
ibm_filepath <- "Results/SimData/ibm_fit_new_sim.rds"
hsgp_filepath <- "Results/SimData/gp _fit_new_sim.rds"

## Gaussian Random Walk
simfit_rw <- readRDS(grw_filepath)
simfit_rw <- rstan::extract(simfit_rw)
rt_samples_rw <- exp(simfit_rw$log_rt)
summary_rw <- summarize_Rt(rt_samples_rw, true_rt = sim_data$Rt,
                           plot = TRUE,
                           plot_title = "Gaussian Random Walk Posterior")

## Ornstein-Uhlenbeck
# As of right now, OU has one chain with issues bad enough
# that it doesn't contain samples. Note 3 chains instead
# of 4.
simfit_ou <- readRDS(ou_filepath)
simfit_ou <- rstan::extract(simfit_ou)
rt_samples_ou <- exp(simfit_ou$log_rt)
summary_ou <- summarize_Rt(rt_samples_ou, true_rt = sim_data$Rt,
                           plot = TRUE,
                           plot_title = "Ornstein-Uhlenbeck Posterior")

## Integrated Brownian Motion
# As of right now, IBM has one chain that doesn't converge,
# so I remove the samples from that chain.
simfit_ibm <- readRDS(ibm_filepath)
simfit_ibm <- rstan::extract(simfit_ibm, permuted = FALSE)
simfit_ibm <- rbind(simfit_ibm[, 1, ], simfit_ibm[, 2, ], simfit_ibm[, 3, ])
rt_samples_ibm <- exp(simfit_ibm[, grepl("log_rt\\[\\d+,2\\]",
                                         colnames(simfit_ibm))])
summary_ibm <- summarize_Rt(rt_samples_ibm, true_rt = sim_data$Rt,
                            plot = TRUE,
                            plot_title = "Integrated Brownian Motion Posterior")

## Approximate GP with Matern 3/2 Kernel
simfit_hsgp <- readRDS(hsgp_filepath)
simfit_hsgp <- rstan::extract(simfit_hsgp)
rt_samples_hsgp <- exp(simfit_hsgp$log_rt)
summary_hsgp <- summarize_Rt(rt_samples_hsgp, true_rt = sim_data$Rt,
                             plot = TRUE,
                             plot_title = "Approx. GP (Matern 3/2) Posterior")

## Combined numeric summary
summary_df <- rbind(summary_rw,
                    summary_ou,
                    summary_ibm,
                    summary_hsgp)
View(summary_df)
