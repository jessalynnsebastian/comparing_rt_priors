# Using a stochastic SEIRS model with some
# time-varying parameters, simulate an outbreak,
# and then layer in the testing:

sim_seir <- function(pop_size = 300000, init_exposed = 0, init_infectious = 10,
                     beta_tvar = function(t) {2},
                     sigma_tvar = function(t) {7 / 4},
                     gamma_tvar = function(t) {7 / 7.5},
                     omega_tvar = function(t) {1 / 16},
                     num_tests_tvar = function(t) {5000},
                     rho_tvar = function(t) {1e-5},
                     kappa = 5,
                     max_length_outbreak = 52) {
  # Initialize states
  curr_time <- 0
  curr_S <- pop_size - init_exposed - init_infectious
  curr_E <- init_exposed
  curr_I <- init_infectious
  curr_R <- 0
  # Track cumulative S2E and E2I
  S2E <- 0
  E2I <- 0
  outbreak_truth <- data.frame(time = 0,
                               S = curr_S,
                               E = curr_E,
                               I = curr_I,
                               R = curr_R,
                               Rt = (beta_tvar(0) / gamma_tvar(0)) * (curr_S / pop_size),
                               S2E = 0,
                               E2I = 0)
  # Assume there are always some infected or infectious people,
  # otherwise end the outbreak.
  # While the outbreak is still ongoing, do:
  while (curr_E > 0 || curr_I > 0) {
    # Get values of parameters at current times
    beta <- beta_tvar(curr_time)
    sigma <- sigma_tvar(curr_time)
    gamma <- gamma_tvar(curr_time)
    omega <- omega_tvar(curr_time)
    # Generate times to transitions
    nxt_times <- suppressWarnings(c(t_next_exposed = rexp(1, beta * curr_I * curr_S / pop_size),
                                    t_next_infectious = rexp(1, sigma * curr_E),
                                    t_next_removed = rexp(1, gamma * curr_I),
                                    t_next_imm_loss = rexp(1, omega * curr_R)))
    # NAs will be generated if curr_E or curr_I is 0, handle those by replacing them
    # with Inf.
    nxt_times[which(is.na(nxt_times))] <- Inf
    # Whichever time is shortest is the next event
    # First, record states if next event puts us into the next simulated timestep
    if (floor(curr_time + min(nxt_times)) != floor(curr_time)) {
      end_of_timestep <- c("time" = ceiling(curr_time),
                           "S" = curr_S,
                           "E" = curr_E,
                           "I" = curr_I,
                           "R" = curr_R,
                           "Rt" = (beta / gamma) * (curr_S / pop_size),
                           "E2I" = E2I,
                           "S2E" = S2E)
      outbreak_truth <- rbind(outbreak_truth, end_of_timestep)
    }
    # Then, jump to next states
    if (which.min(nxt_times) == 1) {
      # Next event is infection, "exposure"
      curr_time <- curr_time + nxt_times[1]
      curr_S <- curr_S - 1
      curr_E <- curr_E + 1
      S2E <- S2E + 1
    } else if (which.min(nxt_times) == 2) {
      # Next event is individual becoming infectious
      curr_time <- curr_time + nxt_times[2]
      curr_E <- curr_E - 1
      curr_I <- curr_I + 1
      E2I <- E2I + 1
    } else if (which.min(nxt_times) == 3) {
      # Next event is recovery
      curr_time <- curr_time + nxt_times[3]
      curr_I <- curr_I - 1
      curr_R <- curr_R + 1
    } else if (which.min(nxt_times) == 4) {
      # Next event is loss of immunity
      curr_time <- curr_time + nxt_times[4]
      curr_R <- curr_R - 1
      curr_S <- curr_S + 1
    }
    # Don't let the time go past the max desired length
    if (curr_time >= max_length_outbreak) break
  }

  # Now generate the case data from the E2I using negative binomial:
  # First get time-varying rho and tests for all times t (except 0)
  rho_t <- sapply(outbreak_truth$time[-1], rho_tvar)
  M_t <- sapply(outbreak_truth$time[-1], num_tests_tvar)
  E2I_diff <- diff(outbreak_truth$E2I)
  negbin_mean <- rho_t * M_t * E2I_diff
  negbin_variance <- negbin_mean + (negbin_mean^2 / kappa)
  size <- negbin_mean^2 / (negbin_variance - negbin_mean)
  prob <- negbin_mean / negbin_variance
  observed_cases <- rnbinom(length(negbin_mean), size = size, prob = prob)
  outbreak_observed <- data.frame(total_tests = M_t,
                                  total_cases = observed_cases)
  return(list(truth = outbreak_truth, observed = outbreak_observed))
}

# test different parameters
beta_func <- function(time) {
  # Want to produce a fairly "wiggly" curve
  beta <- 0.5 + 0.5 * exp(time / 150) * exp(sin(0.18 * (time + 3)))
  return(beta)
}
# trajectory of beta:
plot(sapply(1:52, beta_func))

# run the sim
test_sim <- sim_seir(beta_tvar = beta_func)
#test_sim <- readRDS("Data/test_sim.rds")

# Plot sim data
par(mfrow = c(3, 1))
# Plot S2E
plot(diff(test_sim$truth$S2E), type = "p", pch = 16, col = "black",
     xlab = "Week", ylab = "S2E",
     main = "Incidence", log = "y")
rect(xleft = 18, xright = 44,
     ybottom = 1, ytop = 1e6,
     col = rgb(0.6, 0.47, 1, alpha = 0.2))

# Plot observed total cases
plot(test_sim$observed$total_cases, type = "p", pch = 16, col = "blue",
     xlab = "Week", ylab = "Total Cases",
     main = "Observed Total Cases", log = "y")
rect(xleft = 18, xright = 44,
     ybottom = 0.01, ytop = 10000,
     col = rgb(0.6, 0.47, 1, alpha = 0.2))

# Plot Rt
plot(test_sim$truth$Rt, type = "l", col = "darkgreen",
     xlab = "Week", ylab = "Rt",
     main = "Rt Over Time")
rect(xleft = 18, xright = 44,
     ybottom = 0, ytop = 3,
     col = rgb(0.6, 0.47, 1, alpha = 0.2))

plot(test_sim$truth$S)

#saveRDS(test_sim, "Data/test_sim2.rds")
