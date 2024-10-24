// rt-estim-gamma with Integrated Brownian Motion prior
functions{
 vector my_reverse(vector v) {
   int num = rows(v);
   vector[num] v2;
   for (i in 1:num) {
     v2[i] = v[num - i + 1];
   }
   return(v2);
 }

}
data {
  // Data length (number of time steps)
  int<lower=0> n;
  // Number of discretized gamma values (previous incidence days to use)
  int<lower=0> d; 
  // Observed cases at time t
  array[n] int<lower=0> obs;
  // Tests performed at time t
  array[n] int<lower=0> test;
  // Discretized gamma values
  vector[d] w; 
  
  vector[d+1] delay_weights;
  int<lower=0> prev_vals;

  real log_incid_rate_mean;
  real log_incid_rate_sd;
    
  real log_sigma_mu;
  real log_sigma_sd;
  
  real log_rho_mu;
  real log_rho_sd;
  
  real log_r0_mu;
  real log_r0_sd;
  
  real kappa_mu;
  real kappa_sd;
  
}

// The parameters accepted by the model.
parameters {
  matrix[n-1,2] Z; // std normals to scale to BM and IBM
  real log_rt0_raw; // initial log r0
  real log_incid_rate_raw; // standard deviation for incidence
  array[prev_vals] real<lower=0> seed_incid_raw;
  real <lower=0> log_rho; // raw underreporting parameter for cases as a function of tests
  array[n] real<lower=0> i_raw; //incidence transformed so that the mean is always 1
  real<lower = 0> kappa;  // overdispersion parameter for cases
  real<lower=0> exp_rate;
  real log_sigma; 
}

transformed parameters{
  matrix[n, 2] log_rt; // first column BM: derivative of log_Rt and second column IBM: log_rt
  real<lower=0> sigma = exp(log_sigma_mu + log_sigma_sd * log_sigma);
  real<lower=0> incid_rate = exp(log_incid_rate_raw * log_incid_rate_sd + 
                                 log_incid_rate_mean);
  vector<lower=0>[prev_vals] seed_incid;
  array[n] real<lower=0> incid; //true unobserved incidence, transformed from raw
  array[n] real<lower=0> weighted_sum;
  array[n] real <lower=0> delay_sum;
  vector<lower=0>[prev_vals+n] aug_incid;
  real <lower=0> rho = exp(log_rho * log_rho_sd + log_rho_mu);

  // seed incid
  for (i in 1:prev_vals){
    seed_incid[i] = seed_incid_raw[i]*inv(exp_rate);
  }
  
  // IBM prior on log_rt, transform Z which is a std normal random vector
  log_rt[1,1] = 0 + log_r0_sd * log_rt0_raw;
  log_rt[1,2] = log_r0_mu + log_r0_sd * log_rt0_raw;
  // Cholesky decomp of [[sigma, (sigma)^2/2], [(sigma)^2/2, (sigma)^3/3]] is [[sqrt(sigma), 0], [sigma^(3/2)/2, sigma^(3/2)/(sqrt(12))]]
  for (r in 2:n){
    // log_rt[r] represents the r-th row of the 2d markov process BM/IBM
    // scale standard normals sampled in model block to bivariate normal w/ BM/IBM conditional mean & covariance matrix
    // do Y = Sigma^1/2 * Z + mu
    log_rt[r] = ([[sqrt(sigma), sigma^(1.5)/2], [sigma^(1.5)/2, sigma^(1.5)/(sqrt(12))]]*(Z[r-1])')' + [ log_rt[r-1,1], log_rt[r-1,2] + sigma*log_rt[r-1,1]];
  }
    
  // Latent incidence
  for (a in 1:prev_vals) {
      aug_incid[a] = seed_incid[a];
    }  
  { int c;
    c = prev_vals;
    for (t in 1:n) {
      weighted_sum[t] = 0;
      if (c <= d ) {
        weighted_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(w, c)));

        }
      if (c > d) {
        weighted_sum[t] = dot_product(segment(aug_incid, c - d + 1, d), 
                                      my_reverse(head(w, d)));
          
    
        }
      incid[t] = exp(log_rt[t,2]) * weighted_sum[t] * i_raw[t];
      aug_incid[t + prev_vals] = incid[t];
      c = c + 1;
      }
    }
      { int c;
    // changing the indexing on this so that it allows for weight on the same day as 
    // the observation
    c = prev_vals +1;
    for (t in 1:n) {
      delay_sum[t] = 0;
      if (c <= d ) {
        delay_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(delay_weights, c)));
        }
      if (c > d) {
        delay_sum[t] = dot_product(segment(aug_incid, c - d + 1, d), 
                                      my_reverse(head(delay_weights, d)));
        }
      c = c + 1;
      }
    }
}
// The model to be estimated
model {
  //prior on underreporting
  log_rho ~ std_normal();

  //prior on transformed kappa
  kappa ~ normal(kappa_mu,kappa_sd) T[0,];

  // prior on random walk sigma
  log_sigma ~ std_normal();

  //priors for log R0
  log_rt0_raw ~ std_normal();

  //prior on sd of incidence
  log_incid_rate_raw ~ std_normal();
  
    //seed incidence
  exp_rate ~ exponential(0.3);
  
  for (i in 1:prev_vals) {
      seed_incid_raw[i] ~ exponential(1);
  }
  
  for (i in 1:(n-1)){
    Z[i,1] ~ std_normal();
    Z[i,2] ~ std_normal();
  }
  // log_rt[1,1] ~ normal(0, log_r0_sd);
  // log_rt[1,2] ~ normal(log_r0_mu, log_r0_sd);
  // for (r in 2:n){
  //   log_rt[r] ~ multi_normal( [log_rt[r-1, 1], log_rt[r-1, 2] + sigma*log_rt[r-1,1]], [[sigma, sigma^2/2],[sigma^2/2, sigma^3/3]]);
  // }

  //raw incidence modeled in terms of previous true incidence and seed incidence
  for (t in 1:n) {
    i_raw[t] ~ gamma(exp(log_rt[t,2]) * weighted_sum[t] * incid_rate, 
                exp(log_rt[t,2])* weighted_sum[t] * incid_rate);

  }

  //observed incidence model
      for (t in 1:n) {
    obs[t] ~ neg_binomial_2(rho*test[t] * delay_sum[t], kappa);
  }

}

//generated quantities{
//  array[n] int gen_obs;
//  for (t in 1:n) {
//     gen_obs[t] = neg_binomial_2_rng(rho * test[t] * delay_sum[t], kappa);
//  }
//}

generated quantities {
  array[n] int gen_obs;
  for (t in 1:n) {
    real mean_param = rho * test[t] * delay_sum[t];
    real dispersion_param = kappa;

    // Clamp the mean parameter to a maximum value
    if (mean_param > 1e5) {
      mean_param = 1e5;
    }

    // Ensure the dispersion parameter is within a valid range
    if (dispersion_param > 1e5) {
      dispersion_param = 1e5;
    }

    gen_obs[t] = neg_binomial_2_rng(mean_param, dispersion_param);
  }
}

