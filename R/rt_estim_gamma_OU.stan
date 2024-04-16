//gamma model of incidence
//exponential seeded incidence
//observations as nb with mean as function of tests
//log scale Ornstein Uhlenbeck process as Rt prior
//exponential hyperprior on theta
//lognormal hyperprior on sigma
//normal hyperprior on mu
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
  //data length (number of time steps)
  int<lower=0> n;
  //number of discretized gamma values (previous incidence days to use)
  int<lower=0> d; 
  //observed cases at time t
  int<lower=0> obs[n];
  //tests performed at time t
  int<lower=0> test[n];
  //discretized gamma values
  vector[d] w; 
  
  vector[d+1] delay_weights;
  int<lower=0> prev_vals;

  real log_incid_rate_mean;
  real log_incid_rate_sd;
    
  real log_delta_mu;
  real log_delta_sd;
  
  real log_rho_mu;
  real log_rho_sd;
  
  real log_r0_mu;
  real log_r0_sd;
  
  real kappa_mu;
  real kappa_sd;
  
  real theta_mu;
  
}

// The parameters accepted by the model.
parameters {
  vector[n-1] z;
  real log_rt0_raw; //initial log r0
  real log_incid_rate_raw; // standard deviation for incidence
  real<lower=0> seed_incid_raw[prev_vals];
  real <lower=0> log_rho; // raw underreporting parameter for cases as a function of tests
  real<lower=0> i_raw[n]; //incidence transformed so that the mean is always 1
  real<lower = 0> kappa;  // overdispersion parameter for cases
  real<lower=0> exp_rate;
  // priors for OU parameters
  real theta;
  real log_delta; 
  real mu; // mean that OU process prior reverts to
}

transformed parameters{
  vector[n] log_rt;
  real<lower=0> delta = exp(log_delta_mu + log_delta_sd*log_delta);
  real<lower=0> incid_rate = exp(log_incid_rate_raw*log_incid_rate_sd + 
                                 log_incid_rate_mean);
  vector<lower=0>[prev_vals] seed_incid;
  real<lower=0> incid[n]; //true unobserved incidence, transformed from raw
  real<lower=0> weighted_sum[n];
  real <lower=0> delay_sum[n];
  vector<lower=0>[prev_vals+n] aug_incid;
  real <lower=0> rho = exp(log_rho*log_rho_sd + log_rho_mu);
  real <lower=0> sigma = sqrt(2*theta*delta^2);
  
  
    
  // log_rt prior is OU process
  log_rt[1] = log_r0_mu + log_rt0_raw*log_r0_sd;
  for (r in 2:n) {
    log_rt[r] = log_rt[r-1]*exp(-theta) + mu*(1-exp(-theta)) + z[r-1]*sqrt( ( sigma^2/(2*theta) )*( 1 - exp(-2*theta) ) );
  }


  //seed incid
  for (i in 1:prev_vals){
    seed_incid[i] = seed_incid_raw[i]*inv(exp_rate);
  }

  
  //putting incidence fully into latent land
  //starting with aug_incid first 8 values being the seed_incids
  for (a in 1:prev_vals) {
      aug_incid[a] = seed_incid[a];
    }
  
  //true incidence modeled in terms of previous true incidence and seed incidence
  { int c;
    c = prev_vals;
    for (t in 1:n) {
      weighted_sum[t] = 0;
    
      if (c <= d ) {
        
        weighted_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(w, c)));

        }
    
      if (c > d) {

        weighted_sum[t] = dot_product(segment(aug_incid, c-d +1, d), 
                                      my_reverse(head(w, d)));
          
    
        }
    //create incid and aug incid
      incid[t] = exp(log_rt[t])*weighted_sum[t]*i_raw[t];
      aug_incid[t+prev_vals] = incid[t];
      c = c+1;
      }
    }
    
      { int c;
    //changing the indexing on this so that it allows for weight on the same day as 
    //the observation
    c = prev_vals +1;
    for (t in 1:n) {
      delay_sum[t] = 0;
    
      if (c <= d ) {
        
        delay_sum[t] = dot_product(head(aug_incid, c), 
                                       my_reverse(head(delay_weights, c)));

        }
    
      if (c > d) {

        delay_sum[t] = dot_product(segment(aug_incid, c-d +1, d), 
                                      my_reverse(head(delay_weights, d)));
          
    
        }
      c = c+1;
      }
    }
}
// The model to be estimated
model {

  //prior on underreporting
  log_rho ~ normal(0,1);

  //prior on transformed kappa
  kappa ~ normal(kappa_mu,kappa_sd) T[0,];

  //priors for log R0
  log_rt0_raw ~ normal(0, 1);
  //prior on sd of incidence
  log_incid_rate_raw ~ normal(0, 1);
  
    //seed incidence
  exp_rate ~ exponential(0.3);
  
  for (i in 1:prev_vals) {
      seed_incid_raw[i] ~ exponential(1);

  }
  
  // OU process param priors
  mu ~ normal(0,1);
  
  log_delta ~ normal(0,1);
  
  theta ~ exponential(1/(theta_mu));

  
  // std normal to be scaled to OU process
  for (r in 1:(n-1)) {
    z[r] ~ normal(0,1);
  }
  

  //raw incidence modeled in terms of previous true incidence and seed incidence
  for (t in 1:n) {
    i_raw[t] ~ gamma(exp(log_rt[t])*weighted_sum[t]*incid_rate, 
                exp(log_rt[t])*weighted_sum[t]*incid_rate);

  }


  //observed incidence model
      for (t in 1:n) {
    obs[t]~neg_binomial_2(rho*test[t]*delay_sum[t], kappa);

  }

  


}

generated quantities{
  int gen_obs[n];
  for (t in 1:n) {
      gen_obs[t] = neg_binomial_2_rng(rho*test[t]*delay_sum[t], kappa);
  }

}


