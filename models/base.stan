data {
  int<lower=1> R; // number of recipients
  int<lower=1> N_tot; // total number of observations
  int<lower=1> N[R]; // number of observations per country
  int<lower=0> cum_N[R+1]; // start and end index of country time series
  int<lower=1> max_N; // maximum number of years for a country time series
  int<lower=0> model_start; // year at which to start modeling
  int<lower=model_start> Y; // number of years
  int<lower=0> hiv_new[N_tot]; // hiv_new infections
  vector<lower=0>[N_tot] hiv_new_lower; // lower of plausibility estimates for new hiv infections
  vector<lower=0>[N_tot] hiv_new_upper; // upper of plausibility estimates for new hiv infections
  vector<lower=0>[N_tot] hiv_living; // hiv prevalent infections
  vector<lower=0>[N_tot] hiv_living_lower; // lower of hiv prevalent infections
  vector<lower=0>[N_tot] hiv_living_upper; // upper of hiv prevalent infections
  vector<lower=0>[N_tot] pop; // size of population
  int<lower=-model_start+1> year[N_tot]; // year indicator: 1 is the first year data is modelled (i.e. y0 + model_start)
}

transformed data {
  vector<lower=0>[N_tot] uninfected = 1 - hiv_living .* inv(pop); // proportion of population not infected with HIV
  vector<lower=0>[N_tot] hiv_new_range = 0.5 * (log(hiv_new_upper) - log(hiv_new_lower)); // range of hiv estimates
  vector<lower=0>[N_tot] hiv_living_range = 0.5 * (log(hiv_living_upper) - log(hiv_living_lower)); // range of hiv estimates
  int max_N_target = max_N - model_start; // maximum number of observations per country used for modeling
  int N_tot_target = N_tot - model_start * R; // total number of observations used for modeling
  int N_target[R]; // number of observations by country used for modeling 
  int cum_N_target[R+1]; // start and end index for time series used for modeling
  cum_N_target[1] = 0;
  for (r in 1:R) {
    N_target[r] = N[r] - model_start;
    cum_N_target[r+1] = cum_N[r+1] - model_start * r; 
  }
}

parameters {
  real alpha; // exp(alpha) population-average proportion of new infections 
  vector[R] alpha_r; // country-specific proportion of new infections
  real<lower=0> tau; // between-group (country and year) variation for alpha_r, alpha_y, and beta_r
  vector[N_tot_target] sigma_hiv_new_me; // measurement error in infections
  vector[N_tot_target] sigma_hiv_living_me; // measurement error in living
  real<lower=0> phi_inv; // inverse over-dispersion
}

transformed parameters {
  vector[N_tot_target] mu_hiv_new; // expected number of new infections
  vector[N_tot_target] hiv_new_me; // expected number of new infections with measurement error
  vector[N_tot_target] hiv_living_me; // expected number of living with measurement error
  real phi = inv_square(phi_inv); // over-dispersion 

  // loop model over countries
  {
    int counter_data = 0; // counter for the input data an
    int counter_target = 0; // counter for modeled data
    for (r in 1:R) {
      for (i in 1:N_target[r]) {
        counter_data = cum_N[r] + i + model_start;
        counter_target = cum_N_target[r] + i;
        // estimated living with HIV incorporating measurement error
        hiv_living_me[counter_target] = exp(log(hiv_living[counter_data])+sigma_hiv_living_me[counter_target]);
        // estimate expected number of new infections
        // - with year RE: alpha_y[year[counter_data]] 
        mu_hiv_new[counter_target] = uninfected[counter_data] * hiv_living_me[counter_target] * exp( alpha + alpha_r[r]);
        // ... incorporate measurement error
        hiv_new_me[counter_target] = exp(log(mu_hiv_new[counter_target])+sigma_hiv_new_me[counter_target]);
      }
    }
  }
}


model {
  // priors
  alpha ~ student_t(7., 0., 10.);
  alpha_r ~ normal(0., tau);
  tau ~ student_t(4., 0., 1.);
  phi ~ normal(0., 1.);
  
  // likelihood
  for (r in 1:R) {
    // measurement error
    sigma_hiv_living_me[(cum_N_target[r]+1):(cum_N_target[r+1])] ~ normal(0., hiv_living_range[(cum_N[r]+1+model_start):(cum_N[r+1])]);
    sigma_hiv_new_me[(cum_N_target[r]+1):(cum_N_target[r+1])] ~ normal(0., hiv_new_range[(cum_N[r]+1+model_start):(cum_N[r+1])]);
    // hiv new ~ NegBinom()
    hiv_new[(cum_N[r]+1+model_start):(cum_N[r+1])] ~ neg_binomial_2(hiv_new_me[(cum_N_target[r]+1):(cum_N_target[r+1])], phi);
  }
}


generated quantities {
  // log-likelihood
  vector[N_tot_target] log_lik;
  // expected number of new infections with and without funding
  {
    int counter_data = 0; // counter for the input data an
    int counter_target = 0; // counter for modeled data
    for (r in 1:R) {
      for (i in 1:N_target[r]) {
        counter_data = cum_N[r] + i + model_start;
        counter_target = cum_N_target[r] + i;
         // log-likelihood
        log_lik[counter_target] = neg_binomial_2_log_lpmf(hiv_new[counter_data] | log(hiv_new_me[counter_target]), phi);
      }
    }
  }
}
