functions {
  // Discretization of the cumulative normal distribution function
  real diff_norm(real x, real max_x, real mu, real sigma) {
    {
      if (x == 0) {
        return (normal_cdf(1,  mu, sigma) - normal_cdf(0,  mu, sigma)) / (normal_cdf(max_x, mu, sigma) - normal_cdf(0, mu, sigma));
      } else {
        return (normal_cdf(x+1, mu, sigma) - normal_cdf(x, mu, sigma)) / (normal_cdf(max_x, mu, sigma) - normal_cdf(0, mu, sigma));
      }
    }
  }
}

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
  int<lower=1> D; // number of funding variables
  int<lower=1> K; // number of control variables
  matrix<lower=0>[N_tot,D] funding; // std control funding: [1] foreign and [2] domestic
  matrix[N_tot,K] Z; // control variables: [1] gdp, [2] edu, [3] maternal, [4] pop density
  int<lower=-model_start+1> year[N_tot]; // year indicator: 1 is the first year data is modelled (i.e. y0 + model_start)
  vector<lower=0>[4] p_a_hyper; // activity distribution hyperparameters
  vector[2] p_a_hyper_diff; // difference between activity distribution hyperparameters
}

transformed data {
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
  real mu_p_a; // expected activity delay
  real<lower=0> sigma_p_a; // standard deviation in activity delay
  real diff_mu_p_a; // difference between funding types in activity delay
  real alpha; // exp(alpha) population-average proportion of new infections 
  vector[R] alpha_r; // country-specific proportion of new infections
  vector[Y] alpha_y; // year-specific proportion of new infections
  vector<lower=0>[2] tau; // between-group (country and year) variation for alpha_r, alpha_y, and beta_r
  real theta; // effect of STD control funding
  real theta_share;
  vector[K] zeta; // effect of control variables
  matrix<lower=0>[N_tot_target,D] active_funding; // active funding
  vector[N_tot_target] sigma_hiv_new_me; // measurement error in infections
  vector[N_tot_target] sigma_hiv_living_me; // measurement error in living
  vector<lower=0>[D+1] phi_inv; // inverse over-dispersion
}

transformed parameters {
  vector[D] mu_p_a_d; // expected activity delay by funding type
  real p_a[D,model_start+1]; // activity distribution
  vector[N_tot_target] mu_hiv_new; // expected number of new infections
  vector[N_tot_target] hiv_new_me; // expected number of new infections with measurement error
  vector[N_tot_target] hiv_living_me; // expected number of living with measurement error
  vector[N_tot_target] uninfected_me; // uninfected with measurement error
  matrix[N_tot_target,D] mu_active_funding; // expected active funding
  matrix[N_tot_target,D] sigma_active_funding; // expected active funding
  vector[D+1] phi = inv_square(phi_inv); // over-dispersion 
  if (D > 1) {
    mu_p_a_d[1] = mu_p_a + diff_mu_p_a;
    mu_p_a_d[2] = mu_p_a - diff_mu_p_a;
  } else {
    mu_p_a_d[1] = mu_p_a + diff_mu_p_a;
  }

  // compute activity distribution
  for (r in 1:R) {
    for (d in 1:D) {
      for (s in 0:model_start) p_a[d,s+1] = diff_norm(model_start-s, model_start, mu_p_a_d[d], sigma_p_a);
    }
  }
  
  
  // loop model over countries
  {
    int counter_data = 0; // counter for the input data an
    int counter_target = 0; // counter for modeled data
    for (r in 1:R) {
      for (i in 1:N_target[r]) {
        counter_data = cum_N[r] + i + model_start;
        counter_target = cum_N_target[r] + i;
        // compute active funding
        for (d in 1:D) {
          mu_active_funding[counter_target,d] = dot_product(funding[(counter_data-model_start):counter_data,d], to_vector(p_a[d]));
          sigma_active_funding[counter_target,d] = sqrt(mu_active_funding[counter_target,d] * (1 + mu_active_funding[counter_target,d] / phi[d]));
        }
        // estimated living with HIV incorporating measurement error
        hiv_living_me[counter_target] = exp(log(hiv_living[counter_data])+sigma_hiv_living_me[counter_target]);
        uninfected_me[counter_target] = 1 - hiv_living_me[counter_target] .* inv(pop[counter_data]);
        // estimate expected number of new infections
        mu_hiv_new[counter_target] = uninfected_me[counter_target] * hiv_living_me[counter_target] * exp( alpha + alpha_r[r] + alpha_y[year[counter_data]] + (1 / 1e3 * sum(active_funding[counter_target]) / hiv_living[counter_data]) * (theta) + theta_share * active_funding[counter_target,1] / sum(active_funding[counter_target]) + Z[counter_data] * zeta);
        // ... incorporate measurement error
        hiv_new_me[counter_target] = exp(log(mu_hiv_new[counter_target])+sigma_hiv_new_me[counter_target]);
      }
    }
  }
}


model {
  // priors
  mu_p_a ~ normal(p_a_hyper[1], p_a_hyper[2]);
  sigma_p_a ~ normal(p_a_hyper[3], p_a_hyper[4]);
  diff_mu_p_a ~ normal(p_a_hyper_diff[1], p_a_hyper_diff[2]);
  alpha ~ student_t(7., 0., 10.);
  theta ~ student_t(4., 0., 2.5);
  theta_share ~ student_t(4., 0., 2.5);
  zeta ~ student_t(4., 0., 2.5);
  alpha_r ~ normal(0., tau[1]);
  alpha_y ~ normal(0., tau[2]);
  tau ~ student_t(4., 0., 1.);
  phi ~ normal(0., 1.);
  
  // likelihood
  for (r in 1:R) {
    // estimated active funding
    for (d in 1:D) {
      active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d] ~ normal(mu_active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d], sigma_active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d]);
    }
    // measurement error
    sigma_hiv_living_me[(cum_N_target[r]+1):(cum_N_target[r+1])] ~ normal(0., hiv_living_range[(cum_N[r]+1+model_start):(cum_N[r+1])]);
    sigma_hiv_new_me[(cum_N_target[r]+1):(cum_N_target[r+1])] ~ normal(0., hiv_new_range[(cum_N[r]+1+model_start):(cum_N[r+1])]);
    // hiv new ~ NegBinom()
    hiv_new[(cum_N[r]+1+model_start):(cum_N[r+1])] ~ neg_binomial_2(hiv_new_me[(cum_N_target[r]+1):(cum_N_target[r+1])], phi[D+1]);
  }
}


generated quantities {
  // log-likelihood
  vector[N_tot_target] log_lik;
  // expected number of new infections with and without funding
  matrix[max_N_target,R] mu_hiv_new_none;
  real mu_hiv_new_funding[max_N_target,R] ;
  {
    int counter_data = 0; // counter for the input data an
    int counter_target = 0; // counter for modeled data
    for (r in 1:R) {
      for (i in 1:N_target[r]) {
        counter_data = cum_N[r] + i + model_start;
        counter_target = cum_N_target[r] + i;
         // log-likelihood
        log_lik[counter_target] = neg_binomial_2_log_lpmf(hiv_new[counter_data] | log(hiv_new_me[counter_target]), phi[D+1]);
        // expected number of new infections with/without std
        mu_hiv_new_none[i,r] = uninfected_me[counter_target] * hiv_living_me[counter_target] * exp( alpha + alpha_r[r] + Z[counter_data] * zeta);
        mu_hiv_new_funding[i,r] = uninfected_me[counter_target] * hiv_living_me[counter_target] * exp( alpha + alpha_r[r] + alpha_y[year[counter_data]] + (1 / 1e3 * sum(mu_active_funding[counter_target]) / hiv_living_me[counter_target]) * (theta) + theta_share * mu_active_funding[counter_target,1] / sum(mu_active_funding[counter_target]) + Z[counter_data] * zeta );
      }
    }
  }
}












