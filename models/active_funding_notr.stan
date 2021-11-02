functions {
  // Discretization of the cumulative normal distribution function
  real diff_norm(real x, real max_x, real mu, real sigma) {
    {
      if (x == 0) {
        return normal_cdf(1,  mu, sigma) / normal_cdf(max_x, mu, sigma);
      } else {
        return (normal_cdf(x+1, mu, sigma) - normal_cdf(x, mu, sigma)) / normal_cdf(max_x, mu, sigma);
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
  vector<lower=0>[N_tot] hiv_living; // hiv prevalent infections
  vector<lower=0>[N_tot] pop; // size of population
  int<lower=1> D; // number of funding variables
  int<lower=1> K; // number of control variables
  matrix<lower=0>[N_tot,D] funding; // std control funding: [1] foreign and [2] domestic
  matrix[N_tot,K] Z; // control variables: [1] gdp, [2] edu, [3] maternal, [4] pop density
  int<lower=-model_start+1> year[N_tot]; // year indicator: 1 is the first year data is modelled (i.e. y0 + model_start)
  vector<lower=0>[4] p_a_hyper;
}

transformed data {
  vector<lower=0>[N_tot] uninfected = 1 - hiv_living .* inv(pop); // proportion of population not infected with HIV
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
  real mu_p_a;
  real<lower=0> sigma_p_a;
  real diff_mu_p_a;
  real alpha; // exp(alpha) population-average proportion of new infections 
  vector[R] alpha_r; // country-specific proportion of new infections
  vector[Y] alpha_y; // year-specific proportion of new infections
  vector<lower=0>[2] tau; // between-group (country and year) variation for alpha_r, alpha_y, and beta_r
  real theta; // effect of STD control funding
  //vector[R] theta_r;
  real theta_share;
  vector[K] zeta; // effect of control variables
  matrix<lower=0>[N_tot_target,D] active_funding; // active funding
  vector<lower=0>[D+1] phi_inv; // inverse over-dispersion
}

transformed parameters {
  vector[2] mu_p_a_d;
  real p_a[D,model_start+1]; // activity distribution
  vector[N_tot_target] mu_hiv_new; // expected number of new infections
  matrix[N_tot_target,D] mu_active_funding; // expected active funding
  matrix[N_tot_target,D] sigma_active_funding; // expected active funding
  real phi_hiv_new; // over-dispersion 
  vector[D+1] phi = inv_square(phi_inv); // over-dispersion 
  mu_p_a_d[1] = mu_p_a + diff_mu_p_a;
  mu_p_a_d[2] = mu_p_a - diff_mu_p_a;
  
  // compute activity distribution
  for (r in 1:R) {
    for (d in 1:D) {
      for (s in 0:model_start) p_a[d,s+1] = diff_norm(model_start-s, model_start, mu_p_a, sigma_p_a);
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
        // estimate expected number of new infections
        mu_hiv_new[counter_target] = uninfected[counter_data] * hiv_living[counter_data] * exp( alpha + alpha_r[r] + alpha_y[year[counter_data]] + (1 / 1e3 * sum(active_funding[counter_target]) / hiv_living[counter_data]) * (theta) + theta_share * active_funding[counter_target,1] / sum(active_funding[counter_target]) + Z[counter_data] * zeta);
      }
    }
  }
}


model {
  // priors
  mu_p_a ~ normal(p_a_hyper[1], p_a_hyper[2]);
  sigma_p_a ~ normal(p_a_hyper[3], p_a_hyper[4]);
  diff_mu_p_a ~ normal(0.5, .5);
  alpha ~ student_t(7., 0., 10.);
  theta ~ student_t(4., 0., 2.5);
  theta_share ~ student_t(4., 0., 2.5);
  zeta ~ student_t(4., 0., 2.5);
  alpha_r ~ normal(0., tau[1]);
  alpha_y ~ normal(0., tau[2]);
  //theta_r ~ normal(0., tau[3]);
  tau[1] ~ student_t(4., 0., 1.);
  tau[2] ~ student_t(4., 0., 1.);
  //tau[3] ~ student_t(4., 0., 1.);
  phi ~ normal(0., 1.);
  
  // likelihood
  for (r in 1:R) {
    for (d in 1:D) {
      active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d] ~ normal(mu_active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d], sigma_active_funding[(cum_N_target[r]+1):(cum_N_target[r+1]),d]);
    }
    hiv_new[(cum_N[r]+1+model_start):(cum_N[r+1])] ~ neg_binomial_2(mu_hiv_new[(cum_N_target[r]+1):(cum_N_target[r+1])], phi[D+1]);
  }
}
