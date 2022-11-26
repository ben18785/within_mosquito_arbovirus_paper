functions {
  real logistic_curve(real dilution, real a, real k, real b, real q) {
    
    real log_concentration = log(1.0 / dilution);
    return(a + (k - a)/(1 + q * exp(-b * log_concentration)));
  }
  
  real[] v_pop_ode(real t, // time
                   real[] y, // states
                   real[] theta, // parameters
                   real[] d_r, // data
                   int[] d_i // data
  ){
    real gamma = theta[1]; // per-capita degradation rate parameter
    real k_lm = theta[2]; // flow of virus through lumen parameter
    real a = theta[3]; // flow of virus through lumen parameter
    real alpha_m = theta[4]; //
    real k_m = theta[5];
    real k_mh = theta[6];
    real x_star = theta[7];
    real eta = theta[8];
    real x_0 = theta[9];
    real x_r = theta[10];
    real alpha_h = theta[11];
    real k_h = theta[12];
    
    real t_refeed = d_r[1];
    real x_t; // permeability of the membrane
    real dy_dt[3];
    real l = y[1];
    real m = y[2];
    real h = y[3];
    
    x_t = t<t_refeed? x_star + exp(-eta*t) * (x_0 - x_star): x_star + exp(-eta*(t - t_refeed)) * (x_r - x_star);
    
    dy_dt[1] = -gamma * l - x_t * (k_lm * l^2) / (a^2 + l^2); // lumen
    dy_dt[2] = x_t * (k_lm * l^2) / (a^2 + l^2) + alpha_m * m * (1 - m / k_m) - k_mh * (x_t - x_star) * m; // midgut
    dy_dt[3] = k_mh * (x_t - x_star) * m + alpha_h * h * (1 - h/k_h); // disseminated
    
    return dy_dt;
  }
  
  real [, ] solve_ode_single_dilution(real [] y_init, real t0,
        real [] ts, real [] theta, real [] d_r, int [] d_i, int n_unq_t,
        int n_difeq) {
          real y_hat[n_unq_t, n_difeq];
          y_hat = integrate_ode_bdf(v_pop_ode, y_init, t0, ts, theta, d_r, d_i);
          return(y_hat);
  }
  
    real [] solve_ode_multiple_dilution(real [] y_init, real t0,
        real [] ts, real [] theta, real [] d_r, int [] d_i, int n_unq_t,
        int n_difeq, int n_dilutions, real[] dilutions, real l0, 
        int n_obs, int[] t_ind, int[] difeq_ind, int[] dilution_ind,
        real[] y) {
          
          real y_hat_[n_unq_t, n_difeq,n_dilutions];
          real y_init_[n_difeq,n_dilutions];
          real mu_[n_obs];
          real pdf_out[n_obs];
          
          for(i in 1:n_dilutions){
            y_init_[1, i] = l0 * 1/dilutions[i];
            y_init_[2, i] = 0;
            y_init_[3, i] = 0;
            y_hat_[,,i] = solve_ode_single_dilution(y_init_[,i], t0, ts, theta, d_r, d_i, n_unq_t, n_difeq);
            }
            
            for(i in 1:n_obs){
              mu_[i] = log(y_hat_[t_ind[i], difeq_ind[i], dilution_ind[i]]);
              pdf_out[i] = lognormal_lpdf(y[i]|mu_[i], 1);
              }
            
            return(pdf_out);
    }
}

data {
  // to run the ODE
  int<lower=1> n_unq_t; // number of unique times
  int<lower=1> n_theta; // number of parameters
  int<lower=1> n_difeq; // number of differential equations
  int<lower=1> n_dilutions; // number of dilutions
  
  real t0; // initial time
  real ts[n_unq_t]; // time-points
  real d_r_in; // time of second feed
  real dilutions[n_dilutions]; // dilution values
  
  // dissection data
  int<lower=1> n_obs; // number of data points
  real y[n_obs]; // viral titres
  int t_ind[n_obs]; // indices of the times of dissection
  int difeq_ind[n_obs]; // indices of the mosquito dissection
  int dilution_ind[n_obs]; // indices of the dilutions
  
  // fixed parameters
  real x_0;
  real x_r;
  // need two titer lower bounds since we rescale midgut
  // and leg infection separately
  real titer_lower_bound[2];
  
  // data for binary dissemination outcomes
  int n_unq_t_binary;
  real ts_binary[n_unq_t_binary];
  int n_experiment_types_binary;
  int dilutions_binary[n_experiment_types_binary];
  real t_refeed_binary[n_experiment_types_binary];
  real refeed_amount_binary[n_experiment_types_binary];
  int<lower=1> n_obs_binary;
  int t_ind_binary[n_obs_binary];
  int n_infected_binary[n_obs_binary];
  int n_dissected_binary[n_obs_binary];
  int ind_binary[n_obs_binary];
  int n_refeed_types;
  real t_refeeds_unique[n_refeed_types];
  real refeed_amounts_unique[n_refeed_types];
  int difeq_ind_binary[n_obs_binary];
  
  // switches for selecting likelihoods to include
  int include_continuous;
  int include_continuous_censored;
  int include_binary;
  int include_hurdle;
  
  // chp damage data
  int n_chp;
  real chp[n_chp];
  real time_chp[n_chp];
  
  // posterior predictive simulations
  int n_dilutions_sim;
  int dilutions_sim[n_dilutions_sim];
  int n_g_t; // number of times
  real g_t[n_g_t]; // times
  int n_dilutions_sim_fine; // for day 5 dilutions only
  real dilutions_sim_fine[n_dilutions_sim_fine];
}

transformed data{
  int d_i[0];
  real d_r[1];
  int n_experiment = 2;
  d_r[1] = d_r_in;
}

parameters {
  // all parameters can be positive
  real<lower=0> l0;
  
  real<lower=0> gamma; // per-capita degradation rate parameter
  real<lower=0> k_lm; // flow of virus through lumen parameter
  real<lower=0> a; // flow of virus through lumen parameter
  real<lower=0> alpha_m; //
  real<lower=0> k_m;
  real<lower=0> k_mh;
  real<lower=0> alpha_h;
  real<lower=0> k_h;
  
  real<lower=0> sigma[2]; // standard deviation
  
  // hurdle parameters
  real<lower=0, upper=1> phi_d; // prob(escape midgut)
  
  // experiment effect on initial conditions
  real<lower=0> zeta;
  
  // determine probability of infection as function of concentration
  real<lower=0, upper=1> b1;
  real<lower=0, upper=1> b2;
  real<lower=0> b3;
  real<lower=0> b4;
  real<lower=0, upper=x_0> x_star;
  
  // CHP data parameters
  real<lower=0> eta;
  positive_ordered[2] chp_vals;
  real<lower=0> chp_sigma;
  
}

transformed parameters{
  
  real theta[n_theta];
  real mu[n_obs];
  real mu_binary[n_obs_binary];
  real<lower=0> sigma_in[n_obs];
  
  real y_init[n_difeq, n_dilutions];
  real y_hat[n_unq_t, n_difeq, n_dilutions];
  real y_hat_binary[n_unq_t_binary, n_difeq, n_experiment_types_binary];
  real y_init_binary[3];
  y_init_binary[2] = 0.0;
  y_init_binary[3] = 0.0;
  
  theta[1] = gamma; // per-capita degradation rate parameter
  theta[2] = k_lm ; // flow of virus through lumen parameter
  theta[3] = a; // flow of virus through lumen parameter
  theta[4] = alpha_m;
  theta[5] = k_m;
  theta[6] = k_mh;
  theta[7] = x_star;
  theta[8] = eta;
  theta[9] = x_0;
  theta[10] = x_r;
  theta[11] = alpha_h;
  theta[12] = k_h;
  
  // non-binary experiments
  for(i in 1:n_dilutions){
        y_init[1, i] = l0 * 1/dilutions[i];
        y_init[2, i] = 0;
        y_init[3, i] = 0;
        y_hat[,,i] = solve_ode_single_dilution(
            y_init[,i], t0, ts, theta, d_r, d_i, n_unq_t, n_difeq);
  }
  
  for(i in 1:n_obs){
    mu[i] = y_hat[t_ind[i], difeq_ind[i], dilution_ind[i]] > 0 ? log(y_hat[t_ind[i], difeq_ind[i], dilution_ind[i]]): log(1E-10);
    sigma_in[i] = sigma[(difeq_ind[i]-1)];
  }
  
  // binary experiments
  for(i in 1:n_experiment_types_binary){
    real refeed_amount_temp = refeed_amount_binary[i];
    real theta_temp[n_theta] = theta;
    theta_temp[10] = refeed_amount_temp;
    theta_temp[11] = alpha_h;
    y_init_binary[1] = zeta * l0 * 1 / dilutions_binary[i];
    y_hat_binary[, , i] = solve_ode_single_dilution(
          y_init_binary, t0, ts_binary, theta_temp, {t_refeed_binary[i]}, d_i, n_unq_t_binary, n_difeq);
  }
  
  for(i in 1:n_obs_binary) {
    real y_temp = y_hat_binary[t_ind_binary[i], difeq_ind_binary[i], ind_binary[i]];
    mu_binary[i] = y_temp > 0 ? log(y_temp): log(1E-10);
  }
  
}

model {
  real phi;
  
  // priors
  l0 ~ normal(0, 1);
  gamma ~ normal(0, 10); // per-capita degradation rate parameter
  k_lm ~ normal(0, 10); // flow of virus through lumen parameter
  a ~ normal(0, 10); // flow of virus through lumen parameter
  alpha_m ~ normal(3, 10); //
  k_m ~ normal(1, 10);
  k_mh ~ normal(0, 10);
  alpha_h ~ normal(1.5, 10);
  k_h ~ normal(1, 10);
  sigma[1] ~ cauchy(0, 5);
  sigma[2] ~ cauchy(0, 5);
  zeta ~ normal(1, 1);
  b1 ~ normal(0.1, 0.1);
  b2 ~ normal(0.94, 0.5);
  b3 ~ normal(3, 1);
  b4 ~ normal(0, 0.005);
  eta ~ normal(2, 1);
  x_star ~ normal(0.1, 0.5);
  
  // priors for chp damage
  chp_vals[1] ~ normal(300, 200);
  chp_vals[2] ~ normal(3000, 3000);
  chp_sigma ~ normal(0, 200);
  
  // likelihood for semi-non-binary data
  if(include_continuous == 1) {
    if(include_hurdle == 1) {
      for(i in 1:n_obs) {
        phi = logistic_curve(dilutions[dilution_ind[i]], b1, b2, b3, b4);
        if(difeq_ind[i] == 3) // legs
          phi *= phi_d;
        if(y[i] > 0){
          target += log(phi);
          y[i] ~ lognormal(mu[i], sigma_in[i]) T[titer_lower_bound[difeq_ind[i] - 1],];
        }else {
          real log_prob[2];
          log_prob[1] = log1m(phi);
          log_prob[2] = log(phi) + lognormal_lcdf(titer_lower_bound[difeq_ind[i] - 1]| mu[i], sigma_in[i]);
          target += log_sum_exp(log_prob);
        }
      }
    }else {
      for(i in 1:n_obs) {
        if(y[i] >= 0){
          y[i] ~ lognormal(mu[i], sigma_in[i]);
        } else {
          if(include_continuous_censored == 1)
            target += lognormal_lcdf(titer_lower_bound[difeq_ind[i] - 1]| mu[i], sigma_in[i]);
        }
      }
    }
  }
  
  // likelihood for binary data for disseminated infection
  if(include_binary == 1) {
    for(i in 1:n_obs_binary) {
      int ind_temp = difeq_ind_binary[i] - 1;
      real titer_lowerb_temp = titer_lower_bound[ind_temp];
      real p = 1.0 - lognormal_cdf(titer_lowerb_temp, mu_binary[i], sigma[ind_temp]);
      if(include_hurdle == 1) {
        // include zeta here since this yields dilution for these data
        real phi_m = logistic_curve(1 / zeta * dilutions_binary[ind_binary[i]], b1, b2, b3, b4);
        if(ind_temp == 1) // midgut
          p *= phi_m;
        else if(ind_temp == 2) // legs
          p *= phi_m * phi_d;
      }
      n_infected_binary[i] ~ binomial(n_dissected_binary[i], p);
    }
  }
  
  // likelihood for CHP damage data
  for(i in 1:n_chp) {
    real mu_tmp = chp_vals[1] + (chp_vals[2] - chp_vals[1]) * exp(-eta * time_chp[i]);
    chp[i] ~ lognormal(log(mu_tmp), chp_sigma);
  }
}

generated quantities{
  real pred_y_hat[n_experiment, n_refeed_types, n_g_t, n_difeq, n_dilutions_sim];
  real dose_response_midgut[n_experiment, n_dilutions_sim_fine];
  
  for(k in 1:n_experiment) {
   for(i in 1:n_dilutions_sim){
     for(j in 1:n_refeed_types) {
       real inits[3];
       real theta_temp[n_theta] = theta;
       inits[1] = l0 * 1 / dilutions_sim[i];
       inits[2] = 0;
       inits[3] = 0;
       theta_temp[10] = refeed_amounts_unique[j];
       if(k == 2) {
        inits[1] *= zeta;
       }
       pred_y_hat[k,j,,,i] = integrate_ode_bdf(v_pop_ode, inits, t0, g_t, theta_temp, {t_refeeds_unique[j]}, d_i);
     }
   }
  }
  
  // dose response curve at day 5
  for(k in 1:n_experiment) {
    for(i in 1:n_dilutions_sim_fine) {
      real temp[1, 3];
      real inits[3];
      inits[1] = l0 * 1 / dilutions_sim_fine[i];
      inits[2] = 0;
      inits[3] = 0;
      if(k == 2) {
          inits[1] *= zeta;
      }
      temp = integrate_ode_bdf(v_pop_ode, inits, t0, {5.0}, theta, {100.0}, d_i);
      dose_response_midgut[k, i] = temp[1, 2]; // get midgut result only
    }
  }
  
}
