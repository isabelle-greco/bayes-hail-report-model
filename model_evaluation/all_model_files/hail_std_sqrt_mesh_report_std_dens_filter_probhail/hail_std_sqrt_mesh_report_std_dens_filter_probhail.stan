// Testing new set up

// input data
data {
  int<lower=0> N_train;
  int<lower=0> N_test;
  matrix[N_train, 2] X_mesh_train;
  matrix[N_test, 2] X_mesh_test;
  matrix[N_train, 2] X_report_train;
  matrix[N_test, 2] X_report_test;
  int<lower=0, upper=1> report_train[N_train];
  int<lower=0, upper=1> report_test[N_test];
  vector[2] beta_hail_lower_bounds;
  vector[2] beta_report_lower_bounds;
}
transformed data {
  // calculate statistics
  real mesh_train_mean = mean(X_mesh_train[, 2]);
  real mesh_train_sd = sd(X_mesh_train[, 2]);
  real dens_train_mean = mean(X_report_train[, 2]);
  real dens_train_sd = sd(X_report_train[, 2]);
  
  // transformed matrices
  matrix[N_train, 2] X_mesh_train_std;
  matrix[N_train, 2] X_report_train_std;
  matrix[N_test, 2] X_mesh_test_std;
  matrix[N_test, 2] X_report_test_std;

  // filling in matrices
  X_mesh_train_std[, 1] = X_mesh_train[, 1];
  X_mesh_train_std[, 2] = (X_mesh_train[, 2] - mesh_train_mean) / mesh_train_sd;
  X_mesh_test_std[, 1] = X_mesh_test[, 1];
  X_mesh_test_std[, 2] = (X_mesh_test[, 2] - mesh_train_mean) / mesh_train_sd;

  X_report_train_std[, 1] = X_report_train[, 1];
  X_report_train_std[, 2] = (X_report_train[, 2] - dens_train_mean) / dens_train_sd;
  X_report_test_std[, 1] = X_report_test[, 1];
  X_report_test_std[, 2] = (X_report_test[, 2] - dens_train_mean) / dens_train_sd;
}
// parameters
parameters {
  vector<lower=beta_hail_lower_bounds>[2] beta_hail_raw;
  vector<lower=beta_report_lower_bounds>[2] beta_report_raw;
}
// transformed parameters
transformed parameters {
  vector<lower=0, upper=1>[N_train] prob_report_train;
  vector[2] beta_hail;
  vector[2] beta_report;
  real<lower=0> mesh_at_prob_0p5;
  {
    // can't have constraints in these local blocks
    vector[N_train] prob_hail_train = inv_logit(X_mesh_train_std * beta_hail_raw);
    vector[N_train] prob_report_given_hail_train = inv_logit(X_report_train_std * beta_report_raw);
    prob_report_train = prob_hail_train .* prob_report_given_hail_train;
  }
  beta_hail[1] = beta_hail_raw[1] - (beta_hail_raw[2] * (mesh_train_mean / mesh_train_sd));
  beta_hail[2] = beta_hail_raw[2] / mesh_train_sd;
  beta_report[1] = beta_report_raw[1] - (beta_report_raw[2] * (dens_train_mean / dens_train_sd));
  beta_report[2] = beta_report_raw[2] / dens_train_sd;
  mesh_at_prob_0p5 = square(-beta_hail[1] / beta_hail[2]); // need to square as mesh is sqrt
}
// likelihood and priors
model {
  // likelihoood
  target += bernoulli_lpmf(report_train | prob_report_train);
  // priors
  target += normal_lpdf(beta_hail_raw[1] | -2.0, 8.0);
  target += normal_lpdf(beta_hail_raw[2] | 0.5, 1.0);
  target += normal_lpdf(beta_report_raw[1] | -2.0, 8.0);
  target += normal_lpdf(beta_report_raw[2] | 0.5, 1.0);
  // stronger prior information on mesh leading to 0.5 prob of hail
  target += skew_normal_lpdf(mesh_at_prob_0p5 | 25.0, 8.0, 4.0);
}
// generated quantities (e.g. post predictive check)
generated quantities {
  // delcaration needs to be before assignment later 
  vector[N_test] log_lik_test;
  vector<lower=0, upper=1>[N_test] prob_report_test;
  int<lower=0, upper=1> sim_report_test[N_test];
  // probability of report 
  {
    vector[N_test] prob_hail_test = inv_logit(X_mesh_test_std * beta_hail_raw);
    vector[N_test] prob_report_given_hail_test = inv_logit(X_report_test_std * beta_report_raw);
    prob_report_test = prob_hail_test .* prob_report_given_hail_test;
  }
  // generating posterior predictive samples
  sim_report_test = bernoulli_rng(prob_report_test);
  // calculating log-likelihood
  for (i in 1:N_test){
    log_lik_test[i] = bernoulli_lpmf(report_test[i] | prob_report_test[i]);
  }
}

