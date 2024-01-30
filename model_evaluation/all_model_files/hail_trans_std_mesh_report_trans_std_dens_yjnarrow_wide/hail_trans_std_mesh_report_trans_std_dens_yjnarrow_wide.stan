// Testing new set up

functions {
  // implements the yeo_johnson transform 
  real yeo_johnson(real x, real lambda){
    real eps = 1e-5;
    if (x >= 0.0){
      if (abs(lambda) < eps) {
        return log1p(x);
      } else {
        return (((x + 1.0) .^ lambda) - 1.0) ./ lambda;
      }
    } else {
      if (abs(lambda - 2.0) < eps) {
        return -log1p(-x);
      } else {
        return -(((-x + 1.0) .^ (2.0 - lambda)) - 1.0) ./ (2.0 - lambda);
      }
    }
  }
}
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
  vector[3] beta_hail_lower_bounds;
  vector[3] beta_report_lower_bounds;
}
// parameters
parameters {
  vector<lower=beta_hail_lower_bounds>[3] beta_hail_raw;
  vector<lower=beta_report_lower_bounds>[3] beta_report_raw;
}
// transformed parameters
transformed parameters {
  // unstandardised parameters
  vector[3] beta_hail;
  vector[3] beta_report;

  // standardised and transformed matrices
  matrix[N_train, 2] X_mesh_train_std_trans;
  matrix[N_test, 2] X_mesh_test_std_trans;
  matrix[N_train, 2] X_report_train_std_trans;
  matrix[N_test, 2] X_report_test_std_trans;

  // probability of report
  vector<lower=0, upper=1>[N_train] prob_report_train;
  
  // Yeo-Johnson transformed variables
  vector[N_train] mesh_train_yj;
  vector[N_test] mesh_test_yj;
  vector[N_train] dens_train_yj;
  vector[N_test] dens_test_yj;

  // applying transform
  for (i in 1:N_train) {
    mesh_train_yj[i] = yeo_johnson(X_mesh_train[i, 2], beta_hail_raw[3]);
    dens_train_yj[i] = yeo_johnson(X_report_train[i, 2], beta_report_raw[3]);
  }
  for (i in 1:N_test) {
    mesh_test_yj[i] = yeo_johnson(X_mesh_test[i, 2], beta_hail_raw[3]);
    dens_test_yj[i] = yeo_johnson(X_report_test[i, 2], beta_report_raw[3]);
  }
  
  // calculating statistics
  real mesh_train_yj_mean = mean(mesh_train_yj);
  real mesh_train_yj_sd = sd(mesh_train_yj);
  real dens_train_yj_mean = mean(dens_train_yj);
  real dens_train_yj_sd = sd(dens_train_yj);  

  // making transformed and normalised matrics
  X_mesh_train_std_trans[, 1] = X_mesh_train[, 1];
  X_mesh_train_std_trans[, 2] = (mesh_train_yj - mesh_train_yj_mean) / mesh_train_yj_sd;

  X_mesh_test_std_trans[, 1] = X_mesh_test[, 1];
  X_mesh_test_std_trans[, 2] = (mesh_test_yj - mesh_train_yj_mean) / mesh_train_yj_sd;

  X_report_train_std_trans[, 1] = X_report_train[, 1];
  X_report_train_std_trans[, 2] = (dens_train_yj - dens_train_yj_mean) / dens_train_yj_sd;

  X_report_test_std_trans[, 1] = X_report_test[, 1];
  X_report_test_std_trans[, 2] = (dens_test_yj - dens_train_yj_mean) / dens_train_yj_sd;
  {
    // can't have constraints in these local blocks
    // note the yj parameter is in the last entry of beta_hail
    vector[N_train] prob_hail_train = inv_logit(X_mesh_train_std_trans * beta_hail_raw[1:2]);
    vector[N_train] prob_report_given_hail_train = inv_logit(X_report_train_std_trans * beta_report_raw[1:2]);
    prob_report_train = prob_hail_train .* prob_report_given_hail_train;
  }

  // unstandardised parameters  
  beta_report[1] = beta_report_raw[1] - (beta_report_raw[2] * (dens_train_yj_mean / dens_train_yj_sd));
  beta_report[2] = beta_report_raw[2] / dens_train_yj_sd;
  beta_report[3] = beta_report_raw[3];
  beta_hail[1] = beta_hail_raw[1] - (beta_hail_raw[2] * (mesh_train_yj_mean / mesh_train_yj_sd));
  beta_hail[2] = beta_hail_raw[2] / mesh_train_yj_sd;
  beta_hail[3] = beta_hail_raw[3];
}
// likelihood and priors
model {
  // likelihoood
  target += bernoulli_lpmf(report_train | prob_report_train);
  // priors
  target += normal_lpdf(beta_hail_raw[1] | -2.0, 8.0);
  target += normal_lpdf(beta_hail_raw[2] | 0.5, 1.0);
  target += normal_lpdf(beta_hail_raw[3] | 0.45, 0.08);
  target += normal_lpdf(beta_report_raw[1] | -2.0, 8.0);
  target += normal_lpdf(beta_report_raw[2] | 0.5, 1.0);
  target += normal_lpdf(beta_report_raw[3] | -0.17, 0.10);
}
// generated quantities (e.g. post predictive check)
generated quantities {
  // delcaration needs to be before assignment later 
  vector[N_test] log_lik_test;
  vector<lower=0, upper=1>[N_test] prob_report_test;
  int<lower=0, upper=1> sim_report_test[N_test];
  // probability of report 
  {
    vector[N_test] prob_hail_test = inv_logit(X_mesh_test_std_trans * beta_hail_raw[1:2]);
    vector[N_test] prob_report_given_hail_test = inv_logit(X_report_test_std_trans * beta_report_raw[1:2]);
    prob_report_test = prob_hail_test .* prob_report_given_hail_test;
  }
  // generating posterior predictive samples
  sim_report_test = bernoulli_rng(prob_report_test);
  // calculating log-likelihood
  for (i in 1:N_test){
    log_lik_test[i] = bernoulli_lpmf(report_test[i] | prob_report_test[i]);
  }
}
