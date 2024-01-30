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
  matrix[N_train, 3] X_report_train;
  matrix[N_test, 3] X_report_test;
  int<lower=0, upper=1> report_train[N_train];
  int<lower=0, upper=1> report_test[N_test];
  vector[2] beta_hail_lower_bounds;
  vector[4] beta_report_lower_bounds;
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
  // note that third column filled in transformed parameters 
  
  X_report_test_std[, 1] = X_report_test[, 1];
  X_report_test_std[, 2] = (X_report_test[, 2] - dens_train_mean) / dens_train_sd;
  // again third column filled in transformed parameters
}
// parameters
parameters {
  vector<lower=beta_hail_lower_bounds>[2] beta_hail_raw;
  vector<lower=beta_report_lower_bounds>[4] beta_report_raw;
}
// transformed parameters
transformed parameters {
  // unstandardised parameters
  vector[2] beta_hail;
  vector[4] beta_report;

  // probability of report
  vector<lower=0, upper=1>[N_train] prob_report_train;
  
  // Yeo-Johnson transformed variables
  vector[N_train] mesh_train_yj_report;
  vector[N_test] mesh_test_yj_report;

  // applying transform
  for (i in 1:N_train) {
    mesh_train_yj_report[i] = yeo_johnson(X_report_train[i, 3], beta_report_raw[4]);
  }
  for (i in 1:N_test) {
    mesh_test_yj_report[i] = yeo_johnson(X_report_test[i, 3], beta_report_raw[4]);
  }
  
  // calculating statistics
  real mesh_train_yj_report_mean = mean(mesh_train_yj_report);
  real mesh_train_yj_report_sd = sd(mesh_train_yj_report);

  // making transformed and normalised matrics
  matrix[N_train, 3] X_report_train_std_trans = append_col(X_report_train_std, (mesh_train_yj_report - mesh_train_yj_report_mean) / mesh_train_yj_report_sd);
  matrix[N_test, 3] X_report_test_std_trans = append_col(X_report_test_std, (mesh_test_yj_report - mesh_train_yj_report_mean) / mesh_train_yj_report_sd);
  
  {
    // can't have constraints in these local blocks
    // note the yj parameter is in the last entry of beta_hail
    vector[N_train] prob_hail_train = inv_logit(X_mesh_train_std * beta_hail_raw);
    vector[N_train] prob_report_given_hail_train = inv_logit(X_report_train_std_trans * beta_report_raw[1:3]);
    prob_report_train = prob_hail_train .* prob_report_given_hail_train;
  }

  // unstandardised parameters  
  beta_report[1] = beta_report_raw[1] - (beta_report_raw[2] * (dens_train_mean / dens_train_sd)) - (beta_report_raw[3] * (mesh_train_yj_report_mean / mesh_train_yj_report_sd));
  beta_report[2] = beta_report_raw[2] / dens_train_sd;
  beta_report[3] = beta_report_raw[3] / mesh_train_yj_report_sd;
  beta_report[4] = beta_report_raw[4];
  // normalising not using yj mean/sd
  beta_hail[1] = beta_hail_raw[1] - (beta_hail_raw[2] * (mesh_train_mean / mesh_train_sd));
  beta_hail[2] = beta_hail_raw[2] / mesh_train_sd; 
}
// likelihood and priors
model {
  // likelihoood
  target += bernoulli_lpmf(report_train | prob_report_train);
  // priors
  target += normal_lpdf(beta_hail_raw[1] | -2.0, 4.0);
  target += gamma_lpdf(beta_hail_raw[2] | 6.0, 3.0);
  target += normal_lpdf(beta_report_raw[1] | -2.0, 4.0);
  target += gamma_lpdf(beta_report_raw[2] | 6.0, 3.0);
  target += normal_lpdf(beta_report_raw[3] | 0.0, 1.0);
  target += normal_lpdf(beta_report_raw[4] | 1.0, 2.0 / 3.0);
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
    vector[N_test] prob_report_given_hail_test = inv_logit(X_report_test_std_trans * beta_report_raw[1:3]);
    prob_report_test = prob_hail_test .* prob_report_given_hail_test;
  }
  // generating posterior predictive samples
  sim_report_test = bernoulli_rng(prob_report_test);
  // calculating log-likelihood
  for (i in 1:N_test){
    log_lik_test[i] = bernoulli_lpmf(report_test[i] | prob_report_test[i]);
  }
}
