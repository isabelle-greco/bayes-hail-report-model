// Testing new set up

functions {
  // implements the yeo_johnson transform 
  vector yeo_johnson(vector x, real lambda){
    if (lambda != 0){
      return (((x + 1.0) .^ lambda) - 1.0) ./ lambda;
    } else {
      return log1p(x);
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
  vector[2] beta_report_lower_bounds;
}
// parameters
parameters {
  vector<lower=beta_hail_lower_bounds>[3] beta_hail;
  vector<lower=beta_report_lower_bounds>[2] beta_report;
}
// transformed parameters
transformed parameters {
  matrix[N_train, 2] X_mesh_train_trans;
  matrix[N_test, 2] X_mesh_test_trans;
  vector<lower=0, upper=1>[N_train] prob_report_train;
  vector[N_train] mesh_train_yj = yeo_johnson(X_mesh_train[, 2], beta_hail[3]);
  vector[N_test] mesh_test_yj = yeo_johnson(X_mesh_test[, 2], beta_hail[3]);
  
  // making transformed matrics
  X_mesh_train_trans[, 1] = X_mesh_train[, 1];
  X_mesh_train_trans[, 2] = mesh_train_yj;

  X_mesh_test_trans[, 1] = X_mesh_test[, 1];
  X_mesh_test_trans[, 2] = mesh_test_yj;
  {
    // can't have constraints in these local blocks
    // note the yj parameter is in the last entry of beta_hail
    vector[N_train] prob_hail_train = inv_logit(X_mesh_train_trans * beta_hail[1:2]);
    vector[N_train] prob_report_given_hail_train = inv_logit(X_report_train * beta_report);
    prob_report_train = prob_hail_train .* prob_report_given_hail_train;
  }
}
// likelihood and priors
model {
  // likelihoood
  target += bernoulli_lpmf(report_train | prob_report_train);
  // priors
  target += normal_lpdf(beta_hail[1] | 0.0, 1.0);
  target += normal_lpdf(beta_hail[2] | 0.0, 1.0);
  target += normal_lpdf(beta_hail[3] | 0.0, 1.0);
  target += normal_lpdf(beta_report[1] | 0.0, 1.0);
  target += normal_lpdf(beta_report[2] | 0.0, 1.0);
}
// generated quantities (e.g. post predictive check)
generated quantities {
  // delcaration needs to be before assignment later 
  vector[N_test] log_lik_test;
  vector<lower=0, upper=1>[N_test] prob_report_test;
  int<lower=0, upper=1> sim_report_test[N_test];
  // probability of report 
  {
    vector[N_test] prob_hail_test = inv_logit(X_mesh_test_trans * beta_hail[1:2]);
    vector[N_test] prob_report_given_hail_test = inv_logit(X_report_test * beta_report);
    prob_report_test = prob_hail_test .* prob_report_given_hail_test;
  }
  // generating posterior predictive samples
  sim_report_test = bernoulli_rng(prob_report_test);
  // calculating log-likelihood
  for (i in 1:N_test){
    log_lik_test[i] = bernoulli_lpmf(report_test[i] | prob_report_test[i]);
  }
}

