### Script to perform k-fold validation given a stan file, path to the matrix
### data, and somewhere to write the results
###
### Last Modified: 2023-08-29
### Author: Isabelle Greco


# get args from command line
# 1 = library for r installs if need be
# 2 = path to entire data set 
# 3 = directory for matrix-ified data to which to apply the model
# 4 = stan file with model to use 
# 5 = number of MCMC chains
# 6 = number of iterations per chain (including warm up)
# 7 = directory for output files
# 8 = model name for output files
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 8) {
  stop(paste("Must give library for R package installs, location of all data",
	     "directory for matrix-ified data, stan file with model to use,",
	     "number of chains, number of iterations, directory for output",
	     "files, and name of model"), 
       call. = FALSE)
}


# load libraries + install if need be
for (package in c("tidyverse", "rstan", "janitor", "foreach", "doParallel", "loo", "bridgesampling", "posterior", "bayesplot")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


# set rstan options
rstan_options(auto_write = TRUE)

# set random seed for repeatability
set.seed(20230317)

# read in data
mesh_matrix <- readRDS(paste(args[3], "hail_matrix.rds", sep = "/"))
report_matrix <- readRDS(paste(args[3], "report_matrix.rds", sep = "/"))
beta_hail_lower_bounds <- readRDS(paste(args[3], "beta_hail_lower_bounds.rds", sep = "/"))
beta_report_lower_bounds <- readRDS(paste(args[3], "beta_report_lower_bounds.rds", sep = "/"))
model_data <- read_csv(args[2], col_types = "nnTnnnnnnncnfffffff") # folds as factor

# compiling the stan model 
compiled_model <- stan_model(file = args[4])

# rstan parameters
iter <- as.numeric(args[6])
chains <- as.numeric(args[5])
num_samples <- (iter / 2) * chains

if ("folds" %in% colnames(model_data)){
  print("Starting k-fold cross validation")
  # number of folds to iterate through
  num_folds <- mesh_matrix[, "folds"] %>% 
    max()
  
  # adding identifier 
  prepped_model_data <- model_data %>%
    add_column(unique_id = 1:nrow(.))
  prepped_mesh_matrix <- mesh_matrix %>%
    cbind(unique_id = 1:nrow(.))
  prepped_report_matrix <- report_matrix %>%
    cbind(unique_id = 1:nrow(.))

  # creating cluster
  cl <- makeCluster(num_folds)
  registerDoParallel(cl, cores = parallel::detectCores() - 4)
  
  # looping through folds 
  list_of_output <- foreach(k = 1:num_folds, .packages = c("tidyverse", "rstan", "janitor"), .verbose = TRUE) %dopar% {
    print(paste("Evaluating fold ", k, sep = " "))
    
    # separate data
    train_mesh_matrix <- prepped_mesh_matrix[prepped_mesh_matrix[, "folds"] != k, ]
    train_report_matrix <- prepped_report_matrix[prepped_report_matrix[, "folds"] != k, ]
    test_mesh_matrix <- prepped_mesh_matrix[prepped_mesh_matrix[, "folds"] == k, ]
    test_report_matrix <- prepped_report_matrix[prepped_report_matrix[, "folds"] == k, ]

    # as lists
    columns_to_select_mesh <- ! colnames(train_mesh_matrix) %in% c("folds", "unique_id", "report")
    columns_to_select_report <- ! colnames(train_report_matrix) %in% c("folds", "unique_id", "report")
    stan_data <- list(N_train = train_mesh_matrix %>% nrow(),
                      N_test = test_mesh_matrix %>% nrow(),
                      X_mesh_train = train_mesh_matrix[, columns_to_select_mesh],
                      X_report_train = train_report_matrix[, columns_to_select_report], 
                      X_mesh_test = test_mesh_matrix[, columns_to_select_mesh],
                      X_report_test = test_report_matrix[, columns_to_select_report],
		      # note: reports currently stored as column in the report matrix
		      report_train = train_report_matrix[, "report"],
		      report_test = test_report_matrix[, "report"],
		      beta_hail_lower_bounds = beta_hail_lower_bounds,
		      beta_report_lower_bounds = beta_report_lower_bounds)

    # fit model
    print(paste("Running model", k))
    fit <- sampling(compiled_model, data = stan_data, chains = chains, 
             iter = iter, cores = parallel::detectCores(),
	     pars = c("beta_hail", "beta_report", "log_lik_test", "sim_report_test", "prob_report_test"))

    print(paste("Extracting quantities model", k))
    # MCMC diagnostics
    fit %>%
      posterior::summarise_draws() %>%
      write_csv(paste(args[7],
		     paste(args[8], "_mcmc_diags_cv_", k, ".csv", sep = ""),
		     sep = "/")
      )
    # NUTS diagnostics
    fit %>%
      bayesplot::nuts_params() %>%
      saveRDS(paste(args[7],
		    paste(args[8], "_nuts_diags_cv_", k, ".rds", sep = ""),
		    sep = "/")
      )

    # log likelihood
    num_test <- nrow(test_mesh_matrix)
    log_lik_test <- extract(fit,
			    pars = paste("log_lik_test[", 1:num_test, "]", sep = "")) %>%
      unlist() %>%
      matrix(ncol = num_test) %>%
      t() %>%
      cbind(test_mesh_matrix[, "unique_id"])
    
    # posterior predictive samples 
    post_pred_test <- extract(fit,
		              pars = paste("sim_report_test[", 1:num_test, "]", sep = "")) %>%
      unlist() %>%
      matrix(ncol = num_test) %>%
      t() %>%
      cbind(test_mesh_matrix[, "unique_id"])

    # probability of report
    prob_report_test <- extract(fit,
		                pars = paste("prob_report_test[", 1:num_test, "]", sep = "")) %>%
      unlist() %>%
      matrix(ncol = num_test) %>%
      t() %>%
      cbind(test_mesh_matrix[, "unique_id"])

    list("log_lik_test" = log_lik_test, "post_pred_test" = post_pred_test, "prob_report_test" = prob_report_test)
  }

  # making larger matrices
  print("Making Larger Matrices")
  log_lik <- purrr::map(list_of_output, function(x){purrr::pluck(x, "log_lik_test")}) %>% 
    do.call(rbind, .)
  post_pred <- purrr::map(list_of_output, function(x){purrr::pluck(x, "post_pred_test")}) %>% 
    do.call(rbind, .)
  prob_report <- purrr::map(list_of_output, function(x){purrr::pluck(x, "prob_report_test")}) %>% 
    do.call(rbind, .)
  
  # converting to tibble
  print("Converting to tibble")
  log_lik_tibble <- log_lik %>% 
    as_tibble(.name_repair = "universal") %>% 
    clean_names() %>%
    # name to log_lik_1, log_lik_2
    rename_with(~ gsub("x", "log_lik_", .x)) %>% 
    # last one should be id
    rename(unique_id = paste("log_lik", num_samples + 1, sep = "_"))
  post_pred_tibble <- post_pred %>% 
    as_tibble(.name_repair = "universal") %>% 
    clean_names() %>%
    # name to post_pred_1, post_pred_2
    rename_with(~ gsub("x", "post_pred_", .x)) %>% 
    # last one should be id
    rename(unique_id = paste("post_pred", num_samples + 1, sep = "_"))
  prob_report_tibble <- prob_report %>% 
    as_tibble(.name_repair = "universal") %>% 
    clean_names() %>%
    # name to prob_report_1, prob_report_2
    rename_with(~ gsub("x", "prob_report_", .x)) %>% 
    # last one should be id
    rename(unique_id = paste("prob_report", num_samples + 1, sep = "_"))

  # binding to data
  print("Binding to data")
  log_lik_res <- prepped_model_data %>% 
    full_join(log_lik_tibble, by = "unique_id") %>%
    select(-unique_id)
  post_pred_res <- prepped_model_data %>% 
    full_join(post_pred_tibble, by = "unique_id") %>%
    select(-unique_id)
  prob_report_res <- prepped_model_data %>% 
    full_join(prob_report_tibble, by = "unique_id") %>%
    select(-unique_id)

  # saving
  print("Saving")
  log_lik_res %>%
    saveRDS(paste(args[7], 
                  paste(args[8], "log_likelihood_cv.rds", sep = "_"), 
                  sep = "/")
  )
  post_pred_res %>%
    saveRDS(paste(args[7], 
                  paste(args[8], "posterior_predictive_samples_cv.rds", sep = "_"), 
                  sep = "/")
  )
  prob_report_res %>%
    saveRDS(paste(args[7], 
                  paste(args[8], "probability_report_cv.rds", sep = "_"), 
                  sep = "/")
  ) 
}

# run model on all data for parameter information and full model evaluation
print("Running model on all data")

# as lists
columns_to_select_mesh <- ! colnames(mesh_matrix) %in% c("folds", "unique_id", "report")
columns_to_select_report <- ! colnames(report_matrix) %in% c("folds", "unique_id", "report")
stan_data <- list(N_train = mesh_matrix %>% nrow(),
                  N_test = mesh_matrix %>% nrow(),
                  X_mesh_train = mesh_matrix[, columns_to_select_mesh],
                  X_report_train = report_matrix[, columns_to_select_report], 
                  X_mesh_test = mesh_matrix[, columns_to_select_mesh],
                  X_report_test = report_matrix[, columns_to_select_report],
                  # note: reports currently stored as column in the report matrix
                  report_train = report_matrix[, "report"],
                  report_test = report_matrix[, "report"],
                  beta_hail_lower_bounds = beta_hail_lower_bounds,
                  beta_report_lower_bounds = beta_report_lower_bounds)
  
# fit full model
print("Running full model")
params_to_ignore <- c("prob_report_train", 
		      "X_mesh_train_trans",
		      "X_mesh_test_trans",
		      "X_mesh_train_std_trans",
		      "X_mesh_test_std_trans",
		      "X_report_train_std_trans",
		      "X_report_test_std_trans",
		      "mesh_train_yj",
		      "mesh_test_yj", 
		      "dens_train_yj",
		      "dens_test_yj",
		      "mesh_train_yj_mean",
		      "mesh_train_yj_sd",
		      "dens_train_yj_mean",
		      "dens_train_yj_sd",
		      "mesh_train_untrans_mean",
		      "mesh_train_untrans_sd",
		      "mesh_train_yj_hail",
		      "mesh_train_yj_report",
		      "mesh_train_yj_hail2",
		      "mesh_test_yj_hail",
		      "mesh_test_yj_report",
		      "mesh_test_yj_hail2",
		      "mesh_train_yj_hail_mean",
		      "mesh_train_yj_hail_sd",
		      "mesh_train_yj_hail2_mean",
		      "mesh_train_yj_hail2_sd",
		      "mesh_train_yj_report_mean",
		      "mesh_train_yj_report_sd"
		      )
fit <- sampling(compiled_model, data = stan_data, chains = chains,
		pars = params_to_ignore, include = FALSE,
                iter = iter, cores = parallel::detectCores())
  
print("Extracting quantities full model")

num_obs <- nrow(mesh_matrix)
post_pred <- extract(fit,
                     pars = paste("sim_report_test[", 1:num_obs, "]", sep = "")) %>%
    unlist() %>%
    matrix(ncol = num_obs) %>%
    t()

prob_report <- extract(fit,
                       pars = paste("prob_report_test[", 1:num_obs, "]", sep = "")) %>%
    unlist() %>%
    matrix(ncol = num_obs) %>%
    t()

# converting to tibble
print("Converting full results to tibble")
post_pred_tibble <- post_pred %>% 
  as_tibble(.name_repair = "universal") %>% 
  clean_names() %>%
  # name to post_pred_1, post_pred_2
  rename_with(~ gsub("x", "post_pred_", .x))
prob_report_tibble <- prob_report %>% 
  as_tibble(.name_repair = "universal") %>% 
  clean_names() %>%
  # name to prob_report_1, prob_report_2
  rename_with(~ gsub("x", "prob_report_", .x))

# binding to data
print("Binding full results to data")
post_pred_res <- prepped_model_data %>% 
  cbind(post_pred_tibble)
prob_report_res <- prepped_model_data %>% 
  cbind(prob_report_tibble)

# dealing with log likelihood - needs to be in this format for LOO
log_lik_res <- extract_log_lik(fit, parameter_name = "log_lik_test", merge_chains = FALSE)

# doing bridgesampling for marginal likelihood/Bayes factors while we have the model object
# need named bounds object
lb <- c(beta_hail_lower_bounds, beta_report_lower_bounds)
names(lb) <- c(paste("beta_hail[", 1:length(beta_hail_lower_bounds), "]", sep = ""),
	       paste("beta_report[", 1:length(beta_report_lower_bounds), "]", sep = ""))
# now doing samples
bridge_samples <- bridge_sampler(fit, lb = lb)

# saving
print("Saving full results")
# mcmc diagnostics
fit %>%
  posterior::summarise_draws() %>%
  write_csv(paste(args[7],
	    paste(args[8], "mcmc_diags_full.csv", sep = "_"),
	    sep = "/")
  )
# NUTS diagnostics
fit %>%
  bayesplot::nuts_params() %>%
  saveRDS(paste(args[7],
		paste(args[8], "nuts_diags_full.rds", sep = "_"),
		sep = "/")
  )
# previously calculated quantities
bridge_samples %>%
  saveRDS(paste(args[7], 
                paste(args[8], "bridge_samples_full.rds", sep = "_"), 
                sep = "/")
          )
log_lik_res %>%
  saveRDS(paste(args[7], 
                paste(args[8], "log_likelihood_full.rds", sep = "_"), 
                sep = "/")
          )
post_pred_res %>%
  saveRDS(paste(args[7], 
                paste(args[8], "posterior_predictive_samples_full.rds", sep = "_"), 
                sep = "/")
          )
prob_report_res %>%
  saveRDS(paste(args[7], 
                paste(args[8], "probability_report_full.rds", sep = "_"), 
                sep = "/")
          )
# grab model parameters
# setting it up in this way enables us to get both the raw and the transformed
# parameters from the standardised models 
non_model_params <- c("prob_report_train", 
		      "log_lik_test", 
		      "prob_report_test", 
		      "sim_report_test", 
		      "X_mesh_train_trans",
		      "X_mesh_test_trans",
		      "X_mesh_train_std_trans",
		      "X_mesh_test_std_trans",
		      "X_report_train_std_trans",
		      "X_report_test_std_trans",
		      "mesh_train_yj",
		      "mesh_test_yj", 
		      "dens_train_yj",
		      "dens_test_yj",
		      "mesh_train_yj_mean",
		      "mesh_train_yj_sd",
		      "dens_train_yj_mean",
		      "dens_train_yj_sd",
		      "mesh_train_untrans_mean",
		      "mesh_train_untrans_sd",
		      "mesh_train_yj_hail",
		      "mesh_train_yj_report",
		      "mesh_train_yj_hail2",
		      "mesh_test_yj_hail",
		      "mesh_test_yj_report",
		      "mesh_test_yj_hail2",
		      "mesh_train_yj_hail_mean",
		      "mesh_train_yj_hail_sd",
		      "mesh_train_yj_report_mean",
		      "mesh_train_yj_report_sd",
		      "mesh_train_yj_hail2_mean",
		      "mesh_train_yj_hail2_sd",
		      "lp__")
model_params <- fit %>%
  extract(permuted = FALSE, pars = non_model_params, include = FALSE) %>%
  posterior::as_draws_array()
model_params %>%
  saveRDS(paste(args[7],
		paste(args[8], "model_params_full.rds", sep = "_"),
		sep = "/")
  )

