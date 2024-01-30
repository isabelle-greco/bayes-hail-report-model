### Script to run the bulk simulations on LHHR/HHLR data and given experiments
### in /g/data.
### 
### Last modified: 2023-08-23
### Author: Isabelle Greco

# parse arguments
args = commandArgs(trailingOnly = TRUE)
# assert correct number of args
if (length(args) != 6) {
  stop(paste("Must give library for R package installs, directory in which to",
             "save the data, path to the original data, path to file containing",
	     "names of the models, base directory in which to find stan files,",
	     "and number of cores available."),
       call. = FALSE)
}

# load/install relevant packages
for (package in c("tidyverse", "janitor", "rstan", "foreach", "doParallel", "twinning")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1],
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# reading in the data
# args[3] = path to hte original data 
model_data <- read_csv(args[3], col_types = "nnTnnnnnnncnifffff")

# inverse logit function
inv_logit <- function(x){
  # The inverse logit or sigmoid function, used often in calculating probabilities
  #
  # Inputs:
  #  x  numeric, can be a vector
  #
  # Outputs:
  #  _  same shape as input, will be in (0, 1)

  1 / (1 + exp(-x))
}

# functions to generate the data
func_generate_lhhr <- function(data) {
  # Generate data in the low probability of hail, high probability of reporting
  # paradigm. Parameters are chosen to yield an expected number of reports events
  # close to the number of reports observed.
  #
  # Inputs:
  #  data  tibble
  #	   must have columns 'mesh', 'pop_dens'
  # 
  # Outputs:
  #  _     tibble
  #        same columns as input tibble along with integer 'report' column

  data %>%
    mutate(prob_hail = inv_logit(-7 + 0.0932 * mesh),
           prob_report_given_hail = inv_logit(-2 + 1 * log(pop_dens))) %>%
    mutate(prob_report = prob_hail * prob_report_given_hail) %>%
    add_column(rand = runif(nrow(.))) %>%
    mutate(report = case_when(rand < prob_report ~ 1,
                              TRUE ~ 0)) %>%
    select(-prob_hail, -prob_report_given_hail, -prob_report, -rand)
}

func_generate_hhlr <- function(data) {
  # Generate data in the high probability of hail, low probability of reporting
  # paradigm. Parameters are chosen to yield an expected number of reports events
  # close to the number of reports observed.
  #
  # Inputs:
  #  data  tibble
  #        must have columns 'mesh', 'pop_dens'
  # 
  # Outputs:
  #  _     tibble
  #        same columns as input tibble along with integer 'report' column

  data %>%
    mutate(prob_hail = inv_logit(-8 + 0.3 * mesh),
           prob_report_given_hail = inv_logit(-4.258 + 0.5 * log(pop_dens))) %>%
    mutate(prob_report = prob_hail * prob_report_given_hail) %>%
    add_column(rand = runif(nrow(.))) %>%
    mutate(report = case_when(rand < prob_report ~ 1,
                              TRUE ~ 0)) %>%
    select(-prob_hail, -prob_report_given_hail, -prob_report, -rand)
}


# function to split data - same for lhhr and hhlr 
func_generate_folds <- function(data, num_folds = 4, strategy = 2) {
  # Use the multiplet algorithm (Vakayil & Joseph 2022) to calculate optimal
  # cross validation folds for the input data
  #
  # Inputs:
  #  data       tibble
  #             must have columns `mesh`, `pop_dens`, and `report`
  #  num_folds  integer, default 4
  #             number of folds into which to split the data
  #             to use the default implementation, must be a power of two
  #  strategy   integer, default 2
  #             strategy to use when applying multiplet. see `twinning` 
  #             documentation for more detail
  #
  # Outputs:
  #  _          vector of integers
  #             corresponds to folds into which to split the data

  # selecting relevant columns to split into cv folds
  data_to_split <- data %>%
    mutate(report = as_factor(report)) %>% # multiplet takes care of the encoding
    select(mesh, pop_dens, report) # these simulations need only mesh, pop_dens, report

  # calculating fold labels
  multiplet(as.data.frame(data_to_split), k = num_folds, strategy = strategy)
}


# function data -> matrices, same for all data simulations
func_generate_matrices <- function(data) {
  # Function to take a dataframe and create the matrices necessary to ultimately
  # pass along to stan
  #
  # Inputs:
  #  data  tibble
  #        must have columns `report`, `mesh`, `pop_dens`, and `folds`
  #
  # Outputs:
  #  _     list
  #        named with `matrix_hail`, `matrix_report`, `beta_hail_lower_bounds` 
  #        and `beta_report_lower_bounds`

  # matrix for hail
  model_frame_hail <- data %>%
    model.frame(report ~ mesh + folds, data = .)
  model_matrix_hail <- model.matrix(report ~ mesh + folds, data = model_frame_hail)
  
  # matrix for report
  model_frame_report <- data %>%
    model.frame(~ log(pop_dens) + report + folds, data = .)
  model_matrix_report <- model.matrix(~ log(pop_dens) + report + folds, data = model_frame_report)
  
  # also creating lower bounds
  beta_hail_lower_bounds <- c(-Inf, 0.0)
  beta_report_lower_bounds <- c(-Inf, 0.0)

  # returning values
  list(matrix_hail = model_matrix_hail,
       matrix_report = model_matrix_report,
       beta_hail_lower_bounds = beta_hail_lower_bounds,
       beta_report_lower_bounds = beta_report_lower_bounds)
}

func_generate_matrix_folds <- function(matrices, k = 1) {
  # Function to take matrices (as generated by func_generate_matrices) and convert
  # to a form useful for stan. This function is designed to set this up for k-fold
  # cross-validation.
  #
  # Inputs:
  #  matrices  list
  #            named list with matrices `matrix_hail` and `matrix_report` and vectors
  #            `beta_hail_lower_bounds` and `beta_hail_upper_bounds`
  #  k         integer, default 1
  #            which fold to use. as they are all generated the same way, the results
  #            should be agnostic to this choice so we use the first
  #
  # Outputs:
  #  _         list
  #            designed to be passed to stan with elements `N_train`, `N_test`, 
  #            `X_mesh_train`, `X_report_train, `X_mesh_test`, `X_report_test`,
  #            `report_train`, `report_test`, `beta_hail_lower_bounds`, and
  #            `beta_hail_upper_bounds`

  # separate data
  train_mesh_matrix <- matrices$matrix_hail[matrices$matrix_hail[, "folds"] != k, ]
  train_report_matrix <- matrices$matrix_report[matrices$matrix_report[, "folds"] != k, ]
  test_mesh_matrix <- matrices$matrix_hail[matrices$matrix_hail[, "folds"] == k, ]
  test_report_matrix <- matrices$matrix_report[matrices$matrix_report[, "folds"] == k, ]
  
  # columns not to include in the data
  columns_to_select_mesh <- ! colnames(train_mesh_matrix) %in% c("folds", "report")
  columns_to_select_report <- ! colnames(train_report_matrix) %in% c("folds", "report")
  
  # output the list for stan
  list(N_train = train_mesh_matrix %>% nrow(),
       N_test = test_mesh_matrix %>% nrow(),
       X_mesh_train = train_mesh_matrix[, columns_to_select_mesh],
       X_report_train = train_report_matrix[, columns_to_select_report],
       X_mesh_test = test_mesh_matrix[, columns_to_select_mesh],
       X_report_test = test_report_matrix[, columns_to_select_report],
       # note: reports currently stored as column in the report matrix
       report_train = train_report_matrix[, "report"],
       report_test = test_report_matrix[, "report"],
       beta_hail_lower_bounds = matrices$beta_hail_lower_bounds,
       beta_report_lower_bounds = matrices$beta_report_lower_bounds)
}

func_generate_matrix_all <- function(matrices) {
  # Function to take matrices (as generated by func_generate_matrices) and convert
  # to a form useful for stan. This function is designed to use all the available
  # data to train the model.
  #
  # Inputs:
  #  matrices  list
  #            named list with matrices `matrix_hail` and `matrix_report` and vectors
  #            `beta_hail_lower_bounds` and `beta_hail_upper_bounds`
  #            should be agnostic to this choice so we use the first
  #
  # Outputs:
  #  _         list
  #            designed to be passed to stan with elements `N_train`, `N_test`, 
  #            `X_mesh_train`, `X_report_train, `X_mesh_test`, `X_report_test`,
  #            `report_train`, `report_test`, `beta_hail_lower_bounds`, and
  #            `beta_hail_upper_bounds`
  
  # columns not to incldue in the data
  columns_to_select_mesh <- ! colnames(matrices$matrix_hail) %in% c("folds", "report")
  columns_to_select_report <- ! colnames(matrices$matrix_report) %in% c("folds", "report")
  
  # output the list for stan
  list(N_train = matrices$matrix_hail %>% nrow(),
       N_test = matrices$matrix_hail %>% nrow(),
       X_mesh_train = matrices$matrix_hail[, columns_to_select_mesh],
       X_report_train = matrices$matrix_report[, columns_to_select_report],
       X_mesh_test = matrices$matrix_hail[, columns_to_select_mesh],
       X_report_test = matrices$matrix_report[, columns_to_select_report],
       # note: reports currently stored as column in the report matrix
       report_train = matrices$matrix_report[, "report"],
       report_test = matrices$matrix_report[, "report"],
       beta_hail_lower_bounds = matrices$beta_hail_lower_bounds,
       beta_report_lower_bounds = matrices$beta_report_lower_bounds)
}

# utility quantile function 
function_quantiles_from_draws <- function(x_matrix, par_matrix, name){
  # Function which, given a matrix at which to evaluate the probability,
  # will employ the given simulation draws to calculate quantiles of said
  # probability at each evaluation point. 
  # 
  # Inputs:
  #  x_matrix    matrix (neval x npars)
  #              each row representats a point at which to evaluate
  #              the probability and calculate its statistics
  #  par_matrix  matrix (ndraws x npars)
  #              contains draws of the parameters from stan
  #  name        str
  #              used to name the columns of the resulting tibble
  #
  # Outputs:
  #  _           tibble (neval x 8)
  #              each row corresponds to a row in x_matrix with columns
  #              corresponding to statistics (e.g. <name>_2_5_percent, ...)
  
  # probability at each x location (row) for the number of simulations (col)
  mat_func <- inv_logit(x_matrix %*% t(par_matrix))
  
  apply(mat_func, MARGIN = 1, FUN = quantile, 
	probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)) %>%
    t() %>%
    cbind(mean = rowMeans(mat_func)) %>% 
    as_tibble(.name_repair = "unique") %>% 
    # rename: <NAME>_2_5_percent, <NAME>_5_percent, etc.
    rename_with(~ paste(name, gsub("x", "", make_clean_names(.x)), sep = "_"))
}

n_func_eval <- 100 # number of points at which to evaluat e the function
# defining x_mesh vector and corresponding matrix using min and max in data
x_mesh <- seq(min(model_data$mesh), max(model_data$mesh), length.out = n_func_eval)
x_mesh_matrix <- matrix(c(rep(1, n_func_eval), x_mesh), ncol = 2)
# defining x_pop_dens and corresponding matrix again with min and max in data
# note the logarithm (hence the name it is saved at) 
x_pop_dens <- seq(min(log(model_data$pop_dens)), max(log(model_data$pop_dens)),
                  length.out = n_func_eval)
x_pop_dens_matrix <- matrix(c(rep(1, n_func_eval), x_pop_dens), ncol = 2)
# saving to use later for plotting
# note the saving highlights the logarithm
saveRDS(x_mesh_matrix, paste(args[2], "/x_mesh_matrix.rds", sep = ""))
saveRDS(x_pop_dens_matrix, paste(args[2], "/x_log_pop_dens_matrix.rds", sep = ""))

func_process_model_output <- function(stanfit, name, x_mesh_matrix_arg = x_mesh_matrix, 
				      x_pop_dens_matrix_arg = x_pop_dens_matrix) {
  # Function to process stan model output according to the matrices just defined
  #
  # Inputs:
  #  stanfit                stanfit
  #                         stan model output, must contain parameters `beta_hail`
  #                         and `beta_report` of size npars1 and npars2 respectively
  #  name                   string
  #                         name utilised to label the resulting output tibble's
  #                         columns
  #  x_mesh_matrix_arg      matrix (neval x npars1)
  #                         matrix at which to evaluate the probability of hail
  #                         given MESH and consider its statistics
  #  x_pop_dens_matrix_arg  matrix (neval x npars2)
  #                         matrix at which to evaluate the conditional probability
  #                         of reporting and consider its statistics
  #
  # Outputs:
  #  _                      tibble
  #                         columns will be of form <name>_hail_2_5_percent etc. and
  #                         <name>_report_2_5_percent etc.

  # hail curve
  hail <- function_quantiles_from_draws(
    x_mesh_matrix_arg,
    rstan::extract(stanfit, pars = "beta_hail")$beta_hail,
    paste(name, "hail", sep = "_")
    )
  # reporting curve
  report <- function_quantiles_from_draws(
    x_pop_dens_matrix_arg,
    rstan::extract(stanfit, pars = "beta_report")$beta_report,
    paste(name, "report", sep = "_")
  )
  # bind to return
  cbind(hail, report) %>% # returns matrix
    as_tibble()
}


# vector of model names from the file given to args[4]
model_names <- read_lines(args[4]) 
# base dir to find model script from args[5]
base_dir <- args[5] 

# compiling lhhr models by mapping over vertor to create list
# could parallelise with furrr https://furrr.futureverse.org/
lhhr_models <- map(model_names[str_detect(model_names, "lhhr")], 
		   ~ stan_model(file = paste(base_dir, .x, paste(.x, "stan", sep = "."), sep = "/")))
# labelling the list 
names(lhhr_models) <- model_names[str_detect(model_names, "lhhr")]

# compiling hhlr models by mapping over vector to create list
hhlr_models <- map(model_names[str_detect(model_names, "hhlr")], 
		   ~ stan_model(file = paste(base_dir, .x, paste(.x, "stan", sep = "."), sep = "/")))
# labelling the list
names(hhlr_models) <- model_names[str_detect(model_names, "hhlr")]


# function to pull it all togetehr and be run in parallel
func_one_simulation_iteration <- function(data, sim_number) {
  # Function which, given data and the simulation id, generate s the LHHR and HHLR
  # data, creates cross validation folds, runs the models on all but one of the folds
  # and on all the data, evaluates statistics of the probabilities as relevant points
  # and saves the combined results.
  # 
  # Inputs:
  #  data        tibble
  #              will be passed along to all the previously defined functions and so
  #              must contain `mesh` and `pop_dens` columns
  #  sim_number  int
  #              identifies the simulation when it is saved enabling tracking of 
  #              the overall progress
  # 
  # Outputs:
  #  NULL        NULL
  #              final step of this function is saving a result which returns
  #              the NULL value invisibly
  #		 saved object is list with names `lhhr_res_cv`, `hhlr_res_cv`,
  #              `lhhr_res_all`, and `hhlr_res_all`
  #              saved in args[2] directory was bulk_sims_res_<sim_number>.rds

  # setting up rstan
  # done here so that it is relevant to the worker running the function	
  options(mc.cores = parallel::detectCores()) 

  # generate data
  lhhr <- func_generate_lhhr(data) %>%
    mutate(folds = func_generate_folds(.)) 
  hhlr <- func_generate_hhlr(data) %>%
    mutate(folds = func_generate_folds(.))

  # generate matrices
  lhhr_matrices <- func_generate_matrices(lhhr)
  hhlr_matrices <- func_generate_matrices(hhlr)

  # generate stan list - all
  lhhr_stan_list_folds <- func_generate_matrix_folds(lhhr_matrices)
  hhlr_stan_list_folds <- func_generate_matrix_folds(hhlr_matrices)

  # parameters to save in each model run - ignoring the rest saves space
  # abd tune
  pars <- c("beta_hail", "beta_report")
  # run models - folds
  lhhr_res_cv <- foreach (model = lhhr_models, name = names(lhhr_models), 
			  .final = function(x) setNames(x, paste(names(lhhr_models), "cv", sep = "_"))) %do% {
    sampling(model, data = lhhr_stan_list_folds, pars = pars) %>%
      func_process_model_output(name = paste(name, "cv", sep = "_"))
  }

  hhlr_res_cv <- foreach (model = hhlr_models, name = names(hhlr_models), 
			  .final = function(x) setNames(x, paste(names(hhlr_models), "cv", sep = "_"))) %do% {
    sampling(model, data = hhlr_stan_list_folds, pars = pars) %>%
      func_process_model_output(name = paste(name, "cv", sep = "_"))
  }

  # generate stan list - all
  lhhr_stan_list_all <- func_generate_matrix_all(lhhr_matrices)
  hhlr_stan_list_all <- func_generate_matrix_all(hhlr_matrices)

  # run models - all
  lhhr_res_all <- foreach (model = lhhr_models, name = names(lhhr_models), 
			   .final = function(x) setNames(x, paste(names(lhhr_models), "all", sep = "_"))) %do% {
    sampling(model, data = lhhr_stan_list_all, pars = pars) %>%
      func_process_model_output(name = paste(name, "all", sep = "_"))
  }

  hhlr_res_all <- foreach (model = hhlr_models, name = names(hhlr_models), 
			   .final = function(x) setNames(x, paste(names(hhlr_models), "all", sep = "_"))) %do% {
    sampling(model, data = hhlr_stan_list_all, pars = pars) %>%
      func_process_model_output(name = paste(name, "all", sep = "_"))
  }

  # save results - saving separately and reading into list later
  list(lhhr_res_cv = lhhr_res_cv,
       hhlr_res_cv = hhlr_res_cv,
       lhhr_res_all = lhhr_res_all,
       hhlr_res_all = hhlr_res_all) %>%
    saveRDS(paste(args[2], "/bulk_sims_res_", sim_number, ".rds", sep = ""))
}

# set up parallel process
# this set up gives four cores (for foru chains) per work with some left
# to spare for the main process
cl <- makeCluster(floor((as.numeric(args[6]) - 1) / 4))
registerDoParallel(cl, cores = parallel::detectCores() - 4)

# running all the sims (returns list of null as all is saved)
foreach(i = 1:100, .packages = c("tidyverse", "foreach", "rstan", "janitor", "twinning")) %dopar% {
 func_one_simulation_iteration(model_data, i)
}

# stopping cluster
stopCluster(cl)
