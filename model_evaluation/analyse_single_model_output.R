### Script to evaluate a single model output. Produces analysis/plots for the MCMC
### diagnostics as well as model fit
###
### Last Modified: 2023-08-23
### Author: Isabelle Greco


# get args from command line
# 1 = library for r installs if need be
# 2 = path to entire data set 
# 3 = model name
# 4 = path to model results
# 5 = directory in which to save analyses
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 5) {
  stop(paste("Must give library for R package installs, location of full data",
             "name of the model to analyse, path to model results and directory", 
             "in which to save analyses."), 
       call. = FALSE)
}


# load libraries + install if need be
for (package in c("tidyverse", "posterior", "bayesplot", "janitor", "latex2exp", 
                  "kableExtra", "paletteer", "poibin", "patchwork")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# options for plot size
plot_width = 15
plot_height = 10


# reading in model data (may be sim) from args[2]
model_data <- read_csv(args[2], col_types = "nnTnnnnnnncnifffff")

### mcmc diagnostics ###

# read in full MCMC diagnostics
results_folder <- paste(args[4], "results", sep = "/")
eval_folder <- paste(args[4], "eval", sep = "/")
full_mcmc_diags <- read_csv(paste(results_folder, 
                                       paste(args[3], "mcmc_diags_full.csv", sep = "_"), 
                                       sep = "/"))

# saving table of mesh diagnostics to more easily visualise
full_mcmc_diags %>%
  filter(grepl("beta", variable, fixed = TRUE) | variable == "lp__") %>%
  kbl("markdown", digits = 3) %>% # latex, booktabs = TRUE for latex code
  write_lines(paste(args[5], 
                    paste(args[3], "mcmc_diags_full_table_markdown.txt", sep = "_"), 
                    sep = "/")
             )


# appending same table for each fold to compare
append_cv_results_to_table <- function(i) {
  # Function to append the mcmc diagnostics table of each fold to the table 
  # created above from the entire data set.
  #
  # Inputs:
  #  i     integer
  #        cross validation fold identifier
  #
  # Outputs:
  #  NULL  NULL
  #        table is written to file
   
  # read in csv of the results table  
  read_csv(paste(results_folder, 
                 paste0(args[3], "_mcmc_diags_cv_", i, ".csv"), 
                 sep = "/")) %>%
    # filter only key parameters (beta) and the log posterior
    filter(grepl("beta", variable, fixed = TRUE) | variable == "lp__") %>%
    # kable to make markdown table
    kbl("markdown", digits = 3) %>% # latex, booktabs = TRUE for latex code
    # appendingto previous file
    write_lines(paste(args[5], 
                      paste(args[3], "mcmc_diags_full_table_markdown.txt", sep = "_"), 
                      sep = "/"), 
                append = TRUE)
}
# running this function over the four folds
walk(1:4, append_cv_results_to_table)

# rhat histograms
make_rhat_histogram <- function(mcmc_diags, var_name, plot_title){
  # Function to create a histogrma of Rhat results for a given variable
  #
  # Inputs:
  #  mcmc_diags  tibble
  #              must contain a column called `rhat` with the rhat values (typically
  #              created by posterior)
  #  var_name    str
  #              name of the variable (e.g. prob_report, etc.)
  #  plot_title  str
  #              title for the plot
  #
  # Outputs:
  #  _
    
  mcmc_diags %>%
    # filtering for the variable name
    filter(grepl(var_name, variable, fixed = TRUE)) %>%
    # extracting rhat
    pull(rhat) %>%
    # plotting
    mcmc_rhat_hist() +
    # adding title
    ggtitle(plot_title)
}

save_plot <- function(p, name, width = plot_width, height = plot_height){
  # Generate function to save plots in the right location and right size as a png
  # 
  # Inputs:
  #  p       ggplot object
  #          the plot to be saved
  #  name    str
  #          name of the variable(s)
  #  width   numeric
  #          width of the plot
  #          default given by variable `plot_width`
  #  height  numeric
  #          height of the plot
  #          default given by variable `plot_height`
  # 
  # Outputs:
  #  NULL    NULL  
  #          plot is written to file in a location defined by args[5] with name
  #          defined by args[3]
    
  ggsave(filename = paste(args[5], paste(args[3], name, sep = "_"), sep = "/"), 
         plot = p, 
         device = png, 
         width = width, 
         height = height)
}

# making rhat histograms using previous function
# ...for probability of report
make_rhat_histogram(full_mcmc_diags, "prob_report", 
                    TeX(r'($\widehat{R}$ for the probability of a report in each grid cell)')) %>%
  save_plot("mcmc_prob_report_full_rhat.png")

# ...for the log likelihood
make_rhat_histogram(full_mcmc_diags, "log_lik", 
                    TeX(r'($\widehat{R}$ for the log likelihood of a report in each grid cell)')) %>%
  save_plot("mcmc_log_lik_full_rhat.png")

# neff  histograms
make_bulk_ess_histogram <- function(mcmc_diags, var_name, plot_title){
  # Function to create a histogram of bulk ESS ratio for a given variable
  #
  # Inputs:
  #  mcmc_diags  tibble
  #              must contain a column called `ess_bulk` with the bulk ESS values 
  #              (typically created by posterior)
  #  var_name    str
  #              name of the variable (e.g. prob_report, etc.)
  #  plot_title  str
  #              title for the plot
  #
  # Outputs:
  #  _
    
  mcmc_diags %>%
    # filtering for correct variable
    filter(grepl(var_name, variable, fixed = TRUE)) %>%
    # calculating ratio - assumption of 4 chains with 1000 iterations
    mutate(neff_bulk_ratio = ess_bulk / (4 * 1000)) %>%
    pull(neff_bulk_ratio) %>%
    # creating histogram
    mcmc_neff_hist() +
    # adding title
    ggtitle(plot_title)   
}

make_tail_ess_histogram <- function(mcmc_diags, var_name, plot_title){
  # Function to create a histogram of tail ESS ratio for a given variable
  #
  # Inputs:
  #  mcmc_diags  tibble
  #              must contain a column called `ess_tail` with the tail ESS values 
  #              (typically created by posterior)
  #  var_name    str
  #              name of the variable (e.g. prob_report, etc.)
  #  plot_title  str
  #              title for the plot
  #
  # Outputs:
  #  _
    
  mcmc_diags %>%
    # filtering for correct variable
    filter(grepl(var_name, variable, fixed = TRUE)) %>%
    # calculating ratio under the assumption of 4 chains with 1000 observations
    mutate(neff_tail_ratio = ess_tail / (4 * 1000)) %>%
    pull(neff_tail_ratio) %>%
    # creating histogram
    mcmc_neff_hist() +
    # adding title
    ggtitle(plot_title)   
}

# making ess histograms for prob report
make_bulk_ess_histogram(full_mcmc_diags, "prob_report",
                        "Bulk effective sample size ratio for the probability of a report in each grid cell") %>%
  save_plot("mcmc_prob_report_full_ess_bulk.png")

make_tail_ess_histogram(full_mcmc_diags, "prob_report", 
                        "Tail effective sample size ratio for the probability of a report in each grid cell") %>%
  save_plot("mcmc_prob_report_full_ess_tail.png")

# making ess histograms for log likelihood
make_bulk_ess_histogram(full_mcmc_diags, "log_lik", 
                        "Bulk effective sample size ratio for the log-likelihood of a report in each cell") %>%
  save_plot("mcmc_log_lik_full_ess_bulk.png")

make_tail_ess_histogram(full_mcmc_diags, "log_lik", 
                        "Tail effective sample size ratio for the log-likelihood of a report in each cell") %>%
  save_plot("mcmc_log_lik_full_ess_tail.png")

# reading in full parameters 
full_params <- readRDS(paste(results_folder, 
                             paste(args[3], "model_params_full.rds", sep = "_"), 
                             sep = "/"))

# acf plot
full_params %>%
  mcmc_acf() %>%
  save_plot("mcmc_full_params_acf.png")

# rank hist
full_params %>%
  mcmc_rank_hist(ref_line = TRUE) %>%
  save_plot("mcmc_full_params_rank_hist.png")

# ecdf
full_params %>%
  mcmc_rank_ecdf(prob = 0.99, plot_diff = TRUE) %>%
  save_plot("mcmc_full_params_rank_ecdf.png")


### model fit evaluation ###

# diagnostics from the NUTS sampler
full_nuts_diags <- readRDS(paste(results_folder, 
				  paste(args[3], "nuts_diags_full.rds", sep = "_"), 
                           sep = "/"))

# pairs plot of the beta parameters
full_params %>%
  mcmc_pairs(np = full_nuts_diags) %>%
  save_plot("fit_key_params.png")

# get names of non standardised parameters
# useful for plotting curves and comparing with true values 
non_standardised_params <- dimnames(full_params)$variable[
    str_detect(dimnames(full_params)$variable, "raw", negate = TRUE)
]

# splitting up the model name...
split_model_name <- args[3] %>%
  str_split_1("_")
# ...to check if its a sim run or not 
sim <- "sim" %in% split_model_name

# defining the tru params in case needed
params_list <- list("lhhr_beta_hail" = c(-7, 0.0932),
                    "lhhr_beta_report" = c(-2, 1),
                    "hhlr_beta_hail" = c(-8, 0.3),
                    "hhlr_beta_report" = c(-4.258, 0.5))
# selecting the relevant params if it is a sim run
if (sim) {
  data_set <- if ("lhhr" %in% split_model_name) "lhhr" else "hhlr"
  beta_hail_true <- params_list[[paste(data_set, "beta_hail", sep = "_")]]
  beta_report_true <- params_list[[paste(data_set, "beta_report", sep = "_")]]
}

# recover hist plot - for simulations only
if (sim) {
  # full_params is a draws_array with dimensions iterations x chain x variables
  # selecting only the non-standardised parameters for comparison
  full_params[, , non_standardised_params] %>%
    mcmc_recover_hist(true = c(beta_hail_true, beta_report_true)) %>%
    # saving plot
    save_plot("fit_recover_hist_full.png")
}

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

# only needed for we have a transformed mesh
yeo_johnson <- function(x, lambda) {
  # Yeo-Johnson transform of x using parameter lambda
  #
  # Inputs:
  #  x       numeric
  #          the number (not vector) to be transformed
  #  lambda  numeric
  #          real parameter defining the transform
  # 
  # Outputs:
  #  _       numeric
  #          transformed vector same shape as x
    
  eps <- 1e-5
  if (x >= 0.0){
    if (abs(lambda) < eps){
      log1p(x)
    } else {
      (((x + 1.0) ^ lambda) - 1.0) / lambda
    }
  } else {
    if (abs(lambda - 2.0) < eps) {
      -log1p(-x)
    } else {
      -(((-x + 1.0) ^ (2.0 - lambda)) - 1.0) / (2.0 - lambda)
    }
  }
}

# utility quantile function 
function_quantiles_from_draws <- function(x_matrix, par_matrix, name, hail_trans = FALSE, hail_and_mesh_trans = FALSE,
					  report_trans = FALSE, report_and_mesh_trans = FALSE, hail_sqrt = FALSE, report_nolog = FALSE,
					  report_untrans_mesh_trans = FALSE){
  # Function which, given a matrix at which to evaluate the probability,
  # will employ the given simulation draws to calculate quantiles of said
  # probability at each evaluation point. 
  # 
  # Inputs:
  #  x_matrix      matrix (neval x npars)
  #                each row representats a point at which to evaluate
  #                the probability and calculate its statistics
  #  par_matrix    matrix (ndraws x npars)
  #                contains draws of the parameters from stan
  #  name          str
  #                used to name the columns of the resulting tibble
  #  hail_trans    logical
  #                indicates whether this function is being used for a model where
  #                one of the elements of beta was use to transform mesh
  #  hail_and_mesh_trans logical
  #                      indicates whether this is a hail model in which the primary
  #                      and secondary mesh values are both transformed
  #  report_trans  logical
  #                indicates whether this function is being used for a model where
  #                one of the elements in beta transformed (raw) population density
  #  report_and_mesh_trans  logical
  #                         indicates whether this function is being using for a model
  #                         where both dens and mesh are transformed
  #  hail_sqrt     logical
  #                indicates whether this function is being used for a model where
  #                the sqrt(MESH) is used
  #  report_nolog  logical
  #                indicates whether this function is being used for a model where
  #                raw population density is used
  #  report_untrans_mesh_trans  logical
  #                             indicates whether in the reporting model mesh is
  #                             transformed but report is not
  #
  # Outputs:
  #  _             tibble (neval x 8)
  #                each row corresponds to a row in x_matrix with columns
  #                corresponding to statistics (e.g. <name>_2_5_percent, ...)
  
  # create mat_func: matrix with probability at each x location (row) for 
  # the number of simulations (col)
  
  # firstly considering transformed case
  # assumes 3rd element of beta used to transform second element of x
  # written for mesh but works equally well for population density
  if (report_and_mesh_trans) {    
    # similar to above but needs to use the raw not log pop dens
    raw_mult <- matrix(0, nrow = nrow(x_matrix), ncol = nrow(par_matrix))
    # fill matrix by column
    for (j in 1:ncol(raw_mult)) {
      # transform mesh (second column) using jth parameter in third column
      # note we are putting mesh back into original scale
      dens_trans <- map_dbl(exp(x_matrix[, 2]), \(x) yeo_johnson(x, par_matrix[j, 3]))
      mesh_trans <- map_dbl(x_matrix[, 3], \(x) yeo_johnson(x, par_matrix[j, 5]))
      # create new matrix with this second column 
      if (ncol(x_matrix) == 3) {
        dens_trans_matrix <- cbind(x_matrix[, 1], dens_trans, mesh_trans)
      } else {
        dens_trans_matrix <- cbind(x_matrix[, 1], dens_trans, mesh_trans, x_matrix[, 4:ncol(x_matrix)])
      }
      # doing the matrix multiplication to get that one column
      # should need transpose but R does it itself for a vector
      raw_mult[, j] <- dens_trans_matrix %*% par_matrix[j, c(-3, -5)]
     }
    # inverse logit transform of the matrix to get probabilities 
    mat_func <- inv_logit(raw_mult)
  } else if (hail_and_mesh_trans) {    
    # similar to above but needs to use the raw not log pop dens
    raw_mult <- matrix(0, nrow = nrow(x_matrix), ncol = nrow(par_matrix))
    # fill matrix by column
    for (j in 1:ncol(raw_mult)) {
      # transform mesh (second column) using jth parameter in third column
      # note we are putting mesh back into original scale
      mesh_trans1 <- map_dbl(x_matrix[, 2], \(x) yeo_johnson(x, par_matrix[j, 3]))
      mesh_trans2 <- map_dbl(x_matrix[, 3], \(x) yeo_johnson(x, par_matrix[j, 5]))
      # create new matrix with this second column 
      if (ncol(x_matrix) == 3) {
        x_trans_matrix <- cbind(x_matrix[, 1], mesh_trans1, mesh_trans2)
      } else {
        x_trans_matrix <- cbind(x_matrix[, 1], mesh_trans1, mesh_trans2, x_matrix[, 4:ncol(x_matrix)])
      }
      # doing the matrix multiplication to get that one column
      # should need transpose but R does it itself for a vector
      raw_mult[, j] <- x_trans_matrix %*% par_matrix[j, c(-3, -5)]
     }
    # inverse logit transform of the matrix to get probabilities 
    mat_func <- inv_logit(raw_mult)
  } else if (report_untrans_mesh_trans) {
    # pre-allocate matrix for the raw multiplication results
    raw_mult <- matrix(0, nrow = nrow(x_matrix), ncol = nrow(par_matrix))
    # fill matrix by column
    for (j in 1:ncol(raw_mult)) {
      # transform mesh (third column) using jth parameter in 4th column
      mesh_trans <- map_dbl(x_matrix[, 3], \(x) yeo_johnson(x, par_matrix[j, 4]))
      # create new matrix with this second column 
      if (ncol(x_matrix) == 3) {
        mesh_trans_matrix <- cbind(x_matrix[, 1:2], mesh_trans)
      } else {
        mesh_trans_matrix <- cbind(x_matrix[, 1:2], mesh_trans, x_matrix[, 4:ncol(x_matrix)])
      }
      # doing the matrix multiplication to get that one column
      # should need transpose but R does it itself for a vector
      raw_mult[, j] <- mesh_trans_matrix %*% par_matrix[j, -3] 
    }
      # inverse logit transform of the matrix to get probabilities 
    mat_func <- inv_logit(raw_mult)
  } else if (hail_trans) {
    # pre-allocate matrix for the raw multiplication results
    raw_mult <- matrix(0, nrow = nrow(x_matrix), ncol = nrow(par_matrix))
    # fill matrix by column
    for (j in 1:ncol(raw_mult)) {
      # transform mesh (second column) using jth parameter in third column
      mesh_trans <- map_dbl(x_matrix[, 2], \(x) yeo_johnson(x, par_matrix[j, 3]))
      # create new matrix with this second column 
      if (ncol(x_matrix) == 2) {
        mesh_trans_matrix <- cbind(x_matrix[, 1], mesh_trans)
      } else {
        mesh_trans_matrix <- cbind(x_matrix[, 1], mesh_trans, x_matrix[, 3:ncol(x_matrix)])
      }
      # doing the matrix multiplication to get that one column
      # should need transpose but R does it itself for a vector
      raw_mult[, j] <- mesh_trans_matrix %*% par_matrix[j, -3] 
    }
      # inverse logit transform of the matrix to get probabilities 
    mat_func <- inv_logit(raw_mult)
  } else if (report_trans) {    
    # similar to above but needs to use the raw not log pop dens
    raw_mult <- matrix(0, nrow = nrow(x_matrix), ncol = nrow(par_matrix))
    # fill matrix by column
    for (j in 1:ncol(raw_mult)) {
      # transform mesh (second column) using jth parameter in third column
      # note we are putting mesh back into original scale
      dens_trans <- map_dbl(exp(x_matrix[, 2]), \(x) yeo_johnson(x, par_matrix[j, 3]))
      # create new matrix with this second column 
      if (ncol(x_matrix) == 2) {
        dens_trans_matrix <- cbind(x_matrix[, 1], dens_trans)
      } else {
        dens_trans_matrix <- cbind(x_matrix[, 1], dens_trans, x_matrix[, 3:ncol(x_matrix)])
      }
      # doing the matrix multiplication to get that one column
      # should need transpose but R does it itself for a vector
      raw_mult[, j] <- dens_trans_matrix %*% par_matrix[j, -3]
     }
    # inverse logit transform of the matrix to get probabilities 
    mat_func <- inv_logit(raw_mult)
  } else if (hail_sqrt) {
    # taking square root of mesh, assumed to be in the second column
    # assigning directly as not having to loop
    x_matrix[, 2] <- sqrt(x_matrix[, 2])
    # calculating the function in the transformed matrix
    mat_func <- inv_logit(x_matrix %*% t(par_matrix))
  } else if (report_nolog) {
    # taking exponential of pop density, assumed to be in the second column
    x_matrix[, 2] <- exp(x_matrix[, 2])
    # calculing function in transformed matric
    mat_func <- inv_logit(x_matrix %*% t(par_matrix))
  } else {
    # since no transform can just multiply directly
    mat_func <- inv_logit(x_matrix %*% t(par_matrix))
  }
    
  # using apply to take quantiles over matrix
  apply(mat_func, MARGIN = 1, FUN = quantile, 
    probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)) %>%
    t() %>%
    # adding row means 
    cbind(mean = rowMeans(mat_func)) %>% 
    # converting to tibble
    as_tibble(.name_repair = "unique") %>% 
    # rename: <NAME>_2_5_percent, <NAME>_5_percent, etc.
    rename_with(~ paste(name, gsub("x", "", make_clean_names(.x)), sep = "_"))
}


n_func_eval <- 100 # number of points at which to evaluate the function

# defining x_mesh vector and corresponding matrix using min and max in data
x_mesh <- seq(min(model_data$mesh), max(model_data$mesh), length.out = n_func_eval)
x_mesh_matrix <- matrix(c(rep(1, n_func_eval), x_mesh), ncol = 2)

# defining x_pop_dens and corresponding matrix again with min and max in data
# note the logarithm (hence the name it is saved at) 
x_pop_dens <- seq(min(log(model_data$pop_dens)), max(log(model_data$pop_dens)),
                  length.out = n_func_eval)
x_pop_dens_matrix <- matrix(c(rep(1, n_func_eval), x_pop_dens), ncol = 2)

# as we include extra predicts use cbind conditional upon str_detect output
# to augment these matrices
# TODO: simplify logic
if (str_detect(args[3], "report_trans_std_dens_weekend") | str_detect(args[3], "report_std_dens_weekend")) {
  # look at the weekday version
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 0)
}

if (str_detect(args[3], "report_trans_std_dens_timeofday") | str_detect(args[3], "report_std_dens_timeofday")) {
  # look at the weekday version
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 0)
}

if (str_detect(args[3], "hail_trans_std_mesh_timeofday") | str_detect(args[3], "hail_std_sqrt_mesh_timeofday")) {
  # look at the weekday version
  x_mesh_matrix <- cbind(x_mesh_matrix, 0)
}

if (str_detect(args[3], "report_trans_std_dens_poly_month") | str_detect(args[3], "report_std_dens_poly_month")) {
  # look at the november curve
  months_unique <- c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
  contr_poly_res <- contr.poly(factor(months_unique, levels = months_unique, ordered = TRUE))
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, matrix(rep(contr_poly_res[3, ], n_func_eval), nrow = n_func_eval, byrow = TRUE))
}

if (str_detect(args[3], "report_trans_std_dens_treatment_month") | str_detect(args[3], "report_std_dens_treatment_month")) {
  # look at the november curve
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 0, 1, 0, 0, 0, 0, 0)
}

if (str_detect(args[3], "hail_trans_std_mesh_poly_month") | str_detect(args[3], "hail_std_sqrt_mesh_poly_month")) {
  # look at the november curve
  months_unique <- c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
  contr_poly_res <- contr.poly(factor(months_unique, levels = months_unique, ordered = TRUE))
  x_mesh_matrix <- cbind(x_mesh_matrix, matrix(rep(contr_poly_res[3, ], n_func_eval), nrow = n_func_eval, byrow = TRUE))
}

if (str_detect(args[3], "hail_trans_std_mesh_treatment_month") | str_detect(args[3], "hail_std_sqrt_mesh_treatment_month")) {
  # look at the november curve
  x_mesh_matrix <- cbind(x_mesh_matrix, 0, 1, 0, 0, 0, 0, 0)
}

if (str_detect(args[3], "report_trans_std_dens_poly_stormseason") | str_detect(args[3], "report_std_dens_poly_stormseason")) {
  # look at the 2013/14 curve
  months_unique <- c("2009/10", "2010/11", "2011/12", "2012/13", "2013/14", "2014/15", "2015/16")
  contr_poly_res <- contr.poly(factor(months_unique, levels = months_unique, ordered = TRUE))
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, matrix(rep(contr_poly_res[5, ], n_func_eval), nrow = n_func_eval, byrow = TRUE))
}

if (str_detect(args[3], "hail_trans_std_mesh_poly_stormseason") | str_detect(args[3], "hail_std_sqrt_mesh_poly_stormseason")) {
  # look at the 2013/14 curve
  months_unique <- c("2009/10", "2010/11", "2011/12", "2012/13", "2013/14", "2014/15", "2015/16")
  contr_poly_res <- contr.poly(factor(months_unique, levels = months_unique, ordered = TRUE))
  x_mesh_matrix <- cbind(x_mesh_matrix, matrix(rep(contr_poly_res[5, ], n_func_eval), nrow = n_func_eval, byrow = TRUE))
}

if (str_detect(args[3], "report_trans_std_dens_treatment_stormseason") | str_detect(args[3], "report_std_dens_treatment_stormseason")) {
  # look at the 2013/14 curve
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 0, 0, 0, 1, 0, 0)
}

if (str_detect(args[3], "hail_trans_std_mesh_treatment_stormseason") | str_detect(args[3], "hail_std_sqrt_mesh_treatment_stormseason")) {
  # look at the 2013/14 curve
  x_mesh_matrix <- cbind(x_mesh_matrix, 0, 0, 0, 1, 0, 0)
}

if (str_detect(args[3], "report_trans_std_dens_std_mesh") | str_detect(args[3], "report_std_dens_std_mesh")) {
  # look at the 30mm mesh curve
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 30)
}

if (str_detect(args[3], "report_trans_std_dens_std_sqrt_mesh") | str_detect(args[3], "report_std_dens_std_sqrt_mesh")) {
  # look at the 30mm mesh curve
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, sqrt(30))
}

if (str_detect(args[3], "report_trans_std_dens_trans_std_mesh") | str_detect(args[3], "report_std_dens_trans_std_mesh")) {
  # look at the 30mm mesh curve
  # will need to transform in the function
  x_pop_dens_matrix <- cbind(x_pop_dens_matrix, 30)
}

if (str_detect(args[3], "hail_std_sqrt_mesh_std_mesh") | str_detect(args[3], "hail_trans_std_mesh_std_mesh")) {
  # using just 'mesh' captures all the mesh variations 
  # look at the 30mm mesh_east curve
  x_mesh_matrix <- cbind(x_mesh_matrix, 30)
}

if (str_detect(args[3], "hail_std_sqrt_mesh_std_sqrt_mesh") | str_detect(args[3], "hail_trans_std_mesh_std_sqrt_mesh")) {
  # using just 'mesh' captures all the mesh variations 
  # look at the 30mm mesh_east curve
  x_mesh_matrix <- cbind(x_mesh_matrix, sqrt(30))
}

if (str_detect(args[3], "hail_std_sqrt_mesh_trans_std_mesh")) {
  # look at the 30mm mesh_east curve
  x_mesh_matrix <- cbind(sqrt(x_mesh_matrix), 30)
}

if (str_detect(args[3], "hail_trans_std_mesh_trans_std_mesh")) {
  # look at the 30mm mesh_... curve
  x_mesh_matrix <- cbind(x_mesh_matrix, 30)
}

# if hail transformation then will need to use parameters differently
hail_trans <- str_detect(args[3], "hail_trans_std_mesh")
hail_sqrt <- str_detect(args[3], "hail_std_sqrt_mesh")
report_nolog <- str_detect(args[3], "nolog_dens")
report_trans <- str_detect(args[3], "report_trans_std_dens")
report_and_mesh_trans <- str_detect(args[3], "report_trans_std_dens_trans_std_mesh")
hail_and_mesh_trans <- str_detect(args[3], "hail_trans_std_mesh_trans_std_mesh")
report_untrans_and_mesh_trans <- str_detect(args[3], "report_std_dens_trans_std_mesh")
hail_untrans_and_mesh_trans <- str_detect(args[3], "hail_std_sqrt_mesh_trans_std_mesh") 
# TODO: regex to make it also pick up non-standardised ones?)

# reporting curve
p <- function_quantiles_from_draws(x_pop_dens_matrix, 
                                   # only using the non-standardised parameters
                                   as_draws_df(full_params[, , non_standardised_params]) %>% 
                                     select(starts_with("beta_report")) %>% 
                                     as.matrix(),
                                   "report", # setting name to create
				   report_nolog = report_nolog,
				   report_trans = report_trans,
				   report_and_mesh_trans = report_and_mesh_trans,
				   report_untrans_mesh_trans = report_untrans_and_mesh_trans
                                    ) %>%
  # adding the vector column for plotting
  add_column(x_log_pop_dens = x_pop_dens) %>%
  ggplot(aes(x = exp(x_log_pop_dens))) + 
  # plotting posterior mean
  geom_line(aes(y = report_mean, colour = "Posterior Mean")) + 
  # plotting confidence bands
  geom_ribbon(aes(ymin = report_2_5_percent, ymax = report_97_5_percent, fill = "95%"), alpha = 0.3) + 
  geom_ribbon(aes(ymin = report_25_percent, ymax = report_75_percent, fill = "50%"), alpha = 0.5) +
  # labels
  xlab("Population Density [people / km2]") + 
  ylab("Posterior P(Report = 1| Hail = 1, Population Density)") + 
  labs(fill = "Confidence", color = "") +
  ggtitle("Reporting Posterior") +
  # log scale on x axis
  scale_x_log10() +
  # colors
  scale_fill_paletteer_d("nord::aurora")
# adding true curve if a simulation
if (sim) {
  p <- p +
    geom_line(aes(y = inv_logit(x_pop_dens_matrix %*% beta_report_true), color = "True Probability"))
}
# saving plot
save_plot(p, "fit_report_curve.png")

# hail probability curve
p <- function_quantiles_from_draws(x_mesh_matrix, 
                                   # again only need the standardised params
                                   as_draws_df(full_params[, , non_standardised_params]) %>% 
                                     select(starts_with("beta_hail")) %>% 
                                     as.matrix(),
                                   "hail",
                                   hail_trans = hail_trans, # informing if we need to deal with transform
				   hail_sqrt = hail_sqrt, # informing if we need to deal with sqrt transform
				   report_untrans_mesh_trans = hail_untrans_and_mesh_trans,
				   hail_and_mesh_trans = hail_and_mesh_trans
                                   ) %>%
  # adding mesh column for plotting
  add_column(x_mesh = x_mesh) %>%
  ggplot(aes(x = x_mesh)) + 
  # plotting posterior mean
  geom_line(aes(y = hail_mean, colour = "Posterior Mean")) + 
  # plotting confidence bands
  geom_ribbon(aes(ymin = hail_2_5_percent, ymax = hail_97_5_percent, fill = "95%"), alpha = 0.3) + 
  geom_ribbon(aes(ymin = hail_25_percent, ymax = hail_75_percent, fill = "50%"), alpha = 0.5) +
  # labelling
  xlab("MESH [mm]") + 
  ylab("Posterior P(Hail = 1| MESH)") + 
  labs(fill = "Confidence", color = "") +
  ggtitle("Hail Posterior") +
  # colors
  scale_fill_paletteer_d("nord::aurora")
# adding true curve if a simulation
if (sim) {
  # currently no transformed simulations so don't need to incorporate that
  p <- p +
    geom_line(aes(y = inv_logit(x_mesh_matrix %*% beta_hail_true), color = "True Probability"))
}
# saving plot
save_plot(p, "fit_hail_curve.png")

# Table for LOO output 
loo_obj <- readRDS(paste(eval_folder, "loo_object.rds", sep = "/"))
loo_obj$estimates %>%
  kbl("markdown", digits = 3) %>% # latex, booktabs = TRUE for latex code
  write_lines(paste(args[5], paste(args[3], "fit_loo_table.txt", sep = "_"), sep = "/"))

# using loo output to create the vehtari plot
base_loo_df_for_plot <- model_data %>% 
  # column of loo pointwise contributions
  add_column(pointwise_loo = loo_obj$pointwise[, 1]) %>%
  # report as factor for color
  mutate(report = as_factor(report))
loo_plot_title <- "Pointwise Leave One Out Log Likelihood Estimate"

# plotting against mesh 
p <- base_loo_df_for_plot %>%
  ggplot(aes(x = mesh, y = pointwise_loo, color = report)) + 
  geom_point() +
  xlab("MESH [mm]") +
  ylab(loo_plot_title) +
  labs(color = "Report")
# saving plot
save_plot(p, "fit_pointwise_loo_mesh.png")

# plotting against population density
p <- base_loo_df_for_plot %>%
  ggplot(aes(x = pop_dens, y = pointwise_loo, color = report)) + 
  geom_point() +
  xlab("Population Density [people / km2]") +
  ylab(loo_plot_title) +
  labs(color = "Report") +
  scale_x_log10()
# saving plot
save_plot(p, "fit_pointwise_loo_dens.png")

# log marginal likelihood as a human readable file
# reading in object
bridge_samples <- readRDS(paste(results_folder, 
              paste(args[3], "bridge_samples_full.rds", sep = "_"), 
              sep = "/"))
# extracting and saving log marginal likelihood
bridge_samples$logml %>%
  write_lines(paste(args[5], 
                    paste(args[3], "fit_log_marginal_likelihood.txt", sep = "_"), 
                    sep = "/"))

### processing all the metric csvs ### 

# creating a table with the leaf metrics 
read_csv(paste(eval_folder, "metrics_leaf.csv", sep = "/")) %>%
  # renaming to minimise table width
  rename_with(~ gsub("metric_", "", .x)) %>%
  # kable to create table
  kbl("markdown", digits = 3) %>% # latex, booktabs = TRUE for latex code
  # writing to file
  write_lines(paste(args[5], 
                    paste(args[3], "fit_leaf_table.txt", sep = "_"), 
                    sep = "/"))

# appending the aggregated metrics to previous file
read_csv(paste(eval_folder, "metrics_aggregated.csv", sep = "/")) %>%
  # again cleaning names
  rename_with(~ gsub("metric_", "", .x)) %>%
  kbl("markdown", digits = 3) %>% # latex, booktabs = TRUE for latex code
  # note the append
  write_lines(paste(args[5], paste(args[3], "fit_leaf_table.txt", sep = "_"), 
                    sep = "/"), 
              append = TRUE)

# plots for the grouped metrics 
make_metrics_plot <- function(file_name, plot_xlab){
  # Function to given a metrics csv file create a useful plot of each metric
  # 
  # Inputs:
  #  file_name  str
  #             name of the file (in the eval directory) to read in to 
  #             provide data for the plot
  #  plot_xlab  str
  #             label for the x axis (the group levels)
  # 
  # Outputs
  #  _          ggplot object
  #             plot of each of the metrics within the different levels
  #             of the group
  
  # col name for plotting 
  col_name <- gsub("\\.csv", "", file_name) %>% 
    gsub("metrics_", "", .)
  
  # create the plot by firstly reading in data
  plotting_df <- read_csv(paste(eval_folder, file_name, sep = "/")) %>%
    # pivotting against each of the metrics
    pivot_longer(c(pvalue, starts_with("metric"))) %>%
    # bulk cleaning of names
    mutate(name = sub("metric_", "", name) %>% gsub("_", " ", .) %>% str_to_sentence()) %>%
    # further tidying of harder to clean names 
    mutate(name = case_when(name == "Auc" ~ "AUC",
                            name == "Log" ~ "Log score",
                            name == "Pb calibration" ~ "Poisson-binomial calibration",
                            name == "Pvalue" ~ "P-value",
                            # if not specified leave unchanged
                            TRUE ~ name)) 
  if(col_name == "group_time_dow") {
    plotting_df <- plotting_df %>%
      mutate(!!col_name := factor(!!sym(col_name), 
				      levels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"), 
				      ordered = TRUE))
  }
  if(col_name == "group_time_month") {
    plotting_df <- plotting_df %>%
      mutate(!!col_name := factor(!!sym(col_name), 
                                      levels = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"), 
                                      ordered = TRUE))
  }
  
  plotting_df %>%
    # plotting against the given column name and coloring by metric
    ggplot(aes(x = !!sym(col_name), y = value, color = name, group = name)) + 
    # the grouping provides the line
    geom_line() +
    # sized by number of observations in level (could easily be number of reports)
    geom_point(aes(size = num_obs_in_level)) +
    # labelling
    labs(color = "Metric", size = "Number of observations in level") +
    ylab("Metric Value") + 
    xlab(plot_xlab) 
}


# list: name is file value, value is xlabel
# contains all the files going into the mega metrics plot
metrics_list <- list("metrics_group_density_binned.csv" = "Population density bins [people / km2]",
                     "metrics_group_density_individual.csv" = "Population density [people / km2]",
                     "metrics_group_density_latitude.csv" = "Latitude [degrees north]",
                     "metrics_group_density_longitude.csv" = "Longitude [degrees east]",
                     "metrics_group_density_urban_rural.csv" = "Urban/rural category",
                     "metrics_group_mesh_fine_bins.csv" = "MESH bins [mm]",
                     "metrics_group_mesh_threshold1.csv" = "MESH grouping [mm]",
                     "metrics_group_mesh_threshold2.csv" = "MESH grouping [mm]",
                     "metrics_group_mesh_threshold3.csv" = "MESH grouping [mm]",
                     "metrics_group_time_dow.csv" = "Day of week",
                     "metrics_group_time_month.csv" = "Month",
                     "metrics_group_time_storm_season.csv" = "Storm Season (September - April",
                     "metrics_group_time_time_of_day.csv" = "Time of day (6-hour bins)",
                     "metrics_group_time_weekend.csv" = "Weekday/Weekend")
# creating plot by iterating over list and its names
metric_plot_list <- imap(metrics_list, \(x, y) make_metrics_plot(file_name = y, plot_xlab = x))
# saving the combination plot
save_plot(wrap_plots(metric_plot_list), "fit_metrics_outofsample.png", width = 30, height = 20)

# reading in the raw probabilites of report for reliability diagram
prob_report_cv <- readRDS(paste(results_folder,
                                paste(args[3], "probability_report_cv.rds", sep = "_"), 
                                sep = "/")
                         )

# reliability diagram
p <- prob_report_cv %>%
  # ensuring report is definitely an integer
  mutate(report = as.integer(as.character(report))) %>%
  # grouping rowwise before taking mean
  rowwise() %>%
  mutate(mean_post_prob_report = mean(c_across(starts_with("prob_report")))) %>%
  # droppign all the now unecessary iterations
  select(-starts_with("prob_report")) %>%
  # grouping by this mean probability
  group_by(post_prob_report_bins = cut(mean_post_prob_report, breaks = 10)) %>%
  # empirical probabilites as well as confidence bounds (under the assumption
  # of independence)
  summarise(mean_emp_prob_report = mean(report), 
            mean_mean_post_prob_report = mean(mean_post_prob_report), 
            lower_bound = qpoibin(0.05, pp = mean_post_prob_report) / n(),
            upper_bound = qpoibin(0.95, pp = mean_post_prob_report) / n(),
            n = n()) %>%
  # plotting 
  ggplot(aes(x = mean_mean_post_prob_report, y = mean_emp_prob_report)) + 
  # error bars
  geom_errorbar(aes(ymin = lower_bound, ymax = upper_bound)) +
  # points
  geom_point(aes(size = n)) + 
  # ideal line
  geom_abline(slope = 1, intercept = 0) +
  # scaling size
  scale_size_continuous(trans = "log10") +
  # labels
  xlab("Posterior mean probability of a report") + 
  ylab("Empirical probability of a report") + 
  labs(size = "Number of\nobservations\nin bin") +
  ggtitle("Reliability Diagram")
# saving the reliability diagram
save_plot(p, "fit_reliability_cv.png")
