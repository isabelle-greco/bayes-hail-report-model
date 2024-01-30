### Script to produce plots and evaluate results from the bulk simulations.
###
### Last Modified: 2023-08-28
### Author: Isabelle Greco

# get args from command line
# 1 = library for r installs if need be
# 2 = directory in which to find the bulk sims results
# 3 = number of simulations run
# 4 = directory in which to save analyses
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 4) {
  stop(paste("Must give library for R package installs, directory in which to find",
             "bulk simulation results, number of simulations run directory in which", 
             "to save analyses."), 
       call. = FALSE)
}


# load libraries + install if need be
for (package in c("tidyverse", "paletteer", "patchwork")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# setting plots to large
options(repr.plot.width = 45, repr.plot.height = 30)

### reading data and defining parameters ###

# reading in results from location given by args[2]
sims_res_dir <- args[2]
# number of simulations given by args[3]
bulk_res <- map(1:args[3], 
                \(x) readRDS(paste0(sims_res_dir, "/bulk_sims_res_", x, ".rds"))
               )

# creating matrices for evaluating and plotting etc. 
x_mesh_matrix <- readRDS(
  paste(sims_res_dir, "x_mesh_matrix.rds", sep = "/")
)
x_log_pop_dens_matrix <- readRDS(
  paste(sims_res_dir, "x_log_pop_dens_matrix.rds", sep = "/")
)

# defining the true parameters 
params_list <- list("lhhr_beta_hail" = c(-7, 0.0932),
                    "lhhr_beta_report" = c(-2, 1),
                    "hhlr_beta_hail" = c(-8, 0.3),
                    "hhlr_beta_report" = c(-4.258, 0.5))

### utility functions ###

# inverse logit
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

# function to analyse simulations
analyse_simulations <- function(quantity, variable, model, inner_ci = 0.5,
                                outer_ci = 0.9){
  # Function to analyse a given quantity (e.g. mean, quantile) of a model (e.g.
  # informative, uninformative, misspecified) on a given data set (e.g. lhhr or
  # hhlr) for a given probability curve (e.g. reporting given hail or hail)
  #
  # Inputs:
  #  quantity    str
  #              "mean" or a percentile "0.025", "0.05", "0.25", "0.5", "0.75", 
  #               "0.95", "0.975"
  #  variable    str
  #              "report" (= reporting model) or "hail" (= hail model)
  #  model       str
  #              of form sim_hail_(std)_mesh_report_(std)_dens_(lhhr/hhlr)_
  #              ((mis)(un)informative/std(_small))_(cv/all)
  #  inner_ci    float in [0, 1]
  #              confidence of the inner ci in the plot, must be less than outer_ci
  #  outer_ci    float in [0, 1]
  #              confidence of the outer ci in the plot, must be greater than inner_ci
  #
  # Outputs:
  #  _            A plot of the statistics for that quantity (50% and 90% confidence 
  #               intervals) along side the true value of the function for comparison
  
  # splitting model name
  split_model_name <- model %>%
    str_split_1("_")
  # extracting useful quantities 
  data_set <- if ("lhhr" %in% split_model_name) "lhhr" else "hhlr"
  cv_status <- if ("cv" %in% split_model_name) "cv" else "all"
  # data name
  data_name <- paste(data_set, "res", cv_status, sep = "_")
  # quantity as a string for column name
  quantity_str <- ifelse(quantity == "mean", quantity,
                      paste(gsub(pattern = "\\.", 
                                 replacement = "_", 
                                 x = 100 * as.numeric("0.025")
                                ), 
                            "percent", 
                            sep = "_")
                        )
  # dropping cv/alll
  column_name <- paste(model, variable, quantity_str, sep = "_")
    
  # selecting right x vector
  # second column of the matrices
  x_matrix <- if (variable == "hail") x_mesh_matrix else x_log_pop_dens_matrix
    
  # transforming back to original scale
  # again using oneline if/else to intentionally avoid the vectorisation
  x_trans <- if (variable == "hail") x_matrix[, 2] else exp(x_matrix[, 2])
  # note x_mesh is untransformed

  # selecting true a and b values
  beta_true <- params_list[[paste(data_set, "beta", variable, sep = "_")]]

  # getting the relevant column from each simulation
  extracting_column <- map(bulk_res, 
                           # for each bulk sim grab the lhhr results corresponding to the 
                           # right sim data set and cv/all regime and then select results 
                           # for the appropriate model
                           \(x) pluck(x, data_name, model) %>% 
                           # and then select the columns with the posterior mean 
                           select(all_of(column_name))
                           )
  # bind each one column tibble (with identical names - ignore the warning)
  results_to_plot <- suppressMessages(
    list_cbind(extracting_column, name_repair = "unique")
  ) %>%
    # make the data frame rowwise
    rowwise() %>%
    # get the bulk mean of the posterior mean, and the quantiles for plotting
    summarise(bulk_mean = mean(c_across(everything())),
              bulk_inner_lower = quantile(c_across(everything()), (1 - inner_ci) / 2),
              bulk_inner_upper = quantile(c_across(everything()), (1 + inner_ci) / 2),
              bulk_outer_lower = quantile(c_across(everything()), (1 - outer_ci) / 2),
              bulk_outer_upper = quantile(c_across(everything()), (1 + outer_ci) / 2)
              )
    

  # axis labels
  x_lab <- ifelse(variable == "hail", "MESH [mm]", 
                  "Population Density [people/km2]")
  y_lab <- ifelse(variable == "hail", "P(Hail = 1| MESH)", 
                  "P(Report = 1| Hail = 1, Pop. Dens.)")
  # creating title in pieces
  # quantity to mention
  title_quantity <- ifelse(quantity == "mean", quantity, 
                           paste(100 * as.numeric(quantity), "% percentile", sep = ""))
  # if standardised or not
  standardised <- ifelse("std" %in% split_model_name, "standardised", "unstandardised")
  # what the priors are 
  priors <- case_when("informative" %in% split_model_name ~ "informative",
                      "uninformative" %in% split_model_name ~ "uninformative",
                      "misinformative" %in% split_model_name ~ "misinformative",
                      "naive" %in% split_model_name ~ "naive",
                      "small" %in% split_model_name ~ "N(0, 0.1^2)",
                      "gamma" %in% split_model_name ~ "gamma and normal",
                      TRUE ~ "N(0, 1^2)")
  # putting all together into title 
  title <- paste("Posterior", paste(title_quantity, ":", sep = ""),
                 toupper(data_set), "data,",
                 standardised, "model with",
                 priors, "priors")

  # confidence bounds 
  label_inner <- paste(inner_ci * 100, "%", sep = "")
  label_outer <- paste(outer_ci * 100, "%", sep = "")
    
  # plotting!
  p <- results_to_plot %>%
    # true probability
    add_column(true = inv_logit(x_matrix %*% beta_true)) %>%
    # plotting
    ggplot(aes(x = x_trans)) + 
    # mean for the quantity
    geom_line(aes(y = bulk_mean, color = paste("Posterior", title_quantity))) + 
    # true value
    geom_line(aes(y = true, color = "True probability")) +
    # outer confidence bounds
    geom_ribbon(aes(ymin = bulk_outer_lower, 
                    ymax = bulk_outer_upper, fill = label_outer), 
                alpha = 0.3) + 
    # inner confidence bounds 
    geom_ribbon(aes(ymin = bulk_inner_lower, 
                    ymax = bulk_inner_upper, fill = label_inner), 
                alpha = 0.5) +
    # axis labels
    xlab(x_lab) +
    ylab(y_lab) +
    labs(fill = "Confidence") + 
    # title
    ggtitle(title) + 
    # color
    scale_fill_paletteer_d("nord::aurora") + 
    # remove legend for color 
    labs(color = "")

  if (variable == "report") {
    p + 
      scale_x_log10()
  } else {
    p
  }
}


# function to plot all models for given name and variable 
plot_all_models <- function(name, var, quantity = "mean") {
  # Function to plot all the models for a given data set (e.g. lhhr_res_cv
  # or hhlr_res_all)
  # 
  # Inputs:
  #  name      str
  #            name of thee group of simulations to plot in form 
  #            (lhhr/hhlr)_res_(cv/all)
  #  var       str
  #            variable to plot. one of "hail" or "report"
  #  quantity  str
  #            quantity to plot. can take values "mean", or a percentile
  #            "0.025", "0.05", "0.25", "0.5", "0.75", "0.95", "0.975"
  # 
  # Outputs:
  #  _         ggplot object
  #            each subplot plot is for a particular mode configuration

  # getes the names of the available model results 
  model_names <- names(bulk_res[[1]][[name]]) # could use any list element
  # using purrr iterates over ache name and applys the analysis function
  p <- map(model_names, 
           \(x) analyse_simulations(quantity, var, x)) %>%
    # wraps plots into one figure and collects legends 
    wrap_plots(guides = "collect", nrow = 3)
  filename <- paste("bulk_sims_res", name, var, quantity, "fig.png", sep = "_")
  ggsave(plot = p,
         filename = paste(args[4], filename, sep = "/"),
         device = png,
	 width = 45, 
	 height = 30)
}

# iterate over both variable sand all four model combinations 
map(c("hail", "report"), 
    \(y) map(names(bulk_res[[1]]), \(x) plot_all_models(x, y))
   )
