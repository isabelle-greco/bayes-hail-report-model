### Script to create the relevant matrices for the hail_std_sqrt_mesh_report_std_dens model.
###
### Last modified: 2023-08-28
### Author: Isabelle Greco

# parse arguments
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 3) {
  print(args)
  stop(paste("Must give library for R package installs, path to data to use, and ",
	     "directory in which to save the resulting matrices."),
       call. = FALSE)
}


# load/install relevant package
if (!require("tidyverse", character.only = TRUE)) {
  # if not installed, install in the directory given by first argument
  install.packages("tidyverse", lib = args[1], 
                   repos = "https://cran.csiro.au", dependencies = TRUE)
  library("tidyverse", character.only = TRUE)
}

# read in data
# note we're reading in `folds` as an integer
model_data <- read_csv(args[2], col_types = "nnTnnnnnnncnffffffn") %>%
  # month as an ordered factor
  mutate(month = factor(month, levels = c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"), ordered = TRUE))


# create matrices
# matrix for hail
model_frame_hail <- model_data %>% 
  model.frame(report ~ mesh + folds, data = .)
model_matrix_hail <- model.matrix(report ~ mesh + folds, data = model_frame_hail)

# matrix for report
model_frame_report <- model_data %>% 
  # ensures report is converted from factor correctly
  mutate(report = case_when(report == 0 ~ 0,
                            report == 1 ~ 1)) %>%  
  model.frame(~ pop_dens + month + report + folds, data = .)
model_matrix_report <- model.matrix(~ pop_dens + month + report + folds, data = model_frame_report, 
				    contrasts = list(month = "contr.treatment"))

# saving matrices
saveRDS(model_matrix_hail, file = paste(args[3], "hail_matrix.rds", sep = "/"))
saveRDS(model_matrix_report, file = paste(args[3], "report_matrix.rds", sep = "/"))

# also creating lower bounds
beta_hail_lower_bounds <- c(-Inf, 0.0, -Inf)
beta_report_lower_bounds <- c(-Inf, 0.0, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf)

# saving
saveRDS(beta_hail_lower_bounds, file = paste(args[3], "beta_hail_lower_bounds.rds", sep = "/"))
saveRDS(beta_report_lower_bounds, file = paste(args[3], "beta_report_lower_bounds.rds", sep = "/"))

