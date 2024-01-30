### Script to create simulation data for evaluation of the modelling chain
### 
### Last modified: 2023-08-15
### Author: Isabelle Greco

# parse arguments
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 3) {
  stop(paste("Must give library for R package installs, directory in which to",
	     "save the data, and path to the original data."),
       call. = FALSE)
}


# load/install relevant packages
# note putting tidyverse later stops MASS masking the dplyr::select function
for (package in c("tidyverse", "janitor")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# reading in the data
model_data <- read_csv(args[3], col_types = "nnTnnnnnnncnffffff")

# inverse logit function (sigmoid function)
inv_logit <- function(x){
  1 / (1 + exp(-x))
}

## simulation 1: low hail probability, high reporting probability ##

# parameters chosen to give the same expected number of reports as those observed
lhhr_data <- model_data %>%
  mutate(prob_hail = inv_logit(-7 + 0.0932 * mesh),
	 prob_report_given_hail = inv_logit(-2 + 1 * log(pop_dens))) %>%
  mutate(prob_report = prob_hail * prob_report_given_hail) %>%
  add_column(rand = runif(nrow(.))) %>%
  mutate(report = case_when(rand < prob_report ~ 1,
			    TRUE ~ 0)) %>%
  select(-prob_hail, -prob_report_given_hail, -prob_report, -rand)

# saving
lhhr_data %>%
  write_csv(paste(args[2], "sim_lhhr.csv", sep = "/")) 

## test 2: low hail probability, high reporting probability  ##

# again parameters chosen to give the same expected number of reports as those observed
hhlr_data <- model_data %>%
  mutate(prob_hail = inv_logit(-8 + 0.3 * mesh),
	 prob_report_given_hail = inv_logit(-4.258 + 0.5 * log(pop_dens))) %>%
  mutate(prob_report = prob_hail * prob_report_given_hail) %>%
  add_column(rand = runif(nrow(.))) %>%
  mutate(report = case_when(rand < prob_report ~ 1,
			    TRUE ~ 0)) %>%
  select(-prob_hail, -prob_report_given_hail, -prob_report, -rand)

# saving
hhlr_data %>%
  write_csv(paste(args[2], "sim_hhlr.csv", sep = "/")) 

