### Script to, given data and the number of folds require, create the test and 
### training sets as well as the cross validation folds.
### 
### Last modified: 2023-07-10
### Author: Isabelle Greco

# get args from command line
# 1 = library for r installs if need be
# 2 = location of cleaned data to split
# 3 = directory of output data files
# 4 = number of folds for k-fold validation
# 5 = model name
# 6 = list of predictors (comma-separated list)
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 6) {
  stop(paste("Must give library for R package installs, location of cleaned",
             "data to split, directory for output data files, number of folds",
	     "for k-fold cross validation, model name, and a list of predictors."), 
       call. = FALSE)
}


# load libraries + install if need be
for (package in c("tidyverse", "twinning", "stringr")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}


# NOTE: could set.seed() here to ensure the splits are the same each time
# currently want different folds each time

filtered_data <- read_csv(args[2], col_types = "nnTnnnnnnncnffffff")

# parameters for train/test twinning
num_folds <- as.numeric(args[4])

# get columns to include 
columns_to_split <- str_split_1(args[6], ",")

# may be NAs in mesh_* columns - this line sets them to zero
# will also error if any of the columns_to_split aren't actual columns
no_na_data <- filtered_data %>% 
  mutate(across(all_of(columns_to_split), ~ replace_na(.x, 0)))

# selecting relevant columns to split into cv folds
data_to_split <- no_na_data %>%
  select(all_of(columns_to_split))

# calculating fold labels
folds <- multiplet(as.data.frame(data_to_split), k = num_folds, strategy = 2)

# adding to data
fold_data <- no_na_data %>%
  add_column(folds = folds)

# writing out data
# not using model name as filename gets too long to save
write_csv(fold_data, paste(args[3], paste(args[5], "folds.csv", sep = "_"), sep = "/"))

