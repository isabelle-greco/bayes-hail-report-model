### Script to combine the population density grid with the existing
### radar and SSA grid files.
###
### Author: Isabelle Greco
### Last modified: 2023-03-14

# get args from command line
# 1 = library for r installs if need be
# 2 = location of radar/ssa file
# 3 = location of population density file
# 4 = location of output data file
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 4) {
  stop(paste("Must give library for R package installs, radar/ssa location,", 
	     "population density location, and output location."), 
       call. = FALSE)
}

# load libraries + install if need be
for (package in c("tidyverse", "janitor", "lubridate", "raster", "tabularaster")) {
    if (!require(package, character.only = TRUE)) {
	# if not installed, install in the directory given by first argument
        install.packages(package, lib = args[1], 
                         repos = "https://cran.csiro.au", dependencies = TRUE)
        library(package, character.only = TRUE)
    }
}

# reading in population density data
pop_dens <- raster(args[3]) %>%
  as_tibble(xy = TRUE) %>%
  rename(pop_dens = cellvalue, x_bins = x, y_bins = y) %>%
  dplyr::select(-cellindex)

# reading in ssa/radar data
hail <- read_csv(args[2], col_types = "TnnncnT") %>%
  clean_names()

# merging
joined_data <- pop_dens %>%
  # joining to the hail data
  right_join(y = hail, by = c("x_bins", "y_bins"), multiple = "all")

# saving
joined_data %>%
  write_csv(args[4])
