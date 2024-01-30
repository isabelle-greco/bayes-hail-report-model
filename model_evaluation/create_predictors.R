### Script to create all predictors from the full data set and filter
### the data down to the form used for modelling.
### 
### Last modified: 2023-07-40
### Author: Isabelle Greco

# get args from command line
# 1 = library for r installs if need be
# 2 = location of cleaned data to filter
# 3 = directory of output data files
# 4 = directory of python to use
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 4) {
  stop(paste("Must give library for R package installs, location of cleaned",
             "data to filter, directory for output data files, and location",
	     "of python installation to use."), 
       call. = FALSE)
}


# load libraries + install if need be
for (package in c("tidyverse", "stringr", "lubridate", "reticulate")) {
  if (!require(package, character.only = TRUE)) {
    # if not installed, install in the directory given by first argument
    install.packages(package, lib = args[1], 
                     repos = "https://cran.csiro.au", dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# read in data
full_data <- read_csv(args[2], col_types = "nnnTncnT") %>%
  # making a factor report variable
  mutate(report = case_when(is.na(report_date) ~ 0,
                            TRUE ~ 1)) %>%
  mutate(report = as_factor(report)) %>%
  # selecting only entries with valid MESH and population density
  filter(pop_dens > 0) %>%
  drop_na(mesh) %>%
  # report_date column superseded by time_bins
  select(-report_date)

# then add extra predictors
# get surrounding mesh (four cardinal directions)
# table with the coordinates and mesh to reference
mesh_table <- full_data %>%
  select(x_bins, y_bins, time_bins, mesh)

# calculate bin width - assumes uniform bin width and at least one pair of adjacent bins
# spatial
x_binwidth <- full_data %>%
  pull(x_bins) %>% # could also use latitude bins here
  unique() %>% # get unique entries
  sort() %>% # sort them in ascending order
  diff() %>% # difference between consecutive elements
  min() # smallest difference between entries is binwidth
# very, very unlikely that there does not exist at least one pair of adjacent 
# zonal bins where mesh is measured at least once

y_binwidth <- full_data %>%
  pull(y_bins) %>% # could also use latitude bins here
  unique() %>% # get unique entries
  sort() %>% # sort them in ascending order
  diff() %>% # difference between consecutive elements
  min() # smallest difference between entries is binwidth
# very, very unlikely that there does not exist at least one pair of adjacent 
# meridional bins where mesh is measured at least once

# temporal
temporal_bindwidth <- full_data %>%
  pull(time_bins) %>% # could also use latitude bins here
  unique() %>% # get unique entries
  sort() %>% # sort them in ascending order
  diff() %>% # difference between consecutive elements
  min() # smallest difference between entries is binwidth
# very unlikely that in 6.5 storm seasons there was no day where mesh was 
# measured both in the morning and the afternoon


# displaces table and merges new predictor, keeping all rows in original data 
add_spatial_pred <- mesh_table %>%
  mutate(x_bins = x_bins + x_binwidth) %>%
  rename(mesh_west = mesh) %>%
  right_join(full_data, by = join_by(x_bins, y_bins, time_bins))
# builds on the new data frame, adding each predictor in the cardinal directions
# note: these predictors will have an NA in there if the data isn't available 
# (e.g. edge of domain, no neighbouring measurement)
add_spatial_pred <- mesh_table %>%
  mutate(x_bins = x_bins - x_binwidth) %>%
  rename(mesh_east = mesh) %>%
  right_join(add_spatial_pred, by = join_by(x_bins, y_bins, time_bins))
add_spatial_pred <- mesh_table %>%
  mutate(y_bins = y_bins + y_binwidth) %>%
  rename(mesh_south = mesh) %>%
  right_join(add_spatial_pred, by = join_by(x_bins, y_bins, time_bins))
add_spatial_pred <- mesh_table %>%
  mutate(y_bins = y_bins - y_binwidth) %>%
  rename(mesh_north = mesh) %>%
  right_join(add_spatial_pred, by = join_by(x_bins, y_bins, time_bins))
# not strictly a spatial predictor but a temporally displaced MESH predictor
# calculating this now can use overnight MESH measurements if they exist
add_spatial_pred <- mesh_table %>%
  mutate(time_bins = time_bins + temporal_bindwidth) %>%
  rename(mesh_morning = mesh) %>%
  right_join(add_spatial_pred, by = join_by(x_bins, y_bins, time_bins))

add_time_pred <- add_spatial_pred %>%
  # temporal predictors
  mutate(dow = wday(time_bins, label = TRUE),
         month = month(time_bins, label = TRUE),
         storm_season = as_factor(
           case_when(month(time_bins) >= 9 ~ paste(year(time_bins), year(time_bins) - 1999, sep = "/"),
                     TRUE ~ paste(year(time_bins) - 1, year(time_bins) - 2000, sep = "/"))
         ),
         # putting time of day into aest
         time_of_day = as_factor(hour(time_bins) + 10)) %>%
  mutate(weekend = as_factor(case_when(
    dow == "Sat" | dow == "Sun" ~ "weekend",
    TRUE ~ "weekday")))

# get file name - split string, take last to get rid of directory
filename <- str_split_1(args[2], "/") %>% 
  tail(n = 1)
write_csv(add_time_pred, paste(args[3], paste("allpredictors", filename, sep = "_"), sep = "/"))

# now that we're creating data for model, filter out irrelevant data 
filtered_data <- add_time_pred %>%
  # only data 2010 - 2016
  filter(time_bins < "2016-05-01") %>%
  # select only storm season
  filter(month(time_bins) %in% c(1, 2, 3, 4, 9, 10, 11, 12)) %>%
  # select only afternoon/evening bins
  filter(hour(time_bins) %in% c(3, 9))


# use output from which python as an argument
use_python(args[4])
source_python("../data_cleaning/grid_radar/gridding_util.py")

# remove datetime_
root_name <- paste("filtered", "allpredictors", filename %>% 
                     str_split_1("datetime_") %>% # removing datetime
                     head(1), # taking everything in front of it
                   sep = "_")
# extract radar id
radar_id <- filename %>% 
  str_split_1("radar_") %>% # splitting at radar/dropping radar
  tail(1) %>% # taking the bit after radar 
  str_split_1("_") %>% # splitting with the underscore
  head(1) %>% # taking the first entry - the radar ID
  as.integer # string to integer 
# get binning info
relevant_components <- filename %>% 
  str_split_1(paste("radar", radar_id, "", sep = "_")) %>%
  tail(1) %>%
  str_split_1("_")
# know order is time, x/longitude, y/latitude
axis_names <- relevant_components[3:(which(relevant_components == "binned")-1)]
binned <- rep(TRUE, length(axis_names))
if (length(axis_names >= 1)) {
  t_bins <- seq(min(filtered_data$time_bins) - 0.5 * temporal_bindwidth, 
                max(filtered_data$time_bins) + 0.5 * temporal_bindwidth, 
                by = temporal_bindwidth)
}
if (length(axis_names >= 2)) {
  x_bins <- seq(min(filtered_data$x_bins) - 0.5 * x_binwidth, 
                max(filtered_data$x_bins) + 0.5 * x_binwidth, 
                by = x_binwidth)
}
if (length(axis_names >= 3)) {
  y_bins <- seq(min(filtered_data$y_bins) - 0.5 * y_binwidth, 
                max(filtered_data$y_bins) + 0.5 * y_binwidth, 
                by = y_binwidth)
}
# create list of arrays
np <- import("numpy", convert = FALSE)
t_bins_array <- np$array(t_bins)
x_bins_array <- np$array(x_bins)
y_bins_array <- np$array(y_bins)
# pass to function
grid_params_filename <- create_file_name(binned, axis_names, list(t_bins_array, x_bins_array, y_bins_array), radar_id)


new_filename <- paste(paste(root_name, grid_params_filename, sep = ""), "csv", sep = ".")
write_csv(filtered_data, paste(args[3], new_filename, sep = "/"))

