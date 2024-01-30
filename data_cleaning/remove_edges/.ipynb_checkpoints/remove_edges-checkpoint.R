### Script to filter out the edges of the radar data. 
### Designed to move in a direction towards the coast (if present) else trim
### from all sides.
###
### Author: Isabelle Greco
### Last modified: 2023-03-09

# get args from command line
# 1 = library for r installs if need be
# 2 = location of file to clean
# 3 = location of output data file
args = commandArgs(trailingOnly = TRUE)

# assert correct number of args
if (length(args) != 3) {
  stop(paste("Must give library for R package installs, location of the ",
	     "existing gridded files, and the output location"),
       call. = FALSE)
}

# function to install/load package, and error if version is too low
use <- function(package, install_lib = args[1], version = 0, ...) {
  package <- as.character(substitute(package))
  # gets package version - zero if package doesn't exist
  pver <- tryCatch(packageVersion(package), error = function(cond){0})
  # load package if available of a suitable version else install
  if (pver < as.character(version)){
    install.packages(package, lib = install_lib, repos = "https://cran.csiro.au", dependencies = TRUE)
  }
  library(package, ..., character.only = TRUE)
  # check high enough version again 
  pver <- packageVersion(package)
  if (pver < as.character(version)) {
    stop("Version ", version, " of '", package, 
         "' required, but only ", pver, " is available")
  }
  invisible(pver)
}

# load packages - need dplyr > 1.1.0 for `reframe`
use(dplyr, version = "1.1.0")
use(tidyverse)


# function implementing the algorithm
filter_data <- function(data, x, y, direction = "l2r", width = 1) {
  # Function to algorithmically filter the radar grid in any grid direction
  # (in general towards a coast if present, and in ant direction otherwise)
  # Assumes a contiguous, circular region.
  # 
  # Arguments:
  #     data          tibble with the data
  #     x             centers of bins in the x (horizontal) direction
  #     y             centers of bins in the y (vertical) direction 
  #     direction     the direction in which the filtering will be performed,
  #                   default = "l2r" (left to right)
  #                   other options = "r2l" (right to left), "t2b" (top to 
  #                   bottom), "b2t" (bottom to top), "none" (no directional 
  #                   preference)
  #     width         width of each grid cell, default = 1
  
  # step 1: generate list of cells present in the data
  list_of_cells <- data %>%
    group_by({{x}}, {{y}}) %>%
    summarise(x = mean({{x}}), y = mean({{y}}), .groups = "drop") %>%
    mutate(cells = paste(x, ",", y, sep = ""))
  
  # step 2: delete first column and check valid direction
  if (direction == "l2r" | direction == "none"){
    # treats omnidirectional as the l2r case for simplicity, but could be any
    min_x <- data %>% 
      summarise(m = min({{x}})) %>% 
      pull
    data_clean <- data %>% 
      filter({{x}} > min_x)
  } else if (direction == "r2l") {
    max_x <- data %>% 
      summarise(m = max({{x}})) %>% 
      pull
    data_clean <- data %>% 
      filter({{x}} < max_x)
  } else if (direction == "t2b") {
    max_y <- data %>% 
      summarise(m = max({{y}})) %>% 
      pull
    data_clean <- data %>% 
      filter({{y}} < max_y)
  } else if (direction == "b2t") {
    min_y <- data %>% 
      summarise(m = min({{y}})) %>% 
      pull
    data_clean <- data %>% 
      filter({{y}} > min_y)
  } else {
    stop("Invalid Direction")
  }
 
  # step 3: deal with remaining cells
  if (direction %in% c("l2r", "r2l", "none")) {
      cols <- data_clean %>% 
      reframe(cols = unique({{x}})) %>% 
      pull
    if (direction == "l2r" | direction == "none"){
      final_col <- max(cols)
    } else {
      final_col <- min(cols)
    }
    cols_to_loop <- cols[cols != final_col]
    for (x_val in cols_to_loop) {
      extreme_yval <- list_of_cells %>% 
	filter({{x}} == x_val) %>%
	summarise(min_y = min({{y}}), max_y = max({{y}}))
      cells_to_delete <- list_of_cells %>% 
        filter({{x}} == x_val & 
                 ((! paste(x_val - width, ",", {{y}}, sep = "") %in% cells | 
                    ! paste(x_val + width, ",", {{y}}, sep = "") %in% cells) |
		  ({{y}} == extreme_yval$min_y | {{y}} == extreme_yval$max_y)))
      # delete selected cells 
      data_clean <- data_clean %>%
        filter({{x}} != x_val | 
                 ({{x}} == x_val & ! {{y}} %in% (cells_to_delete %>% 
                                                   pull({{y}}))
                  )
               )
    }
    # step 4 deal with last column
    extreme_yval <- list_of_cells %>%
      filter({{x}} == final_col) %>%
      summarise(min_y = min({{y}}), max_y = max({{y}}))
    if (direction == "l2r") {
      cells_to_delete <- list_of_cells %>% 
        filter({{x}} == final_col & 
                 ((! paste(final_col - width, ",", {{y}}, sep = "") %in% cells) | 
		  ({{y}} == extreme_yval$min_y | {{y}} == extreme_yval$max_y))) 
    } else if (direction == "r2l") {
      cells_to_delete <- list_of_cells %>% 
        filter({{x}} == final_col & 
                 ((! paste(final_col + width, ",", {{y}}, sep = "") %in% cells) | 
		  ({{y}} == extreme_yval$min_y | {{y}} == extreme_yval$max_y))) 
    } else {
      cells_to_delete <- list_of_cells %>% 
        filter({{x}} == final_col)
    }
    # delete selected cells 
    data_clean <- data_clean %>%
      filter({{x}} != final_col | 
               ({{x}} == final_col & ! {{y}} %in% (cells_to_delete %>% 
                                                     pull({{y}}))
                )
             )
  } else {
    # step 3: deal with remaining cells
    rows <- data_clean %>% 
      reframe(rows = unique({{y}})) %>% 
      pull
    if (direction == "t2b"){
      final_row <- min(rows)
    } else {
      final_row <- max(rows)
    }
    rows_to_loop <- rows[rows != final_row]
    for (y_val in rows_to_loop) {
      extreme_xval <- list_of_cells %>% 
	filter({{y}} == y_val) %>%
	summarise(min_x = min({{x}}), max_x = max({{x}}))
      cells_to_delete <- list_of_cells %>% 
        filter({{y}} == y_val & 
                 ((! paste({{x}}, ",", y_val - width, sep = "") %in% cells | 
                    ! paste({{x}}, ",", y_val + width, sep = "") %in% cells) |
		  ({{x}} == extreme_xval$min_x | {{x}} == extreme_xval$max_x)))
      # delete selected cells 
      data_clean <- data_clean %>%
        filter({{y}} != y_val | 
                 ({{y}} == y_val & ! {{x}} %in% (cells_to_delete %>% 
                                                   pull({{x}}))
                  )
               )
    }
    # step 4 deal with last row
    extreme_xval <- list_of_cells %>%
      filter({{y}} == final_row) %>%
      summarise(min_x = min({{x}}), max_x = max({{x}}))
    if (direction == "t2b") {
      cells_to_delete <- list_of_cells %>% 
        filter({{y}} == final_row & 
                 ((! paste({{x}}, ",", final_row + width, sep = "") %in% cells) | 
		  ({{x}} == extreme_xval$min_x | {{x}} == extreme_xval$max_x))) 
    } else if (direction == "b2t") {
      cells_to_delete <- list_of_cells %>% 
        filter({{y}} == final_row & 
                 ((! paste({{x}}, ",", final_row - width, sep = "") %in% cells) | 
		  ({{x}} == extreme_xval$min_x | {{x}} == extreme_xval$max_x))) 
    }
    # delete selected cells 
    data_clean <- data_clean %>%
      filter({{y}} != final_row | 
               ({{y}} == final_row & ! {{x}} %in% (cells_to_delete %>% 
                                                     pull({{x}}))
                )
             )
  }
  data_clean
}

# reading in data 
hail <- read_csv(args[2], col_types = "nnnTncnT")

# applying the filtering
no_edges <- hail %>%
  # removing cells where no MESH ever recorded
  group_by(x_bins, y_bins) %>%
  mutate(num_mesh = sum(!is.na(mesh))) %>%
  ungroup() %>%
  filter(num_mesh > 0) %>%
  select(-num_mesh) %>%
  # now doing the filtering
  filter_data(x_bins, y_bins, width = 0.25)

# saving
no_edges %>%
  write_csv(args[3])
