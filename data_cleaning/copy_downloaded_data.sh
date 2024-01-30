#!/bin/bash

### Script to copy downloaded data into a folder on Gadi for further processing.
### 
### Author: Isabelle Greco
### Last modified: 2024-01-25

# variables - desintation directory names
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
DATA_DOWNLOAD_DIRECTORY="$DATA_BASE_DIRECTORY/downloaded_data"
DATA_INTERMEDIATE_DIRECTORY="$DATA_BASE_DIRECTORY/intermediate_data"

# variables - origin directory
DATA_CURRENT_DIRECTORY="downloaded_data"

# make dirs if don't exist
mkdir -p $DATA_BASE_DIRECTORY
mkdir -p $DATA_DOWNLOAD_DIRECTORY
mkdir -p $DATA_INTERMEDIATE_DIRECTORY

# copy into destination directory
# population density - uncomment if the data is downloaded
# cp $DATA_CURRENT_DIRECTORY/ESRI_GRID_2016 $DATA_DOWNLOAD_DIRECTORY -r # entire directory
# radar site list
cp $DATA_CURRENT_DIRECTORY/radar_site_list.csv $DATA_DOWNLOAD_DIRECTORY
# ssa reports
cp $DATA_CURRENT_DIRECTORY/ssa_hail_060101-230125.csv $DATA_DOWNLOAD_DIRECTORY
# processed population density - uncomment if intermediate files created
# cp $DATA_CURRENT_DIRECTORY/population_density_new_crs.tif $DATA_INTERMEDIATE_DIRECTORY
# cp $DATA_CURRENT_DIRECTORY/population_density_zero.tif $DATA_INTERMEDIATE_DIRECTORY
cp $DATA_CURRENT_DIRECTORY/population_density_coarse_grid.tif $DATA_INTERMEDIATE_DIRECTORY

