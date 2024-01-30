#!/bin/bash

### Driving script to create the predictors from the full data
###
### Last modified: 2023-07-07
### Author: Isabelle Greco

# variables - R modules
export R_MODULE_DIRECTORY="/g/data/dk92/apps/Modules/modulefiles"
export R_MODULE="NCI-data-analysis/2023.02"
export R_LIB_PATH="/scratch/$PROJECT/$USER/rlibs"

# variables - python modules
export PYTHON_MODULE_DIRECTORY="/g/data/hh5/public/modules"
export PYTHON_MODULE="conda/analysis3"

# purging any currently loaded modules
module purge

# loading R
module use $R_MODULE_DIRECTORY
module load $R_MODULE

# add module directory to path
export R_LIBS=$R_LIBS:$R_LIB_PATH

# loading python
module use $PYTHON_MODULE_DIRECTORY
module load $PYTHON_MODULE

# getting python path
PYTHON_PATH=$(which python)

# directories where data is located
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
DATA_CLEAN_DIRECTORY="$DATA_BASE_DIRECTORY/clean_data"
# make directory for split data - will do nothign if already exists
DATA_MODELLING_DIRECTORY="$DATA_BASE_DIRECTORY/modelling_data"
mkdir -p $DATA_MODELLING_DIRECTORY

# get file name of cleaned data 
# continues assumption that only one set of files in the folder
CLEAN_DATA_FNAME=$(ls $DATA_CLEAN_DIRECTORY | grep noedges)
CLEAN_DATA_PATH="$DATA_CLEAN_DIRECTORY/$CLEAN_DATA_FNAME"

### running the script to separate the data

# running R script
Rscript --vanilla create_predictors.R $R_LIB_PATH $CLEAN_DATA_PATH \
	$DATA_MODELLING_DIRECTORY $PYTHON_PATH
