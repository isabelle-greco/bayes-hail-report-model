#!/bin/bash

### Driving script to analyse output from a single model
###
### Last modified: 2023-08-23
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=4
#PBS -l mem=15GB
#PBS -l jobfs=30GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/w42+gdata/dk92+scratch/w42
#PBS -l wd

# variables - R modules
export R_MODULE_DIRECTORY="/g/data/dk92/apps/Modules/modulefiles"
export R_MODULE="NCI-data-analysis/2023.02"
export R_LIB_PATH="/scratch/$PROJECT/$USER/rlibs"

# loading R
module use $R_MODULE_DIRECTORY
module load $R_MODULE

# add module directory to path
export R_LIBS=$R_LIBS:$R_LIB_PATH

# getting arguments we need 

# getting path to results
BASE_DIR="/g/data/$PROJECT/$USER/bayesian_paper_data"
# get model name from environment into MODELNAME
PATH_TO_MODEL_RESULTS="$BASE_DIR/model_eval/$MODEL_NAME"

# getting path to data 
DATA_DIR="$PATH_TO_MODEL_RESULTS/data"
DATA_PATH="$DATA_DIR/$(ls $DATA_DIR | grep .csv)"

# saving 
SAVE_DIR="$PATH_TO_MODEL_RESULTS/analysis"
mkdir -p $SAVE_DIR

Rscript analyse_single_model_output.R $R_LIB_PATH $DATA_PATH $MODEL_NAME $PATH_TO_MODEL_RESULTS $SAVE_DIR
