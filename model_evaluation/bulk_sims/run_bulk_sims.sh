#! /bin/bash
  
### Script to run the bulk simulationsions. 
###
### Last modified: 2023-08-28
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=144
#PBS -l mem=570GB
#PBS -l jobfs=1200GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=14:00:00
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

# directories where data is located
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
DATA_MODELLING_DIRECTORY="$DATA_BASE_DIRECTORY/modelling_data"

# file name
DATA_FILE_NAME="$(ls $DATA_MODELLING_DIRECTORY | grep filtered)"
DATA_FILE_PATH="$DATA_MODELLING_DIRECTORY/$DATA_FILE_NAME"

# directory for bulk sims results 
BULK_SIMS_DIRECTORY="$DATA_BASE_DIRECTORY/bulk_sims"
BULK_SIMS_RESULTS_DIRECTORY="$BULK_SIMS_DIRECTORY/results"
mkdir -p $BULK_SIMS_DIRECTORY
mkdir -p $BULK_SIMS_RESULTS_DIRECTORY

# base dir for stan files
BASE_DIR_STAN="$DATA_BASE_DIRECTORY/model_eval"

# model names
MODEL_NAMES_FILE_PATH="$BULK_SIMS_DIRECTORY/model_names.txt"
ls $BASE_DIR_STAN | grep sim > $MODEL_NAMES_FILE_PATH

# run all the simulations with the given parameters (echoed for debugging)
echo $R_LIB_PATH
echo $BULK_SIMS_RESULTS_DIRECTORY
echo $DATA_FILE_PATH
echo $MODEL_NAMES_FILE_PATH
echo $BASE_DIR_STAN
echo $PBS_NCPUS
Rscript run_bulk_sims.R $R_LIB_PATH $BULK_SIMS_RESULTS_DIRECTORY $DATA_FILE_PATH $MODEL_NAMES_FILE_PATH $BASE_DIR_STAN $PBS_NCPUS 
