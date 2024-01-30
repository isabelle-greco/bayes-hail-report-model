#!/bin/bash

### Driving script run and evaluate a specified model. Pass name of the model
### via the first argument.
### Assumes that there is a folder in model_eval with the name of this model
### which also contains the relevant .stan file.
###
### To be used as `qsub -v "MODEL_NAME=<your model name>" -N <your_model_name>  run_and_evaluate_model_parallel.sh`
### This command passes the variable MODEL_NAME to the script and sets the job
### name to <your_model_name> to facilitate easy access later.
###
### Last modified: 2023-08-29
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l jobfs=400GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=04:00:00
#PBS -l storage=gdata/w42+gdata/dk92+scratch/w42
#PBS -l wd

# variables - R modules
export R_MODULE_DIRECTORY="/g/data/dk92/apps/Modules/modulefiles"
export R_MODULE="NCI-data-analysis/2023.02"
export R_LIB_PATH="/scratch/$PROJECT/$USER/rlibs"

# purging any currently loaded modules
module purge

# loading R
module use $R_MODULE_DIRECTORY
module load $R_MODULE

# add module directory to path
export R_LIBS=$R_LIBS:$R_LIB_PATH


# directories where data is located
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
if [[ -z "${DATA_UNSPLIT_PATH+set}" ]]
then
  echo "Using observed data:"
  DATA_UNSPLIT_DIRECTORY="$DATA_BASE_DIRECTORY/modelling_data"
  DATA_UNSPLIT_PATH="$DATA_UNSPLIT_DIRECTORY/$(ls $DATA_UNSPLIT_DIRECTORY | grep filtered)"
else
  echo "Using simulated data:"
fi
echo $DATA_UNSPLIT_PATH


# model name and associated directory
# note that the variable MODEL_NAME must be given via -v to qsub
echo "Preparing data, running, and evaluating model $MODEL_NAME"
MODEL_EVAL_DIRECTORY="$DATA_BASE_DIRECTORY/model_eval/$MODEL_NAME"

# check stan model exists
STAN_MODEL="$MODEL_EVAL_DIRECTORY/$MODEL_NAME.stan"
echo $STAN_MODEL
if [ ! -f "$STAN_MODEL" ]; then
  printf '%s\n' "Error: stan model not found." >&2
  exit 1 # error if can't find the model
fi

### create cross-validation folds for the data ###

DATA_SPLIT_DIRECTORY="$MODEL_EVAL_DIRECTORY/data"
mkdir -p $DATA_SPLIT_DIRECTORY

# number of CV folds
NUM_CV_FOLDS=4

# variables of interest
MODEL_VARS=$(cat $MODEL_EVAL_DIRECTORY/${MODEL_NAME}_vars.txt)

# splitting
Rscript cross_validation_split.R $R_LIB_PATH $DATA_UNSPLIT_PATH $DATA_SPLIT_DIRECTORY $NUM_CV_FOLDS $MODEL_NAME $MODEL_VARS

# data path
DATA_SPLIT_PATH="$DATA_SPLIT_DIRECTORY/$(ls $DATA_SPLIT_DIRECTORY | grep $MODEL_NAME | grep .csv)"
echo "Split data:"
echo $DATA_SPLIT_PATH

### fit model ###

MODEL_RESULTS_DIRECTORY="$MODEL_EVAL_DIRECTORY/results"
mkdir -p $MODEL_RESULTS_DIRECTORY

# create model matrix
echo "Preparing data for $MODEL_NAME"
Rscript $MODEL_EVAL_DIRECTORY/${MODEL_NAME}_matrix.R $R_LIB_PATH $DATA_SPLIT_PATH $DATA_SPLIT_DIRECTORY
# hail_matrix.rds and report_matrix.rds now in $DATA_SPLIT_DIRECTORY
# as are the bounds (beta_hail_lower_bounds.rds and beta_report_lower_bounds.rds)

CHAINS=4
ITER=2000
echo "Running model $STAN_PATH"
Rscript run_model_parallel.R $R_LIB_PATH $DATA_SPLIT_PATH $DATA_SPLIT_DIRECTORY $STAN_MODEL $CHAINS $ITER $MODEL_RESULTS_DIRECTORY $MODEL_NAME

### evaluate model ###

MODEL_EVAL_RESULTS_DIRECTORY="$MODEL_EVAL_DIRECTORY/eval"
mkdir -p $MODEL_EVAL_RESULTS_DIRECTORY

echo "Evaluating results in $MODEL_EVAL_RESULTS_DIRECTORY"
Rscript --vanilla calculate_and_aggregate.R $R_LIB_PATH $MODEL_EVAL_RESULTS_DIRECTORY $MODEL_RESULTS_DIRECTORY $MODEL_NAME

### analyse model ###
MODEL_ANALYSIS_DIR="$MODEL_EVAL_DIRECTORY/analysis"
mkdir -p $MODEL_ANALYSIS_DIR

echo "Analysing results in $MODEL_ANALYSIS_DIR"
Rscript analyse_single_model_output.R $R_LIB_PATH $DATA_SPLIT_PATH $MODEL_NAME $MODEL_EVAL_DIRECTORY $MODEL_ANALYSIS_DIR

