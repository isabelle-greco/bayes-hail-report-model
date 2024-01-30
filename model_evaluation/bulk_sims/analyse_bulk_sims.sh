#! /bin/bash

### Script to analyse the output of bulk simuluations run by run_bulk_sims.sh
###
### Last modified: 2023-08-28
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=4
#PBS -l mem=15GB
#PBS -l jobfs=33GB
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

# defining directories for script arguments
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
BULK_SIMS_DIRECTORY="$DATA_BASE_DIRECTORY/bulk_sims"
BULK_SIMS_RES_DIRECTORY="$BULK_SIMS_DIRECTORY/results"

BULK_SIMS_ANALYSIS_DIRECTORY="$BULK_SIMS_DIRECTORY/analysis"
mkdir -p $BULK_SIMS_ANALYSIS_DIRECTORY

# defining number of sims to analyse
NUM_SIMS=100

# calling script
Rscript analyse_bulk_sims.R $R_LIB_PATH $BULK_SIMS_RES_DIRECTORY $NUM_SIMS $BULK_SIMS_ANALYSIS_DIRECTORY

