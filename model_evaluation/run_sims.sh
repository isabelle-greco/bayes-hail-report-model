#! /bin/bash

### Script to create simulated data and submit the simulation runs to queues
###
### Last modified: 2023-08-25
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=4
#PBS -l mem=15GB
#PBS -l jobfs=30GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=00:30:00
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

# directory for fake data 
SIM_DATA_DIRECTORY="$DATA_MODELLING_DIRECTORY/sim_data"
mkdir -p $SIM_DATA_DIRECTORY

# create perfect data
Rscript create_sim_data.R $R_LIB_PATH $SIM_DATA_DIRECTORY $DATA_FILE_PATH 

## submit jobs to run models and evaluate ##

# lhhr #

# lhhr with the naive N(0, 1) priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_lhhr_naive,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_mesh_report_dens_lhhr_naive run_and_evaluate_model_parallel.sh

# lhhr with informative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_lhhr_informative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_mesh_report_dens_lhhr_informative run_and_evaluate_model_parallel.sh

# lhhr with uninformative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_lhhr_uninformative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_mesh_report_dens_lhhr_uninformative run_and_evaluate_model_parallel.sh

# lhhr with uninformative gamma priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_lhhr_uninformative_gamma,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_mesh_report_dens_lhhr_uninformative_gamma run_and_evaluate_model_parallel.sh

# lhhr with misinformative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_lhhr_misinformative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_mesh_report_dens_lhhr_misinformative run_and_evaluate_model_parallel.sh

# std lhhr with N(0, 1) priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_lhhr_std,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_std_mesh_report_std_dens_lhhr_std run_and_evaluate_model_parallel.sh

# std lhhr with N(0, 0.1) priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_lhhr_std_small,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_std_mesh_report_std_dens_lhhr_std_small run_and_evaluate_model_parallel.sh

# std lhhr with gamma priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_lhhr_gamma,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_lhhr.csv -N sim_hail_std_mesh_report_std_dens_lhhr_gamma run_and_evaluate_model_parallel.sh

# hhlr # 

# hhlr with naive N(0, 1) priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_hhlr_naive,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_mesh_report_dens_hhlr_naive run_and_evaluate_model_parallel.sh

# hhlr with informative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_hhlr_informative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_mesh_report_dens_hhlr_informative run_and_evaluate_model_parallel.sh

# hhlr with uninformative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_hhlr_uninformative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_mesh_report_dens_hhlr_uninformative run_and_evaluate_model_parallel.sh

# hhlr with uninformative priors gamma
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_hhlr_uninformative_gamma,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_mesh_report_dens_hhlr_uninformative_gamma run_and_evaluate_model_parallel.sh

# hhlr with misinformative priors
qsub -v MODEL_NAME=sim_hail_mesh_report_dens_hhlr_misinformative,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_mesh_report_dens_hhlr_misinformative run_and_evaluate_model_parallel.sh

# std hhlr with N(0, 1) priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_hhlr_std,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_std_mesh_report_std_dens_hhlr_std run_and_evaluate_model_parallel.sh

# std hhlr with N(0, 0.1) priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_hhlr_std_small,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_std_mesh_report_std_dens_hhlr_std_small run_and_evaluate_model_parallel.sh

# std hhlr with gamma priors
qsub -v MODEL_NAME=sim_hail_std_mesh_report_std_dens_hhlr_gamma,DATA_UNSPLIT_PATH=$SIM_DATA_DIRECTORY/sim_hhlr.csv -N sim_hail_std_mesh_report_std_dens_hhlr_gamma run_and_evaluate_model_parallel.sh

