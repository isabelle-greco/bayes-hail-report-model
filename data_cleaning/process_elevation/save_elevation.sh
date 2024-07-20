### Script to save elevation on the hugemem queue
###
### Last modified: 2024-07-20
### Author: Isabelle Greco

#!/bin/bash

# Parameters for the PBS job
#PBS -l ncpus=7
#PBS -l mem=256GB
#PBS -l jobfs=50GB
#PBS -q hugemem
#PBS -P w42
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/w42+gdata/hh5+gdata/rr1+scratch/w42
#PBS -l wd

# changing directory into radar directory
cd process_elevation

# make directory in /g/data
ELEVATION_DIR="/g/data/$PROJECT/$USER/bayesian_paper_data/intermediate_data"
mkdir -p $ELEVATION_DIR

# define file name
ELEVATION_OUTPUT_FILE="$ELEVATION_DIR/brisbane_elevation.csv"

# Loading modules
PYTHON_MODULE_DIRECTORY="/g/data/hh5/public/modules"
PYTHON_MODULE="conda/analysis3"

module use $PYTHON_MODULE_DIRECTORY
module load $PYTHON_MODULE

# pass filename to script to save data 
python save_elevation.py -f $ELEVATION_OUTPUT_FILE
