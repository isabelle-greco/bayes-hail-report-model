### Script to grid the radar data on Gadi's normal queue.
###
### Last modified: 2023-07-03
### Author: Isabelle Greco

#!/bin/bash

# Parameters for the PBS job
#PBS -l ncpus=96
#PBS -l mem=380GB
#PBS -l jobfs=100GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=01:00:00
#PBS -l storage=gdata/w42+gdata/hh5+gdata/rq0+scratch/w42
#PBS -l wd

# changing directory into radar directory
cd grid_radar

# user parameters - read in from command line
RADAR_METADATA_FILENAME=$(ls $DATA_DOWNLOAD_DIRECTORY | grep radar)
RADAR_METADATA_FILEPATH="$DATA_DOWNLOAD_DIRECTORY/$RADAR_METADATA_FILENAME"

# Loading modules
module use $PYTHON_MODULE_DIRECTORY
module load $PYTHON_MODULE

# make a file to pass gridding between scripts
FILE_NAME=$TEMP_DIR/$PBS_JOBID.txt
touch $FILE_NAME

# do the mesh gridding which outputs the file name
bash $COMMAND_FILE > $FILE_NAME

# TODO: remove check
cat $FILE_NAME

# do the ssa gridding - saves data in same folder as mesh data (i.e. $DATA_CLEAN_DIRECTORY)
python include_ssa.py -r $(cat $FILE_NAME | tail -n 1) -s $SSA_CLEAN_PATH -d $DATA_CLEAN_DIRECTORY

