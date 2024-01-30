#!/bin/bash

### Driver script to perform all the data cleaning.
### Checks for file existence along the way in order to skip
### steps that do not need to be re-done.
###
### Author: Isabelle Greco
### Last modified: 2024-01-25

# variables - R modules
R_MODULE_DIRECTORY="/g/data/dk92/apps/Modules/modulefiles"
R_MODULE="NCI-data-analysis/2023.02"
R_LIB_PATH="/scratch/$PROJECT/$USER/rlibs"

# variables - python modules
export PYTHON_MODULE_DIRECTORY="/g/data/hh5/public/modules"
export PYTHON_MODULE="conda/analysis3"

# variables - data directories
export DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
export DATA_DOWNLOAD_DIRECTORY="$DATA_BASE_DIRECTORY/downloaded_data"
export DATA_INTERMEDIATE_DIRECTORY="$DATA_BASE_DIRECTORY/intermediate_data"
export DATA_CLEAN_DIRECTORY="$DATA_BASE_DIRECTORY/clean_data"

# purging any currently loaded modules
module purge

### gridding ssa
# ssa file name
SSA_FILE_NAME=$(ls $DATA_DOWNLOAD_DIRECTORY | grep ssa)

# current and future data file paths
SSA_DOWNLOAD_PATH="$DATA_DOWNLOAD_DIRECTORY/$SSA_FILE_NAME"
SSA_INTERMEDIATE_PATH="$DATA_INTERMEDIATE_DIRECTORY/$SSA_FILE_NAME"
export SSA_CLEAN_PATH="$DATA_CLEAN_DIRECTORY/$SSA_FILE_NAME"

if [ ! -f $SSA_CLEAN_PATH ]; then
	# load R
	module use $R_MODULE_DIRECTORY
	module load $R_MODULE

	# add module directory to path
	export R_LIBS=$R_LIBS:/$R_LIB_PATH

	# clean the data file
	cat $SSA_DOWNLOAD_PATH | LC_ALL=C sed 's/NULL//g' > $SSA_INTERMEDIATE_PATH

	# runs the r script and saves the data
	# give r library path, data path, and output path as args
	Rscript --vanilla cleaning_ssa/cleaning_ssa.R $R_LIB_PATH $SSA_INTERMEDIATE_PATH \
	$SSA_CLEAN_PATH
fi

# purging currently loaded modules
module purge

### gridding radar
# load python
module use $PYTHON_MODULE_DIRECTORY
module load $PYTHON_MODULE

# variables
TIME_START="2010/01/01"
# end date matches end of SSA data
TIME_END="2023/01/25"
# time step in hours
TIME_STEP=6
# will pick closest radar
RADAR_LOCATION="-27.7 153.2"
# grids in min max step
LATITUDE_GRID="-29.5 -25.5 0.25" 
LONGITUDE_GRID="151 153.5 0.25" 

# make temp dir
export TEMP_DIR=$(mktemp -d)

# check existence 
RADAR_METADATA_FILEPATH="$DATA_DOWNLOAD_DIRECTORY/radar_site_list.csv"
RADAR_GRID_EXIST=$(python grid_radar/check_file_exist.py --startdate $TIME_START \
	--enddate $TIME_END --timestep $TIME_STEP --radarlocation $RADAR_LOCATION \
	--latitude $LATITUDE_GRID --longitude $LONGITUDE_GRID \
	--radarmetadata $RADAR_METADATA_FILEPATH --finaldir $DATA_CLEAN_DIRECTORY)

if [[ $RADAR_GRID_EXIST == 0 ]]; then
	# create gridding command
    DATA_DAILY_FILES_DIRECTORY="$DATA_INTERMEDIATE_DIRECTORY/daily_files"
    mkdir $DATA_DAILY_FILES_DIRECTORY -p
	GRID_RADAR_COMMAND="python grid_max_mesh.py --startdate $TIME_START \
		--enddate $TIME_END \
		--timestep $TIME_STEP \
		--radarlocation $RADAR_LOCATION \
		--latitude $LATITUDE_GRID \
		--longitude $LONGITUDE_GRID \
        --radarmetadata $RADAR_METADATA_FILEPATH \
        --intermediatedir $DATA_DAILY_FILES_DIRECTORY \
        --finaldir $DATA_INTERMEDIATE_DIRECTORY"

	# save gridding command to file
	export COMMAND_FILE="$TEMP_DIR/grid_radar_command.txt"
	echo $GRID_RADAR_COMMAND > $COMMAND_FILE
	# submit script to gadi - script reads command from file
	# pass the file name and the temp dir name
	VARS_TO_PASS="DATA_DOWNLOAD_DIRECTORY,PYTHON_MODULE_DIRECTORY,PYTHON_MODULE,"
	VARS_TO_PASS+="TEMP_DIR,COMMAND_FILE,DATA_CLEAN_DIRECTORY,SSA_CLEAN_PATH"
	qsub -v $VARS_TO_PASS grid_radar/gridding_workflow.sh
fi

# need to sleep in here until file is created
while true ; do
	if [[ $RADAR_GRID_EXIST == 1 ]]; then
		break
	fi
	# wait 2 minutes
	echo "Waiting for radar gridding to process"
	sleep 120
	# recalculate RADAR_GRID_EXIST
	RADAR_GRID_EXIST=$(python grid_radar/check_file_exist.py --startdate $TIME_START \
		--enddate $TIME_END --timestep $TIME_STEP --radarlocation $RADAR_LOCATION \
		--latitude $LATITUDE_GRID --longitude $LONGITUDE_GRID \
		--radarmetadata $RADAR_METADATA_FILEPATH --finaldir $DATA_CLEAN_DIRECTORY)
	# check on process
	qstat
done 


### grid population density data
# while python is loaded do the transition to csv
RADAR_FILE_NAME=$(ls $DATA_CLEAN_DIRECTORY | grep ssa_variable | grep nc | grep -v popdens)
RADAR_FILE_PATH="$DATA_CLEAN_DIRECTORY/$RADAR_FILE_NAME"
OUTPUT_FILE_NAME_END="${RADAR_FILE_NAME%.*}"
RADAR_CSV_FILE_PATH="$DATA_CLEAN_DIRECTORY/$OUTPUT_FILE_NAME_END.csv"
if [ ! -f $RADAR_CSV_FILE_PATH ]; then
	python -c "import xarray as xr; x = xr.open_dataset('$RADAR_FILE_PATH'); x.to_dataframe().to_csv('$RADAR_CSV_FILE_PATH')"
fi

# clear loaded modules
module purge

# load R (again)
module use $R_MODULE_DIRECTORY
module load $R_MODULE

# add module directory to path
export R_LIBS=$R_LIBS:/$R_LIB_PATH

# runs the r script and saves the data
# give r library path, radar data path, pop density data path, and
# output path as args
POP_DENS_FILE_NAME=$(ls $DATA_INTERMEDIATE_DIRECTORY | grep population_density_coarse)
POP_DENS_FILE_PATH="$DATA_INTERMEDIATE_DIRECTORY/$POP_DENS_FILE_NAME"
OUTPUT_FILE_NAME="popdens_$OUTPUT_FILE_NAME_END.csv"
OUTPUT_FILE_PATH="$DATA_CLEAN_DIRECTORY/$OUTPUT_FILE_NAME"
if [ ! -f $OUTPUT_FILE_PATH ]; then
	Rscript --vanilla grid_population/combine_population.R $R_LIB_PATH $RADAR_CSV_FILE_PATH \
	$POP_DENS_FILE_PATH $OUTPUT_FILE_PATH
fi

### filter out the edges of the radar data
NO_EDGES_FILE_PATH="$DATA_CLEAN_DIRECTORY/noedges_$OUTPUT_FILE_NAME"
if [ ! -f $NO_EDGES_FILE_PATH ]; then
	Rscript --vanilla remove_edges/remove_edges.R $R_LIB_PATH $OUTPUT_FILE_PATH \
	$NO_EDGES_FILE_PATH
fi

# delete temp dir
trap "rm -rf $TEMP_DIR" EXIT
