#!/bin/bash

### Script to copy sim model files (.stan, .R, and .txt) to Gadi where they
### can be run using run_and_evaluate_model_parallel.sh
### 
### Author: Isabelle Greco
### Last modified: 2024-01-26

# variables - desintation directory names
DATA_BASE_DIRECTORY="/g/data/$PROJECT/$USER/bayesian_paper_data"
MODEL_RUNNING_DIRECTORY="$DATA_BASE_DIRECTORY/model_eval"

# variables - origin directory
MODEL_CURRENT_DIRECTORY="sim_model_files"

# make dirs if don't exist 
mkdir -p $DATA_BASE_DIRECTORY # should have no impact as data should 
			      # already have been created in this directory
mkdir -p $MODEL_RUNNING_DIRECTORY

# copy over the stan files after first making the directory
for F in $MODEL_CURRENT_DIRECTORY/*/*.stan; do
	# extracts the model name by finding hte directory name and removing the first directory
	MODEL_NAME=$(F2=$(dirname $F) && echo ${F2#$MODEL_CURRENT_DIRECTORY/})
	# makes the model directory and then copies in the stan file
	mkdir -p $MODEL_RUNNING_DIRECTORY/$MODEL_NAME && cp "$F" "$_/$MODEL_NAME.stan"
done

# copy over the .R files
# note that we do not first create the directoroy as this has already been done
for F in $MODEL_CURRENT_DIRECTORY/*/*_matrix.R; do
	MODEL_NAME=$(F2=$(dirname $F) && echo ${F2#$MODEL_CURRENT_DIRECTORY/})
	cp "$F" "$MODEL_RUNNING_DIRECTORY/$MODEL_NAME/${MODEL_NAME}_matrix.R"
done

# copy over the .txt files
for F in $MODEL_CURRENT_DIRECTORY/*/*_vars.txt; do
	MODEL_NAME=$(F2=$(dirname $F) && echo ${F2#$MODEL_CURRENT_DIRECTORY/})
	cp "$F" "$MODEL_RUNNING_DIRECTORY/$MODEL_NAME/${MODEL_NAME}_vars.txt"
done

