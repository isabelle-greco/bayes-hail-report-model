#! /bin/bash

### Script to create simulated data and submit the simulation runs to queues
###
### Last modified: 2023-08-29
### Author: Isabelle Greco

# Parameters for the PBS job
#PBS -l ncpus=1
#PBS -l mem=1GB
#PBS -l jobfs=3GB
#PBS -q normal
#PBS -P w42
#PBS -l walltime=00:10:00
#PBS -l storage=gdata/w42+gdata/dk92+scratch/w42
#PBS -l wd

## log pop_dens models ##

# untransformed mesh
echo "hail_std_mesh_report_std_dens"
qsub -v MODEL_NAME=hail_std_mesh_report_std_dens -N hail_std_mesh_report_std_dens run_and_evaluate_model_parallel.sh

# sqrt mesh
echo "hail_std_sqrt_mesh_report_std_dens"
qsub -v MODEL_NAME=hail_std_sqrt_mesh_report_std_dens -N hail_std_sqrt_mesh_report_std_dens run_and_evaluate_model_parallel.sh

# trans mesh
# ...2 separates it from earlier unsatisfactory verions
echo "hail_trans_std_mesh_report_std_dens2"
qsub -v MODEL_NAME=hail_trans_std_mesh_report_std_dens2 -N hail_trans_std_mesh_report_std_dens2 run_and_evaluate_model_parallel.sh

## no log pop_dens models ##

# untransformed mesh
echo "hail_std_mesh_report_std_nolog_dens"
qsub -v MODEL_NAME=hail_std_mesh_report_std_nolog_dens -N hail_std_mesh_report_std_nolog_dens run_and_evaluate_model_parallel.sh

# sqrt mesh
echo "hail_std_sqrt_mesh_report_std_nolog_dens"
qsub -v MODEL_NAME=hail_std_sqrt_mesh_report_std_nolog_dens -N hail_std_sqrt_mesh_report_std_nolog_dens run_and_evaluate_model_parallel.sh

# trans mesh
# same meaning for number 2 as above
echo "hail_trans_std_mesh_report_std_nolog_dens2"
qsub -v MODEL_NAME=hail_trans_std_mesh_report_std_nolog_dens2 -N hail_trans_std_mesh_report_std_nolog_dens2 run_and_evaluate_model_parallel.sh

## trans pop_dens models ##

# untransformed mesh
echo "hail_std_mesh_report_trans_std_dens"
qsub -v MODEL_NAME=hail_std_mesh_report_trans_std_dens -N hail_std_mesh_report_trans_std_dens run_and_evaluate_model_parallel.sh

# sqrt mesh
echo "hail_std_sqrt_mesh_report_trans_std_dens"
qsub -v MODEL_NAME=hail_std_sqrt_mesh_report_trans_std_dens -N hail_std_sqrt_mesh_report_trans_std_dens run_and_evaluate_model_parallel.sh

# trans mesh
echo "hail_trans_std_mesh_report_trans_std_dens"
qsub -v MODEL_NAME=hail_trans_std_mesh_report_trans_std_dens -N hail_trans_std_mesh_report_trans_std_dens run_and_evaluate_model_parallel.sh

