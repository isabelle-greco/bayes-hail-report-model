# Model evaluation

Running code in this directory assumes that code in the [data\_cleaning](../data_cleaning) directory has been run already.

The first step before running any models is to execute the final step of data cleaning (equivalently the first step of evaluating the models) with the [create\_predictors.sh](create_predictors.sh) script.
This script creates all necessary extra predictors from the data and filters the data to the storm seasons on which this study focusses.

Most of the remaining files are called by [run\_and\_evaluate\_model\_parallel.sh](run_and_evaluate_model_parallel.sh) which, unsurprisingly, compiles the model, applies it to the cleaned data (both in its entirety and using cross validation), evaluates it, and generates a range of diagnostics and plots.
Recommended usage is by submitting a job on Gadi via:
```
qsub -v "MODEL_NAME=<your_model_name>" -N <descriptive_name_for_log_files> run_and_evaluate_model_parallel.sh
```

Each model has its own directory which must first be copied to `/g/data` on Gadi before the model is used (the script [copy\_model\_files.sh](copy_model_files) copies all the models in [all\_model\_files](all_model_files) into a `model_evaluation` directory in `/g/data`).
The naming system aims to uniquely identify models: 
```
hail_<transform>_<hail predictor1>_..._report_<transform>_<report_predictor1>_..._<prior_on_transform_parameters>_<priors_on_coefficients_of_first_predictors>_<other_priors>
```
The absence of an indication denotes the default for that category.
For example, the model examinined in detail in Section 4 used a mix of normal (consdidered as the default) and gamma priors along with a very weakly informative prior on the Yeo-Johnson transform parameters.
Hence this model was termed `hail_trans_std_mesh_report_trans_std_dens_yjwide_gamma`.
When the inclusion of (standardised) MESH in the reporting model was tested it became `hail_trans_std_mesh_report_trans_std_dens_std_mesh_yjwide_gamma`.
Note that the two models followed by '2' indicate simply that the first version of the model (not included) had an error and the second is correct.
Several scripts ([run\_base\_models.sh](run_base_models.sh), [run\_perturbed\_prior\_models.sh](run_perturbed_prior_models.sh)) are provided which submit multiple jobs at a time to run multiple models that one may be interested in comparing.
The [analyse\_single\_model\_output.sh](analyse_single_model_output.sh), employed in the same way as [run\_and\_evaluate\_model\_parallel.sh](run_and_evaluate_model_parallel.sh), is provided to run the analysis script ([analyse\_single\_model\_output.R](analyse_single_model_output.R)) over existing model output.

In each appropriately named directory are three files essential to running and evaluating the model:

1. The `.stan` file contains the description of the Bayesian model for use with `RStan` or other pathways to `Stan` (e.g. [PyStan](https://pystan.readthedocs.io/en/latest/)). The structure of these files depends primarily upon whether or not the Yeo-Johnson transform is employed. Generally the differences between different model files is only small variations in the exact priors used.
2. The `.R` file descibes the construction of the predictor matrices. These again change very little between models, varying mostly with the use (or not) of the Yeo-Johnson transform and the inclusion of additional predictors.
3. The `.txt` file lists the variables required for the model (including the indicator of reports) in plain text. This file onlychanges when additional predictors are included.

Note that the variable naming in the stan file does not perfectly align with the notation introduced in the paper due to changes in written notation after the analysis was complete.
Notably, the $\alpha$ parameters are named `beta_hail_raw` and the $\beta$ parameters named `beta_report_raw`.
Further, `beta_hail` and `beta_report` are defined as coefficients of the unstandardised Yeo-Johnson predictors to facilitate the analysis in [figures](../figures) directory.
Despite these semantic differences and additional calculations, the models are identical mathematically to those described in the paper.

Once the models are run, these directories are populated with the complied Stan model (a non-human readable `.rds` file used to save R variables) along with four subdirectories:

1. `data` which contains a `.csv` of the data with assigned cross validation folds,`.rds` files with the predictors matrices, and `.rds` files giving bounds on the parameters which are also passed to the stan model,
2. `results` which contains samples output directly from the MCMC sampling, again as `.rds` files,
3. `eval` which contains `.csv` files giving the results of the different evaluation metrics applied over different groupings of the data, and
4. `analysis` which contains a mix of plots and formatted `.txt` tables summarising the less accessible output in `results` and `eval`. Approximately half the files focus upon evaluataing the MCMC algorithm and the other upon the quality of model fit.

The `analysis` folder is primarily for visual human inspection, whilst the figures and plots in TODO:FOLDER mainly leverage the `results` and `eval` folders.
One will notice that many of the metrics calculated in `eval` and several other model evaluation metrics (e.g. the marginal likelihood, a customised aggregated score) are calculauted but not mentioned in the paper.
All the evaluation metrics and approaches we considered have been retained in this code in case they interest other users, however we do not vouch for their utility.

Also in this directory are corresponding scripts and folders for the simulation experiments.
The script [copy\_sim\_model\_files.sh](copy_sim_model_files.sh) copies sixteen simulated models (of which only six were discussed in the paper) into `/g/data` again from [sim\_model\_files](sim_model_files).
The script [run\_sims.sh](run_sims.sh) then submits each of these using [run\_and\_evaluate\_model\_parallel.sh](run_and_evaluate_model_parallel.sh) after creating data with [create\_sim\_data.R](create_sim_data.R).
The bulk simulations are run independently of the main scripts (but do require the copying script [copy\_sim\_model\_files.sh](copy_sim_model_files.sh) to be run first) with running and analysis scripts located in [bulk\_sims](bulk_sims).
Finally, the simulations investigating parametric misalignment are explored in the notebook [parametric\_misalignment.ipynb](parametric_misalignment.ipynb) and instructions to explore and experiment are provided therein.
