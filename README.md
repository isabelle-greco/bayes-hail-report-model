# bayes-hail-report-model

Data and analysis scripts used for the publication:

Greco, I. C., Sherwood, S. C., Raupach, T. H., & Abramowitz G. (submitted 2024). *A Bayesian framework for the probabilistic interpretation of radar observations and severe hailstorm reports*. Weather & Forecasting.

## Navigating this repository

This repository is broken into three directories. 

* [data\_cleaning](data_cleaning): This directory contains both the necessary downloaded data and the cleaning and processing scripts to prepare data for our model. The radar data is very large and so is not included here but its location on the [National Computing Infrastructure (NCI)](https://nci.org.au/) [is provided](https://dapds00.nci.org.au/thredds/catalog/rq0/level_2/catalog.html).

* [model\_evaluation](model_evaluation): In this directory are the necessary files to run and evaluate all the models and simulations discussed in the paper (and some that are not!).

* [figures](figures): As it says on the tin, the code to generate all (supplementary) figures is included in this folder.

## Software and packages

The code in this repository was designed to be run on [Gadi](https://nci.org.au/our-systems/hpc-systems) to easily access the [radar data](https://dapds00.nci.org.au/thredds/catalog/rq0/level_2/catalog.html) and to obtain adequate computing power to efficiently run the models.
Python used the ARC Centre of Excellence for Climate Extremes' conda installation in project [hh5](http://climate-cms.wikis.unsw.edu.au/Conda) and we are immensely grateful for the support of the Computational Modelling Support team. 
All python packages used are included in this installation .
R code used the [NCI-data-analysis](https://opus.nci.org.au/pages/viewpage.action?pageId=134742126) installation in project `dk92`.
Several additional packages were employed:

* [ozmaps](https://mdsumner.github.io/ozmaps/),
* [ggspatial](https://paleolimbot.github.io/ggspatial/),
* [latex2exp](https://cran.r-project.org/web/packages/latex2exp/index.html),
* [posterior](https://mc-stan.org/posterior/reference/posterior-package.html),
* [janitor](https://www.rdocumentation.org/packages/janitor/versions/2.2.0),
* [tabularaster](https://cran.r-project.org/web/packages/tabularaster/index.html),
* [ggridges](https://cran.r-project.org/web/packages/ggridges/index.html)
* [bayesplot](https://mc-stan.org/bayesplot/),
* [twinning](https://cran.r-project.org/web/packages/twinning/index.html),
* [poibin](https://cran.r-project.org/web/packages/poibin/index.html),
* [bigstatsr](https://privefl.github.io/bigstatsr/),
* [bridgesampling](https://cran.r-project.org/web/packages/bridgesampling/index.html).