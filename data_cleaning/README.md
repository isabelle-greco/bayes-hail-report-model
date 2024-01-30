# Data processing

The main purpose of this folder is to facilitate the pre-processing of the observational data: the radar observations, severe hail reports, and population density estimates.

To do so, first run [copy\_downloaded\_data.sh](copy_downloaded_data.sh) which will copy the raw downloaded data in [downloaded\_data](downloaded_data) to the `g/data` system on Gadi.
This downloaded data includes the observations from the Bureau of Meteorology's [Severe Storms Archive](http://www.bom.gov.au/australia/stormarchive/), the [radar site list](https://www.openradar.io/operational-network), and [population density estimates](https://www.abs.gov.au/statistics/people/population/regional-population/2021#data-download) from the Australian Bureau of Statistics.
Note that the raw ESRI grid and two of the intermediate files is too large to be included here and hence must be downloaded using the link given and the intermediate files created using QGIS.
Alternatively, the authors can provide them upon request.
The radar data is large and not included here, but can be found on Gadi in project `rq0`. 

From this point run [clean\_all\_data.sh](clean_all_data.sh) which processes the reports, grids the radar observations of MESH, and joins these data sets together with the population data before removing the edges of the radar domain.

Note that the [downloaded\_data](downloaded_data) directory contains both the raw [ESRI grids](https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-grid-format.htm) of the population density data and the processed csv.
[QGIS](https://qgis.org/en/site/) was used (offline) to make this conversion as it cannot (currently) be run on Gadi.
The code used to do so is given in [grid\_population/grid\_population.txt](grid_population/grid_population.txt).

This could need only (and should only) be run once.
Be aware of intermediate files potentially causing issues if this code is run multiple times without changing the grid in any way.

