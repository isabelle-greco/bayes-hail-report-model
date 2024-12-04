### Script to process the elevation data from Geoscience Australia 
### into a csv that can be read into R.
### Submit to Gadi using save_elevation.sh
###
### Author: Isabelle Greco
### Date modified: 2024-07-20

# imports
import xarray as xr
import numpy as np
import os
import argparse

# create parser
parser = argparse.ArgumentParser(
        description = "Given a file name, extract elevation data around the study area."
    )
parser.add_argument("-f", "--filename", type = str, required = True,
                    help = "Path for file to which to save elevation data as csv")

# parse arguments
args = parser.parse_args()

# define areas of interes
lat = np.arange(-30, -25)
lon = np.arange(150, 155)

# create a grid of lat/lon pairs
lat_grid, lon_grid = np.meshgrid(lat, lon)

# flatten into array for iterating
lat_grid = lat_grid.flatten()
lon_grid = lon_grid.flatten()

# base file name
base_name = "/g/data/rr1/Elevation/NetCDF/1secSRTM_DEMs_v1.0/DEM-S/Elevation_1secSRTM_DEMs_v1.0_DEM-S_Tiles_"

# iterate over lat/lon pairs to obtain file names of interest
fnames = [
    base_name + "e" + str(int(lon)) + "s" + str(int(-1 * lat)) + "dems.nc" for lat, lon in zip(lat_grid, lon_grid)
]

# new list of valid file names
fnames_valid = [f for f in fnames if os.path.isfile(f)]

# open the files into an xarray object
elev = xr.open_mfdataset(fnames_valid)

# coarsen and save data
# coarsening reduces resolution (~750m from ~30m)
elev.elevation.coarsen(
    {"lat" : 25, "lon" : 25}
    ).mean().to_dataframe().to_csv(args.filename)
