### Script to take the relevant MESH file and include the SSA data
### Inspired by stackexchange
###     https://stackoverflow.com/questions/41286677/place-x-y-coordinates-into-bins
###
### Last modified: 2023-07-03
### Author: Isabelle Greco

### Package Imports ###

import pandas as pd
import xarray as xr
from datetime import datetime
import numpy as np
import argparse
import subprocess

### Set up argument parser ###

parser = argparse.ArgumentParser(
    description = ("For a given radar file, will align with hail reports from the "
                   "severe storms archive.")
)
parser.add_argument("-r", "--radarfilepath", type = str, required = True,
                    help = "File path to the radar file of interest.")
parser.add_argument("-s", "--ssafilepath", type = str, required = True,
                    help = "File path to the ssa file to merge with the radar file")
parser.add_argument("-d", "--finaldir", type = str, required = True,
                    help = "Directory to in which to save the merged data")
parser.add_argument("--sep", type = str, required = False, 
                    default = "_", help = "Separator for file names")

# parse args
args = parser.parse_args()

### Read in MESH data ###

mesh_file_path = args.radarfilepath
# remove MESH directory from file name
broken_up_file_path = mesh_file_path.split(f"/")
# need file path later when naming the new file
mesh_file_name = broken_up_file_path[-1] # filename is last part of path

# opening file with mesh data
mesh = xr.open_mfdataset(mesh_file_path, parallel = True)

### Read in SSA data ###

# getting SSA filename
ssa_file_path = args.ssafilepath

# reading in from data directory
#ssa = pd.read_csv(f"{ssa_dir}/{ssa_filename}")
ssa = pd.read_csv(ssa_file_path)
# cleaning ssa times
ssa.date_time = ssa.date_time.apply(lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ"))

### Extract bins from MESH data ###

# extracting the relevant bins
x_bins = mesh.x_bins.values
y_bins = mesh.y_bins.values
time_bins = mesh.time_bins.values

# calculate bin width (assuming equal spacing)
x_width = x_bins[1] - x_bins[0]
y_width = y_bins[1] - y_bins[0]
t_width = time_bins[1] - time_bins[0]

# calculating pandas cut objects
x_cut = pd.cut(
    ssa.longitude, # dimension along which to cut
    # going from bin centres to bin edges
    bins = np.arange(x_bins.min() - 0.5*x_width, x_bins.max() + x_width, x_width),
    labels = x_bins # same labels as SHI data
)
y_cut = pd.cut(
    ssa.latitude,
    bins = np.arange(y_bins.min() - 0.5*y_width, y_bins.max() + y_width, y_width),
    labels = y_bins
)
t_cut = pd.cut(
    ssa.date_time,
    np.arange(time_bins.min() - 0.5*t_width, time_bins.max() + t_width, t_width),
    # can't use time_bins directly as there is missing data
    labels = np.arange(time_bins.min(), time_bins.max() + t_width, t_width) 
)

### Grouping the SSA data ###

# three step process
#    1. Group the ssa_bin by the cuts calculated above
#    2. Grab maximum hail size/corresponding comment and time
#    3. Concatenate these into one pandas table
ssa_grouped = ssa.groupby([x_cut, y_cut, t_cut])
ssa_bin = pd.concat([ssa_grouped.apply(
                         lambda x: x.loc[x.hail_size.idxmax(), 'comments_clean']
                     ).to_frame("comments_clean"),
                     ssa_grouped.hail_size.max(),
                     ssa_grouped.apply(
                         lambda x: x.loc[x.hail_size.idxmax(), 'date_time']
                     ).to_frame("report_date")],
                    axis = 1
                   )

# rename and re order index to make merge easier later
ssa_bin.index.names = ["x_bins", "y_bins", "time_bins"]
ssa_bin = ssa_bin.reorder_levels(["y_bins", "x_bins", "time_bins"])

# convert pandas to xarray
ssa_xr = ssa_bin.to_xarray()

# add attributes to xarray
# dataset-level attributes
ssa_xr.attrs = {"summary" : ("Gridded data drawn from the hail section of the "
                             "Australian severe storms archive. For each grid box "
                             "(with values representing bin centres), the "
                             "corresponding comment and report time are given."),
                "history" : ("Data downloaded in Jan 2023 and gridding performed by "
                             f"Isabelle Greco at {datetime.today().strftime('%Y-%m-%d %H:%M')}."),
                "acknowledgements" : "Supported by the CCRC at UNSW and CLEX",
                "original_data" : ("Data publically available from "
                                   "http://www.bom.gov.au/australia/stormarchive/"),
                "creator_email" : "i.greco@unsw.edu.au",
                "creator_name" : "Isabelle Greco"}

# hail size attributes
ssa_xr.hail_size.attrs = {"units": "cm",
                          "long_name" : "reported hail diameter",
                          "description" : "as per the severe storms archive, the reported hail diameter",
                          "comments" : "likely to be an estimate at best"}

# comments attributes
ssa_xr.comments_clean.attrs = {"long_name" : "comments on reported hail",
                               "description" : "any comments associated with the hail report",
                               "comments" : "may give further information on report's reliability"}

# report date attributes
ssa_xr.report_date.attrs = {"long_name" : "date and time of hail report",
                            "comments" : "likely to be approximate"}

### Merging SSA and SHI data ###

# using xarray merge function as indices are the same
mesh_ssa = xr.merge([mesh.drop_vars(names = ["latitude", "longitude", "isfile"]), ssa_xr], 
                    combine_attrs = "drop_conflicts")

# updating attributes
mesh_ssa.attrs["summary"] = ("Time, latitude, and longitude binning of MESH data from the Australian radar "
                            "network and of hail reports from the severe storms archive.")
old_history = mesh.attrs["history"] 
mesh_ssa.attrs["history"] = (f"{old_history}. Severe Storms Archive data maintained by the Bureau of "
                            f"Meteorology and gridded by Isabelle Greco {datetime.today().strftime('%Y-%m-%d %H:%M')}.")
del mesh_ssa.attrs["acknowledgements"] # duplicated acknowledgement

# file name prepends to radar grid name
merged_file_path = f"{args.finaldir}/{args.sep.join(['ssa', 'variable', 'diameter', 'comment', 'datetime'])}{args.sep}{mesh_file_name}"
# saving the final merged product
mesh_ssa.to_netcdf(path = merged_file_path, engine = "netcdf4")

