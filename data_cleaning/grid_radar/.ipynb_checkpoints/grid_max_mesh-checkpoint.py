### Script to read in MESH fields (doi: 10.25914/JJWZ-0F13) and grid them onto a user
### specified grid, taking the maximum for each cell. Designed to run using Gadi's
### parallel processing capabilities.
### 
### Last modified: 2023-07-03
### Author: Isabelle Greco

### Package Imports ###

import gridding_util # custom utility functions
from datetime import datetime, timedelta # handling date times
import numpy as np
import netCDF4 # writing file
import argparse # argument parser
from multiprocessing import Pool # for the parallel processing
import xarray as xr 

if __name__ == "__main__":

    ### Create argument parser ###
    
    parser = argparse.ArgumentParser(
        description = "Given a radar location and grid parameters, create maximum MESH grid."
    )
    parser.add_argument("-s", "--startdate", type = str, 
                        help = "Starting date of gridding (YYYY/MM/DD)", required = True)
    parser.add_argument("-e", "--enddate", type = str, 
                        help = "Ending date of gridding (YYYY/MM/DD)", required = True)
    parser.add_argument("-r", "--radarlocation", type = float, required = True, nargs = 2,
                        help = "Location near which to seek the radar (lat, lon)")
    parser.add_argument("--latitude", type = float, required = True, nargs = 3,
                        help = "Latitude grid parameters (min, max, step)")
    parser.add_argument("--longitude", type = float, required = True, nargs = 3,
                        help = "Longitude grid parameters (min, max, step)")
    parser.add_argument("-t", "--timestep", type = float, required = True,
                        help = "Time step of the grid in hours")
    parser.add_argument("-m", "--radarmetadata", type = str, required = True, 
                        help = "File path to a csv with radar metadata")
    parser.add_argument("-i", "--intermediatedir", type = str, required = True,
                        help = "Directory in which to save intermediate gridding files")
    parser.add_argument("-f", "--finaldir", type = str, required = True,
                        help = "Directory in which to save the final gridded file")
    args = parser.parse_args()

    ### User-defined parameters ###

    # radar of interest
    radar_lat = args.radarlocation[0]
    radar_lon = args.radarlocation[1]

    # time boundaries
    time_start_str = args.startdate
    time_end_str = args.enddate

    # grid dimensions
    lat_min = args.latitude[0]
    lat_max = args.latitude[1]
    lon_min = args.longitude[0]
    lon_max = args.longitude[1]

    # grid step sizes
    lon_step = args.longitude[2]
    lat_step = args.latitude[2]
    time_step = timedelta(hours = args.timestep)

    # assert valid interval and step size
    assert lat_min < lat_max
    assert lon_min < lon_max
    assert lat_step < lat_max - lat_min
    assert lon_step < lon_max - lon_min

    ### Formatting Dates ###

    # input date format
    input_date_format = "%Y/%m/%d"

    # grid dates
    time_start = datetime.strptime(time_start_str, input_date_format)
    time_end = datetime.strptime(time_end_str, input_date_format)

    # assert valid time interval and valid step size
    assert time_start < time_end
    assert time_step < time_end - time_start

    ### Finding Closest Relevant Radar ###

    # operational radars
    radar_op = gridding_util.find_operational_radars(time_start = time_start, time_end = time_end,
                                                     radar_list_file_path = args.radarmetadata)

    # corresponding closest radar
    closest_idx = np.argmin(
        (radar_op.site_lat - radar_lat)**2 + (radar_op.site_lon - radar_lon)**2
    )
    closest_id = radar_op.iloc[closest_idx, ].id

    # info about closest radar
    radar_info = radar_op[radar_op.id == closest_id]

    # assert radar in grid (i.e. radar actually near point of interest)
    assert np.all(radar_info.site_lat > lat_min)
    assert np.all(radar_info.site_lat < lat_max)
    assert np.all(radar_info.site_lon > lon_min)
    assert np.all(radar_info.site_lon < lon_max)


    # find available files
    avail_files = gridding_util.find_available_radar_files(
        closest_id, time_start, time_end
    )

    ### Gridding and saving data ###

    # definine grid edges
    lon_bins = gridding_util.calculate_bin_edges(lon_min, lon_max, lon_step)
    lat_bins = gridding_util.calculate_bin_edges(lat_min, lat_max, lat_step)
    t_bins = gridding_util.calculate_bin_edges(time_start, time_end, time_step)

    # check file exists before doing all the processing
    f_name = gridding_util.check_final_file_exists(binned = [True, True, True], 
                                                  names = ["t", "x", "y"], 
                                                  grids = [t_bins, lon_bins, lat_bins], 
                                                  radar_id = closest_id,
                                                  data_dir = args.finaldir) 
    
    def process_file(file):
        gridding_util.grid_one_file(file, lon_bins, lat_bins, args.timestep, 
                                    closest_id, args.intermediatedir)
    
    # start parallel processing
    if f_name:
        print(f_name)
    else:
        # creates local cluster on gadi
        with Pool() as pool:    
            pool.map(process_file, avail_files)
        
        # read in as one object
        all_data = xr.open_mfdataset(f"{args.intermediatedir}/*", parallel = True)
        # create file name
        grid_fname = gridding_util.create_file_name(binned = [True, True, True], # all axes binned in the grid_one_file function
                                      names = ["t", "x", "y"], # default names
                                      grids = [t_bins, lon_bins, lat_bins], 
                                      radar_id = closest_id)
        # make full file path
        full_fpath = f"{args.finaldir}/{grid_fname}.nc" # should be equal to f_name
        # to data frame and then to grid
        all_data.to_netcdf(path = full_fpath, engine = "netcdf4")
        # returns the path
        print(full_fpath)
