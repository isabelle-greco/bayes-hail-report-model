### Script to check whether a gridded radar file with the correct 
### specifications already exists and gridding should be skipped
### or not.
###
### Author: Isabelle Greco
### Date modified: 2023-07-03

# imports
from gridding_util import check_final_file_exists, find_operational_radars, calculate_bin_edges
import argparse
from datetime import timedelta, datetime
import numpy as np

if __name__ == "__main__":
    # arg parser to get specifications of the file
    parser = argparse.ArgumentParser(
        description = "Verify whether a file exists or not", 
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
    parser.add_argument("-f", "--finaldir", type = str, required = True,
                        help = "Directory in which to save the final gridded file")
    parser.add_argument("-m", "--radarmetadata", type = str, required = True, 
                        help = "File path to a csv with radar metadata")
    args = parser.parse_args()

    # radar of interest
    radar_lat = args.radarlocation[0]
    radar_lon = args.radarlocation[1]

    # time boundaries
    time_start_str = args.startdate
    time_end_str = args.enddate
    time_start = datetime.strptime(time_start_str, "%Y/%m/%d")
    time_end = datetime.strptime(time_end_str, "%Y/%m/%d")
    
    # grid dimensions
    lat_min = args.latitude[0]
    lat_max = args.latitude[1]
    lon_min = args.longitude[0]
    lon_max = args.longitude[1]

    # grid step sizes
    lon_step = args.longitude[2]
    lat_step = args.latitude[2]
    time_step = timedelta(hours = args.timestep)

    # lists of params
    binned = [True, True, True]
    names = ["t", "x", "y"]

    # grid boundaries
    lon_bins = calculate_bin_edges(lon_min, lon_max, lon_step)
    lat_bins = calculate_bin_edges(lat_min, lat_max, lat_step)
    t_bins = calculate_bin_edges(time_start, time_end, time_step)
    grids = [t_bins, lon_bins, lat_bins]

    # radar id
    radar_op = find_operational_radars(time_start = time_start, time_end = time_end,
                                       radar_list_file_path = args.radarmetadata)
    # corresponding closest radar
    closest_idx = np.argmin(
        (radar_op.site_lat - radar_lat)**2 + (radar_op.site_lon - radar_lon)**2
    )
    closest_id = radar_op.iloc[closest_idx, ].id

    # check existence
    check_exist = check_final_file_exists(binned = binned, names = names, grids = grids, 
                                          radar_id = closest_id, data_dir = args.finaldir,
                                          ssa_check = True)

    # pass out
    if check_exist:
        print(1)
    else:
        print(0)
