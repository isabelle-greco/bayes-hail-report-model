### Utility functions written to support gridding_max_mesh.py which take care
### of some of the more annoying code.
###
### Last modified: 2023-07-03
### Author: Isabelle Greco


### Package Imports ###

import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
import subprocess
import dask
from pathlib import Path
from skimage import morphology

### Classes ###

class Error(Exception):
    """
    Base class for other exceptions
    """
    pass


class NoAvailableRadarError(Error):
    """
    Raised when there is no available radar on the day(s) of interest.
    """
    pass

### Functions ###

def _try_str_to_date(strg : str):
    """
    Attempts to coerce a string from a specified date format to a date.
    If unavailable, returns today's date.

    Input:
        strg    string for conversion

    Output:
        _       datetime object
    """
    # attempt conversion
    try:
        return datetime.strptime(strg, "%d/%m/%Y")
    except:
        # returns today's date if cannot convert
        return datetime.today().date()


def _bin_edge_to_centre(grid : np.ndarray):
    """
    Given the grid edges returns the centre of the grid points.
    Assumes that the grid is equally spaced.

    Input:
        grid       1 x n numpy array with bin edges

    Output:
        centres    1 x (n-1) numpy array giving the centres
    """
    # calculate width of grid box under equality assumption
    bin_width = (grid.max() - grid.min())/(len(grid) - 1)

    # check that this width works for the first two elements
    # if this fails then the gird is not equally spaced and
    # this function won't work
    assert bin_width == grid[1] - grid[0]

    # adds half bin width to left edges
    centres = grid[:-1] + 0.5 * bin_width

    # check correct length
    assert len(centres) == len(grid) - 1

    # check that the centres are contained in the grid points
    assert centres.min() > grid.min()
    assert centres.max() < grid.max()
    assert grid[0] < centres[0] and centres[0] < grid[1]

    return centres


def _format_string_file_name(x_grid : np.ndarray, name : str, sep : str = "_",
                            decimal_point : str = "."):
    """
    Given a grid, formats it into an identifiable portion of a file name.
    Assumes an equally spaced grid

    Inputs:
        x_grid           the grid that we are describing
        name             label of the input
        sep              character(s) used to separate portions of input
        bullet_points    character(s) used to represent or replace decimal points

    Output:
        _                name, min, max, and step separated by sep and formated
                         to two decimal places
    """
    if x_grid.dtype.type == np.datetime64:
        return (f"{name}{sep}{np.datetime_as_string(x_grid.min(), 'D')}{sep}"
                f"{np.datetime_as_string(x_grid.max(), 'D')}{sep}"
                f"{np.timedelta64(x_grid[1] - x_grid[0], 'h').astype(int)}").replace(
                    ".", decimal_point)
    elif x_grid.dtype.type == np.int64 or x_grid.dtype.type == np.float64:
        return (f"{name}{sep}{x_grid.min():.2f}{sep}{x_grid.max():.2f}{sep}"
                f"{x_grid[1] - x_grid[0]:.2f}").replace(
                    ".", decimal_point)
    else:
        raise TypeError("Cannot format grids of this type")


def create_file_name(binned : list, names : list, grids : list, radar_id : int,
                      sep : str = "_", sep_within_grid : str = "_",
                      decimal_point : str = ".", var : str = "mesh"):
    """
    Given the binning status of some variables, their names, the corresponding grids
    and the seperator, returns a file name.

    Inputs:
        binned             list of True/False indicators of whether the corresponding
                           variable has been binned
        names              list of strings given short names of variables
        grids              list np.ndarrays representing the relevant girds
        radar_id           id of the radar whose data is being gridded
        sep                string to use to separate each component
        sep_within_grid    seperator to use within information for each grid
        decimal_point      how to replace decimal points
        var                string name of variable being gridded

    Output:
        _                  string describing the bounds of the data binned variables
    """
    # begin file with name of radar and name of variable
    comps = ["radar", str(radar_id), "variable", var]

    # iterate throug binning status to determine prefixes
    for yes_binned, name in zip(binned, names):
        if yes_binned:
            comps.append(name)

    # if we created a prefix then add the binned indicated
    if any(binned):
        comps.append("binned")

    # use names and grids to get file names reflecting grid bounds
    for name, grid in zip(names, grids):
        comps.append(_format_string_file_name(grid, name, sep = sep_within_grid, 
                                              decimal_point = decimal_point))

    # join the components together separated by the given separator
    return sep.join(comps)


def _despeckle(x):
    """
    As per Josh Soderholm's advice, despeckle the data by removing groups
    of contiguous pixles smaller than 9 pixels in size.
    
    Inputs:
        x   np.ndarray, the input array of booleans where we have readings.
    
    Outputs:
        _   np.ndarray (same size), with small objects removed.
    """
    return morphology.remove_small_objects(x, min_size = 9)


def check_final_file_exists(binned : list, names : list, grids : list, radar_id : int,
                            data_dir : str, sep : str = "_", sep_within_grid : str = "_",
                            decimal_point : str = ".", var : str = "mesh", ssa_check : bool = False):
    """
    Given the binning parameters, checks if the file is already available.

    Inputs:
        binned             list of True/False indicators of whether the corresponding
                           variable has been binned
        names              list of strings given short names of variables
        grids              list np.ndarrays representing the relevant girds
        radar_id           id of the radar whose data is being gridded
        data_dir           directory for gridded files
        sep                string to use to separate each component
        sep_within_grid    seperator to use within information for each grid
        decimal_point      how to replace decimal points
        var                string name of variable being gridded
        ssa_check          check for a file starting with the ssa information

    Output:
        _                  the file path if the file exists, else None
    """
    file_base = create_file_name(binned = binned, names = names, grids = grids, 
                                  radar_id = radar_id, sep = sep, 
                                  sep_within_grid = sep_within_grid, 
                                  decimal_point = decimal_point, var = var)
    if ssa_check:
        # prepend the ssa information to file_base
        file_base = f"{sep.join(['ssa', 'variable', 'diameter', 'comment', 'datetime'])}{sep}{file_base}"
    file_path = f"{data_dir}/{file_base}.nc"
    if Path(file_path).is_file():
        return file_path
    else:
        return None


def calculate_bin_edges(x_min : float, x_max : float, x_step : float):
    """
    Given the min, max, and step of a grid axis return bin edges

    Input:
        x_min     float, minimum value of grid (i.e. left-most edge)
        x_max     float, maximum value of grid (i.e. right-most edge)
        x_step    float, width of each grid cell

    Output:
        _         numpy array with bin edges
    """
    return np.arange(start = x_min, stop = x_max + x_step, step = x_step)


def find_operational_radars(time_start : datetime, radar_list_file_path : str,
                            time_end : datetime = datetime.today()):
    """
    Returns a list of information about all operational Australian radars in the time
    period of interest.

    Input:
        time_start              datetime, start of time period of interest
        radar_list_file_path    str, file path of the csv with the radar information
        time_end                datetime, end of time period of interest
                                defaults to the current date and time

    Outpu:
        radar_op                pandas csv with the radar information
    """
    # check that inputs are valid
    assert time_start < time_end
    assert radar_list_file_path[-4:] == ".csv"

    # read radar list
    radar_list = pd.read_csv(radar_list_file_path)

    # cleaning
    radar_list_clean = radar_list.assign(
        start_date = lambda dataframe: dataframe.postchange_start.map(_try_str_to_date)
    ).assign(
        end_date = lambda dataframe: dataframe.prechange_end.map(_try_str_to_date)
    ).drop(
        ["postchange_start", "prechange_end"], axis = 1
    )

    # stopped operating during period of interest
    cond1 = radar_list_clean.end_date > time_start

    # started operating at some point
    cond2 = np.logical_not(np.isnat(radar_list_clean.start_date))

    # started operating before time period ended
    cond3 = radar_list_clean.start_date < time_end

    # full condition:
    #    definitely started AND
    #    (stopped after the time period of interest started OR
    #      started before the time period stopped)
    full_cond = np.logical_and(cond2, np.logical_and(cond1, cond3))

    # find operational radars
    radar_op = radar_list_clean[full_cond]

    # return operational radar information
    return radar_op


def find_available_radar_files(radar_id : int,
                               time_start : datetime,
                               time_end : datetime ,
                               base_path : str = "/g/data/rq0/level_2",
                               file_date_format : str = '%Y%m%d'
                              ):
    """
    Compares the desired days with the data directory to return a list of available files.

    Input:
        radar_id            int, ID of the radar
        time_start          datetime, start of time period of interest
        time_end            datetime, end of time period of interest
        base_path           str, the data directory for the radar data
                            defaults to the level 2 radar data (project rq0)
        file_date_format    str, the date format using for naming the radar files
                            defaults to %Y%m%d (e.g. 220609)

    Output:
        file_location       list of full file name strings.
                            if there are no available days, raises an error.

    """
    # ls the directory to find available files
    proc = subprocess.Popen(["ls", f"{base_path}/{radar_id}/SHI"], stdout = subprocess.PIPE)

    # cleaning the stdout
    avail = proc.stdout.read().decode().strip("\n").split("\n")

    # check how many days are requested
    one_day = time_start.strftime(file_date_format) == time_end.strftime(file_date_format)

    if one_day:
        # only one day requested
        file_name = f"{radar_id}_{time_start.strftime(file_date_format)}_shi.nc"

        if not file_name in avail:
            raise NoAvailableRadarError((f"Radar {radar_id} has no files available on "
                                         f"{time_start.strftime('%Y-%m-%d')}"))


        # if available, format the full file puat
        file_location  = ['/'.join([base_path, str(radar_id), 'SHI', file_name])]
    else:
        # multiple days requested
        days = np.arange(start = time_start.date(),
                         stop = time_end.date() + timedelta(days = 1),
                         step = timedelta(days = 1)).astype(datetime)

        # turning the datetimes into appropriate strings
        days_str = [x.strftime(file_date_format) for x in days]

        # formatting into file name
        file_names = [f"{radar_id}_{day}_shi.nc" for day in days_str]

        # intersection of the available files and the desired file names
        files_avail = list(set(avail).intersection(file_names))

        # check if there are any available
        if not files_avail:
            tstart = time_start.strftime('%Y-%m-%d')
            tend = time_end.strftime('%Y-%m-%d')
            raise NoAvailableRadarError(
                f"Radar {radar_id} has no files available between {tstart} and {tend}"
            )

        # create list of full file paths
        file_location = [
            '/'.join([base_path, str(radar_id), 'SHI', name]) for name in files_avail
        ]

    # return file location list or string
    return file_location


def grid_one_file(file_name : str, x_grid : np.ndarray, y_grid : np.ndarray, t_int : float, 
                  radar_id : int, int_data_dir : str, x_name : str = "x",
                  y_name : str = "y", t_name : str = "t"):
    """
    For a given grid, radar ID, and filename will despeckle and calculate maximum MESH and save the result.

    Input:
        file            str, file to query for MESH
        x_grid          ndarray, bin edges for the x (longitude) direction
        y_grid          ndarray, bin edges for the y (latitude) direction
        t_int           float, time interval for the temporal grid (hours)
        radar_id        int, id of the radar being used for labelling
        int_data_dir    str, directory to which intermediate files are saved
        x_name          str, name to use for x in the file name
        y_name          str, name to use for y in the file name
        t_name          str, name to use for time in the file name
    
    Output:
        no output
    """
    # TODO: ensure integration now that file is in pieces and NO OUTPUT PROVIDED!
    # TODO: int_data_dir needs to be the daily_files folder

    # set binning status
    t_binned = True
    x_binned = True
    y_binned = True
    # get file ID (the date)
    f_id = file_name.split("/")[-1].split("_")[1]
    # create time grid from the date
    t_grid = np.arange(datetime.strptime(f_id, "%Y%m%d"), 
                       datetime.strptime(f_id, "%Y%m%d") + timedelta(hours = 24 + t_int), 
                       step = timedelta(hours = t_int))
    # create file name
    f_name = create_file_name([t_binned, x_binned, y_binned], 
                              [t_name, x_name, y_name], 
                              [t_grid, x_grid, y_grid], 
                              radar_id) # will be same for each
    path_y = f"{int_data_dir}/{f_name}.nc"

    # check if file already exists
    if not Path(path_y).is_file():
        # load all space and chunk via time (10/file) as this is where we get biggest reduction
        with xr.open_dataset(file_name) as file:
            with dask.config.set(**{"array.slicing.split_large_chunks": True}):
                # replace projection with latitude and longitude
                file.coords["x"] = file.longitude.sel(y = 0, method = "nearest").values
                file.coords["y"] = file.latitude.sel(x = 0, method = "nearest").values
                
                # speckle filtering 
                file["valid"] = (file.shi >= 0).astype('bool')
                file["valid_despeckle"] = xr.apply_ufunc(_despeckle, file.valid, 
                                                         dask = "parallelized")
 
                # adding MESH where valid as per Witt et al. (1998)
                # where SHI unavailable leave as NA (default behaviour)
                file["mesh"] = 2.54 * (file.shi.where(file.valid_despeckle) ** 0.5)
                
                ### time gridding ### 
                
                # do the gridding and take maximum at each step
                grid_time = file.groupby_bins(
                    "time", bins = t_grid
                ).max(keep_attrs = True)

                # calculate + set centre labels
                t_label = _bin_edge_to_centre(t_grid)
                grid_time.coords["time_bins"] = t_label

                # edit high level metadata
                grid_time.attrs["summary"] = ("Time-binned Level 2 dataset from the Australian radar "
                                              "network:  maximum expected size of hail at the surface "
                                              "using Witt et al. 1998")
                old_history = grid_time.attrs["history"]
                # updating history
                grid_time.attrs["history"] = ("time binning by Isabelle Greco at "
                                              f"{datetime.today().strftime('%Y-%m-%d %H:%M')} based on "
                                              f"data {old_history}")
                grid_time.attrs["acknowledgement"] = "Supported by the CCRC at UNSW and CLEX"
                grid_time.attrs["institution"] = "UNSW"
                grid_time.attrs["creator_email"] = "i.greco@unsw.edu.au"
                grid_time.attrs["creator_name"] = "Isabelle Greco"
                grid_time.attrs["references"] = "doi: 10.25914/JJWZ-0F13"
                del grid_time.attrs["creator_url"]

                # adding meta data to mesh
                grid_time.mesh.attrs["long_name"] = "maximum expected size of hail"
                grid_time.mesh.attrs["units"] = "millimetres"
                grid_time.mesh.attrs["description"] = ("De-speckled MESH as developped by Witt et al. 1998" ,
                                                       "doi:10.1175/1520-0434(1998)013<0286:AEHDAF>2.0.CO;2")
                grid_time.mesh.attrs["valid_min"] = 0
                grid_time.mesh.attrs["profile_source"] = grid_time.shi.attrs["profile_source"]

                # assign meta data
                grid_time.longitude.attrs["long_name"] = "longitude_degrees_east" # correcting meta data
                bin_width = (t_grid[1] - t_grid[0]).astype('timedelta64[h]').astype(int) # assumes equal grid
                grid_time.time_bins.attrs["long_name"] = f"time in the centre of {bin_width}-hour bins"

                # removing unnecessary variables
                grid_time = grid_time.drop(["valid", "valid_despeckle", "shi"])
                
                ### zonal gridding ### 
                
                grid_x = grid_time.groupby_bins(
                    "x", bins = x_grid
                ).max(keep_attrs = True)

                # calculate + set centre lables
                x_label = _bin_edge_to_centre(x_grid)
                grid_x.coords["x_bins"] = x_label

                # edit overall metadata
                grid_x.attrs["summary"] = ("Time- and longitude-binned Level 2 dataset from the "
                                           "Australian radar network:  maximum expected hail size "
                                           "(MESH) at the surface using Witt et al. 1998")
                old_history = grid_x.attrs["history"]
                grid_x.attrs["history"] = ("longitude binning by Isabelle Greco at "
                                           f"{datetime.today().strftime('%Y-%m-%d %H:%M')} based on {old_history}")

                # assign new metadata to bins
                grid_x.x_bins.attrs["long_name"] = f"centre of {x_grid[1] - x_grid[0]}-degree bins"
                grid_x.x_bins.attrs["units"] = "degrees east"

                ### meridional gridding ### 
                
                gridded = grid_x.groupby_bins(
                    "y", bins = y_grid
                ).max(keep_attrs = True)

                # calculate + set centre labels
                y_label = _bin_edge_to_centre(y_grid)
                gridded.coords["y_bins"] = y_label

                # edit overall meta data
                gridded.attrs["summary"] = ("Time-, latitude-, and longitude-binned Level 2 dataset "
                                            "from the Australian radar network: maximum expected hail "
                                            "size (MESH) at the surface using Witt et al. 1998")
                old_history = gridded.attrs["history"]
                gridded.attrs["history"] = ("latitude binning by Isabelle Greco at "
                                            f"{datetime.today().strftime('%Y-%m-%d %H:%M')} "
                                            f"based on {old_history}")

                # assign new metadata
                gridded.y_bins.attrs["long_name"] = f"centre of {y_grid[1] - y_grid[0]}-degree bins"
                gridded.y_bins.attrs["units"] = "degrees north"

                # write out
                gridded.to_netcdf(path = path_y, engine = "netcdf4")
            
