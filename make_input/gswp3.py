import sys
sys.path.append('../libs/')
sys.path.append('libs/')

from read_variable_from_netcdf import *
from plot_maps import *
import os

from   io     import StringIO
import numpy  as np
import math

import matplotlib.pyplot as plt
from pdb import set_trace
from iris.experimental.equalise_cubes import equalise_attributes


def read_variable_from_netcdf_stack(filenames, extensions = '', *args, **kw):
    
    cubes = [read_variable_from_netcdf(file, *args, **kw) for file in filenames]
    cubes = [cube for cube in cubes if cube is not None]
    cubes = iris.cube.CubeList(cubes)
    equalise_attributes(cubes)
    
    cubes = cubes.concatenate_cube()
    return cubes


if __name__=="__main__":
    dir = "/hpc//data/d00/hadea/isimip3a/InputData/climate/atmosphere/obsclim/GSWP3-W5E5/gswp3-w5e5_obsclimfill_"
    
    file_years = ["1991_2000", "2001_2010", "2011_2019"]
    filenames = {"tas": "tas_global_daily_",
                 "tas_range": "tas_range_global_daily_"}

    years = [[2000, 2009], [2010, 2019]]

    def make_variables_for_year_range(year):
        def open_variable(varname):
            set_trace()
            filename = filenames[varname]
            filenames = [filename + ext + '.nc' for ext in extensions]
            tas =  read_variable_from_netcdf_stack(filenames, file_years, dir,
                                                   subset_function = sub_year_range, 
                                                   subset_function_args = {'year_range': year})
        tas = open_variable('tas')
        tas_range = open_variable('tas_range')
        set_trace()

    make_variables_for_year_range(years[0])
