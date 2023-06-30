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

def read_variable_from_netcdf_stack(filenames, extensions = '', *args, **kw):
    
    cubes = [read_variable_from_netcdf(file, *args, **kw) for file in filenames]
    cubes = [cube for cube in cubes if cube is not None]
    cubes = iris.cube.CubeList(cubes)
    
    cubes = cubes.concatenate_cube()
    return cubes


if __name__=="__main__":
    dir = "/hpc//data/d00/hadea/isimip3a/InputData/climate/atmosphere/obsclim/GSWP3-W5E5/gswp3-w5e5_obsclimfill_"
    
    file_years = ["1991_2000", "2001_2010", "2011_2019"]
    filenames = {"tas": "tas_global_daily_",
                 "tas_range": "tas_range_global_daily_",
                 "pr": "pr_global_daily_"}

    years = [[2000, 2009], [2010, 2019]]

    subset_functions = [constrain_natural_earth]
    subset_function_argss = [{'Country': 'Brazil'}]

    output_dir = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_"
    output_years = ['2000_2009', '2010_2019']

    def make_variables_for_year_range(year, output_year):
        def open_variable(varname):
            filename = filenames[varname]
            filename = [filename + ext + '.nc' for ext in file_years]
            #set_trace()
            sbs_funs = [sub_year_range] + subset_functions 
            sbs_args = [{'year_range': year}] + subset_function_argss
            
            return read_variable_from_netcdf_stack(filename, file_years, dir,
                                                   subset_function = sbs_funs, 
                                                   subset_function_args = sbs_args)
        tas = open_variable('tas')
        tas_range = open_variable('tas_range')

        tas_monthly = tas.aggregated_by(['year', 'month'], iris.analysis.MEAN)

        tas_max = tas.copy()
        tas_max.data  = tas_max.data + 0.5 * tas_range.data
        tas_max_monthly = tas.aggregated_by(['year', 'month'], iris.analysis.MAX)

        def save_ncdf(cube, varname):
            iris.save(cube, output_dir + output_year + '/' + varname +  '.nc')

        save_ncdf(tas_monthly, 'tas_mean')
        save_ncdf(tas_max_monthly, 'tas_max')
        
        
    [make_variables_for_year_range(year, output_year) for year, output_year in \
            zip(years, output_years)]

