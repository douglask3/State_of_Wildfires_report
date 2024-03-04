import sys
sys.path.append('../libs/')
sys.path.append('libs/')

from read_variable_from_netcdf import *
from combine_path_and_make_dir import *
from plot_maps import *
import os
from pathlib import Path

from   io     import StringIO
import numpy  as np
import math

import matplotlib.pyplot as plt
from pdb import set_trace

def read_variable_from_netcdf_stack(filenames, example_cube = None, 
                                    *args, **kw):
    
    cubes = [read_variable_from_netcdf(file, *args, **kw) for file in filenames]
    cubes = [cube for cube in cubes if cube is not None]
    cubes = iris.cube.CubeList(cubes)
    cubes = cubes.concatenate_cube()
    
    if example_cube is not None:
        example_cube = iris.load_cube(example_cube)
        
        cubes = cubes.regrid(example_cube, iris.analysis.Linear())
        cubes.data = np.ma.masked_greater(cubes.data, 9E9)
        cubes.data[cubes.data > 9E9] = np.nan
    return cubes


def make_variables_for_year_range(year, output_year):
    def test_if_process(var): return any(i == var for i in process)

    def open_variable(varname, plusMinusYr = False):
        filename = filenames[varname]

        yeari = year.copy()
        if plusMinusYr:
            yeari[0] = yeari[0] - 1
            yeari[1] = yeari[1] + 1
        
        file_years_ranges = [(int(pair.split('_')[0]), int(pair.split('_')[1])) \
                             for pair in file_years]
        # Find overlapping elements
        overlapping_years = []
        for i, range_pair in enumerate(file_years_ranges):
            if yeari[1] >= range_pair[0] and yeari[0] <= range_pair[1]:
                overlapping_years.append(file_years[i])
        #set_trace()
        filename = [filename + ext + '.nc' for ext in overlapping_years]
        
        sbs_funs = [sub_year_range] + subset_functions 
        sbs_args = [{'year_range': yeari}] + subset_function_argss
        
        out =  read_variable_from_netcdf_stack(filename, example_cube, dir,
                                               subset_function = sbs_funs, 
                                               subset_function_args = sbs_args)
        
        return out

    def monthly_mean(cube, fun = iris.analysis.MEAN):
        return cube.aggregated_by(['year', 'month'], fun)

    def save_ncdf(cube, varname):        
        out_dir = output_dir + '/' + \
                    subset_function_argss[0][next(iter(subset_function_argss[0]))] + '/' + \
                    dataset_name + '/period_' + output_year + '/'
        
        if not os.path.exists(out_dir): Path(out_dir).mkdir(parents=True)
        iris.save(cube, out_dir + '/' + varname +  '.nc')

    def standard_Monthly_mean(var, fun):
        if test_if_process(var):
            dat = open_variable(var)
            mdat = monthly_mean(dat, fun)
            save_ncdf(mdat, var + '_mean')    

    for var, fun in zip(process_standard, process_function):
        standard_Monthly_mean(var, fun) 
   
    if test_if_process('tas') or test_if_process('vpd'):
        tas = open_variable('tas')
        tas_range = open_variable('tas_range')
        
        tas_max = tas.copy()
        tas_max.data  = tas_max.data + 0.5 * tas_range.data

        if test_if_process('vpd'):
            def SVP(temp):
                svp = temp.copy()
                
                svp = svp - 273.16
                svp.data =  610.78 * np.exp(svp.data / (svp.data +237.3) * 17.2694)
                return svp

            rh = open_variable("hurs")
            svp = SVP(tas_max)
            vpd = svp*(1.0-rh*0.01)
            
            vpd.data[vpd.data<0.0] = 0.0            
            vpd_mean_monthly = monthly_mean(vpd)   
            vpd_max_monthly = monthly_mean(vpd, iris.analysis.MAX)
            
            save_ncdf(vpd_max_monthly, 'vpd_max')
            save_ncdf(vpd_mean_monthly, 'vpd_mean')
            
        if test_if_process('tas'):
            
            tas_monthly = monthly_mean(tas)
            tas_max_monthly = tas.aggregated_by(['year', 'month'], iris.analysis.MAX)

            save_ncdf(tas_monthly, 'tas_mean')
            save_ncdf(tas_max_monthly, 'tas_max')
            
    if test_if_process('pr'):
        pr = open_variable('pr', True)
        pr_mean = monthly_mean(sub_year_range(pr, year))
        
        dry_days = pr.copy()
        dry_days.data[dry_days.data < (1.0/86400.0)] = 0.0
        dry_days.data[dry_days.data > 0.0] = 1.0
        dry_days.data = 1.0 - dry_days.data
        dry_days.units = '1'

        def cummDry(data):
            if np.isnan(data[0]): return data
            try: 
                if data.mask[0]: return data
            except:
                pass
                
            out = data.copy()
            out[:] = 0.0
            for i in range(1, len(out)):
                if data[i] > 0.0: out[i] = out[i-1] + 1
            return out

        consec_dry = dry_days.copy()
        consec_dry.data = np.apply_along_axis(cummDry, 0, dry_days.data)
        
        dry_days_mean = monthly_mean(sub_year_range(dry_days, year))
        consec_dry_mean = monthly_mean(sub_year_range(consec_dry, year), iris.analysis.MAX)
            
        save_ncdf(pr_mean, 'pr_mean')
        save_ncdf(dry_days_mean, 'dry_days')
        save_ncdf(consec_dry_mean, 'consec_dry_mean')

filenames = {"tas": "tas_global_daily_",
             "tas_range": "tas_range_global_daily_",
             "pr": "pr_global_daily_",
             "prsn": "ps_global_daily_",
             "hurs": "hurs_global_daily_",
             "huss": "hurs_global_daily_",
             "sfcwind": "sfcwind_global_daily_",
             "ps": "ps_global_daily_"}

process = ['vpd', 'tas', 'tas_range', 'pr']
process_standard = ['prsn', "hurs", "hurs", "huss", "huss", "sfcwind"]
process_function = [iris.analysis.MEAN, 
                    iris.analysis.MEAN, iris.analysis.MAX,
                    iris.analysis.MEAN, iris.analysis.MAX,
                    iris.analysis.MAX]

example_cube = None#'../../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/GFED4.1s_Burned_Fraction.nc'

subset_functions = [ar6_region]
subset_function_argss = [{'region_code': 'NWN'}]
output_dir = "../data/data/driving_data/"

if __name__=="__main__":
    dir = "/hpc//data/d00/hadea/isimip3a/InputData/climate/atmosphere/obsclim/GSWP3-W5E5/gswp3-w5e5_obsclimfill_"
    
    file_years = ["1901_1910", "1911_1920", "1991_2000", "2001_2010", "2011_2019"]
    
    years = [[2010, 2012], [1901, 1920], [2000, 2019], [2002, 2019]]
    dataset_name = 'isimp3a/GSWP3-W5E5'
    
    output_years = ['2010_2012', '1901_1920', '2000_2019', '2002_2019']

    [make_variables_for_year_range(year, output_year) for year, output_year in \
            zip(years, output_years)]



