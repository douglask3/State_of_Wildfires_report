import sys
sys.path.append('../')

import warnings
warnings.filterwarnings('ignore')

import os, fnmatch
from   io     import StringIO
import numpy  as np
import pandas as pd
import csv

import iris
import iris.coords
from iris.coords import DimCoord
import iris.coord_categorisation as cat
import cf_units

import matplotlib.pyplot as plt
import numpy.ma as ma
import cartopy.crs as ccrs
from   libs.plot_maps    import *
from scipy.stats import norm
from scipy.stats import skewnorm
from scipy.stats import lognorm

from pdb import set_trace as browser

import glob

from ConFire import *
from ConFire import ConFire

import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt

input_dir = "isimip3a/driving_data/GSWP3-W5E5/Global/"
output_dir = "isimip3a/outputs/test/"
identfier = "historic_TS"
experiments = ["obsclim/", "counterclim/"]

param_file = "isimip3a/driving_data/GSWP3-W5E5/params-for_sampling/with_ancils-GSWP5.csv-2-1.csv"
periods = fnmatch.filter(os.listdir(input_dir), 'hist*')
periods = sorted(periods, reverse = True)

output = ["controls", "standard", "potential", "sensitivity", "Bootstraps"]
output = ["controls", "standard", "Bootstraps"]

n_posterior_sample = 10
nBootstrap = 10

qs = np.array([1, 5, 10, 25, 50, 75, 90, 95, 99])

ignore_lock = False

################################################################################
## Run period and experiment                                                  ##
################################################################################
def mkDir(dir):
    try: os.mkdir(dir)  
    except: pass


output = np.array(output)
mkDir(output_dir)
run_potential = np.any(np.array(output) == 'potential')
run_sensitivity = np.any(np.array(output) == 'sensitivity')

def makeSampleOutputFile(cube, dir, varName):
    iris.save(cube, dir + varName + '.nc')

def runSample(paramLoc, output_dir, params, input_data):
    print(paramLoc)
    output_dir_sample = output_dir + '/no_' + str(paramLoc) + '/'
    mkDir(output_dir_sample)

    lock_file = output_dir_sample + 'lock.txt'
    if ignore_lock or not os.path.exists(lock_file):
        model = ConFire(input_data, params.loc[paramLoc], 
                        run_potential = run_potential, run_sensitivity = run_sensitivity, nBootstrap = nBootstrap)
        
        makeSampleOutputFile(model.burnt_area_mode, output_dir_sample, 'burnt_area_mode')
        if np.any(np.array(output) == 'controls'):
            makeSampleOutputFile(model.fuel, output_dir_sample, 'fuel')
            makeSampleOutputFile(model.moisture, output_dir_sample, 'moisture')
            makeSampleOutputFile(model.suppression, output_dir_sample, 'suppression')
            makeSampleOutputFile(model.ignitions, output_dir_sample, 'ignitions')
        if np.any(np.array(output) == 'standard'):
            makeSampleOutputFile(model.standard_fuel, output_dir_sample, 'standard_fuel')
            makeSampleOutputFile(model.standard_moisture, output_dir_sample, 'standard_moisture')
            makeSampleOutputFile(model.standard_suppression, output_dir_sample, 'standard_suppression')
            makeSampleOutputFile(model.standard_ignitions, output_dir_sample, 'standard_ignitions')
        if run_potential:
            makeSampleOutputFile(model.potential_fuel, output_dir_sample, 'potential_fuel')
            makeSampleOutputFile(model.potential_moisture, output_dir_sample, 'potential_moisture')
            makeSampleOutputFile(model.potential_suppression, output_dir_sample, 'potential_suppression')
            makeSampleOutputFile(model.potential_ignitions, output_dir_sample, 'potential_ignitions')
        if run_sensitivity:
            makeSampleOutputFile(model.sensitivity_fuel, output_dir_sample, 'sensitivity_fuel')
            makeSampleOutputFile(model.sensitivity_moisture, output_dir_sample, 'sensitivity_moisture')
            makeSampleOutputFile(model.sensitivity_suppression, output_dir_sample, 'sensitivity_suppression')
            makeSampleOutputFile(model.sensitivity_ignitions, output_dir_sample, 'sensitivity_ignitions')

        if np.any(np.array(output) == 'Bootstraps'): 
            output_dir_boots = output_dir_sample  + '/bootstraps/'  
            mkDir(output_dir_boots)
            nboots = model.burnt_area_bootstraps.shape[0]
            for i in range(nboots):
                makeSampleOutputFile(model.burnt_area_bootstraps[i], 
                                     output_dir_boots, 'bootstrap_' + str(i))
        open(lock_file, 'a').close()

def listFiles_recursive(PATH, fname = '', ext = '.nc'):
    return([os.path.join(dp, f) for dp, dn, filenames in os.walk(PATH) for f in filenames if fname in f and os.path.splitext(f)[1] == ext])

def npLogistic(x):
    return(1/(1+np.exp(-x)))

def build_distribution_from_boots(ensembles_dir, output_dir, var, period, output_period = "monthly", transform = None):
    ## setup directories and path
    if var == 'bootstrap_': output_dir_var = output_dir + "burnt_area_full_posterior"
    else: output_dir_var = output_dir + var
    output_dir_var = output_dir_var + '-' + output_period + '/'
    mkDir(output_dir_var)
    outFile = output_dir_var + period + '.nc'

    if ignore_lock or not os.path.exists(outFile):
        ## open and orgaise cube
        files = listFiles_recursive(ensembles_dir, var)
        cubes = iris.load(files)
        
        try:
            mls = np.arange(0.0, len(cubes[1].coord('model_level_number').points))
            for cube in cubes:  
                cube.coord('model_level_number').points = mls
                mls = mls + len(cubes[1].coord('model_level_number').points)
            cubes = cubes.concatenate()[0]
        except:
            for i in range(len(cubes)):
                coord = DimCoord(np.arange(1) + i, "model_level_number")
                cubes[i].add_aux_coord(coord)
            cubes = cubes.merge()[0]
        
        if transform is not None: 
            ifunc = iris.analysis.maths.IFunc(transform, lambda cube: cf_units.Unit('1'))
            cubes = ifunc(cubes)
        if output_period == "annual": 
            iris.coord_categorisation.add_year(cubes, 'time')
            cubes = cubes.aggregated_by('year',  iris.analysis.MEAN)
    
        ## find quanules
        def quantile_time(tstep):
            print(tstep)
            out = cubes[:,[tstep],:,:].collapsed('model_level_number', iris.analysis.PERCENTILE, percent = qs)
            return(out)
        
        cube_quantile = [quantile_time(tstep) for tstep in range(cubes.shape[1])]
        cube_quantile = iris.cube.CubeList(cube_quantile).concatenate()[0]
        iris.save(cube_quantile, outFile)
        

def run_for_period(period, experiment):
    print("experiment: " + experiment + "; period: " + period)
    dir_run = input_dir + period + '/' + experiment
    files = fnmatch.filter(os.listdir(dir_run), '*nc')
    input_data = { file[:-3] : iris.load_cube(dir_run + '/' + file) for file in files }

    params = pd.read_csv(param_file)
    
    output_dir_exp = output_dir + '/' + experiment + '/'
    mkDir(output_dir_exp)
    
    output_dir_period = output_dir_exp + '/' + period + '/'
    mkDir(output_dir_period)

    output_dir_sample = output_dir_period + '/ensembles/'
    mkDir(output_dir_sample)
    
    ngap =int(params.shape[0]/n_posterior_sample)
    sample_nos = range(0, params.shape[0], ngap)
    for i in sample_nos: runSample(i, output_dir_sample, params, input_data)

    build_distribution_from_boots(output_dir_sample, output_dir_exp, "burnt_area_mode", period, "annual")
    build_distribution_from_boots(output_dir_sample, output_dir_exp, "burnt_area_mode", period)
    if np.any(np.array(output) == 'Bootstraps'):
        build_distribution_from_boots(output_dir_sample, output_dir_exp, 'bootstrap_', period, "annual", transform = npLogistic)
        build_distribution_from_boots(output_dir_sample, output_dir_exp, 'bootstrap_', period, transform = npLogistic)
    

for period in periods:
    for experiment in experiments:
        run_for_period(period, experiment)




