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
output_dir = "isimip3a/outputs/test-alldat-noMT-DE/"
identfier = "historic_TS"
experiments = ["counterclim/", "obsclim/"]

param_file = "isimip3a/driving_data/GSWP3-W5E5/params-for_sampling/with_ancils-GSWP5.csv-2-1.csv"
param_file = "isimip3a/driving_data/GSWP3-W5E5/params-for_sampling/with_ancils_alldats-newVPD-pow-GSWP5.csv-2-0.1.csv"

param_file = "isimip3a/driving_data/GSWP3-W5E5/params-for_sampling/with_ancils_alldats-newVPD-pow-normal-GSWP5.csv-2-0.2.csv"
param_file = "isimip3a/driving_data/GSWP3-W5E5/params-for_sampling/with_ancils_alldats-newVPD-pow-normal-noMoistT-DE-GSWP5.csv-2-0.2.csv"
periods = fnmatch.filter(os.listdir(input_dir), 'hist*')
periods = sorted(periods, reverse = True)

output = ["controls", "standard", "potential", "sensitivity", "Bootstraps"]
output = ["controls", "standard"]#, "Bootstraps"]

n_posterior_sample = 10
nBootstrap = 10

qs = np.array([1, 5, 10, 25, 50, 75, 90, 95, 99])

ignore_lock_samples = False
ignore_lock_summery = False
################################################################################
## Run period and experiment                                                  ##
################################################################################
def mkDir(dir):
    try: os.mkdir(dir)  
    except: pass


output = np.array(output)
mkDir(output_dir)
run_controls = np.any(np.array(output) == 'controls')
run_standard = np.any(np.array(output) == 'standard')
run_potential = np.any(np.array(output) == 'potential')
run_sensitivity = np.any(np.array(output) == 'sensitivity')
run_Bootstraps = np.any(np.array(output) == 'Bootstraps')

def makeSampleOutputFile(cube, varName, dir, period, ens_no):
    var_dir = dir + '/' + varName + '/'
    mkDir(var_dir)
    ens_dir = var_dir + 'ensemble_' + str(ens_no) + '/'
    mkDir(ens_dir)
    iris.save(cube, ens_dir + period + '.nc')

def runSample(paramLoc, output_dir, params, input_data, period):
    print(paramLoc)
    #output_dir_sample = output_dir + '/no_' + str(paramLoc) + '/'
    #mkDir(output_dir_sample)
    def makeSampleOutputFileLocal(cube, varname, ens_no = paramLoc):
        makeSampleOutputFile(cube, varname, output_dir, period, ens_no)
    
    output_dir_lock = output_dir + '/lockFiles/'
    mkDir(output_dir_lock)
    lock_file = output_dir_lock + period + str(paramLoc) + '_lock.txt'
    if ignore_lock_samples or not os.path.exists(lock_file):
        model = ConFire(input_data, params.loc[paramLoc], 
                        run_potential = run_potential, run_sensitivity = run_sensitivity, nBootstrap = nBootstrap)
        
        makeSampleOutputFileLocal(model.burnt_area_mode, 'burnt_area_mode')
        if run_controls:
            makeSampleOutputFileLocal(model.fuel, 'fuel')
            makeSampleOutputFileLocal(model.moisture, 'moisture')
            makeSampleOutputFileLocal(model.suppression, 'suppression')
            makeSampleOutputFileLocal(model.ignitions, 'ignitions')
        if run_standard:
            makeSampleOutputFileLocal(model.standard_fuel, 'standard_fuel')
            makeSampleOutputFileLocal(model.standard_moisture, 'standard_moisture')
            makeSampleOutputFileLocal(model.standard_suppression, 'standard_suppression')
            makeSampleOutputFileLocal(model.standard_ignitions, 'standard_ignitions')
        if run_potential:
            makeSampleOutputFileLocal(model.potential_fuel, output_dir_sample, 'potential_fuel')
            makeSampleOutputFileLocal(model.potential_moisture, 'potential_moisture')
            makeSampleOutputFileLocal(model.potential_suppression, 'potential_suppression')
            makeSampleOutputFileLocal(model.potential_ignitions, 'potential_ignitions')
        if run_sensitivity:
            makeSampleOutputFileLocal(model.sensitivity_fuel, 'sensitivity_fuel')
            makeSampleOutputFileLocal(model.sensitivity_moisture, 'sensitivity_moisture')
            makeSampleOutputFileLocal(model.sensitivity_suppression, 'sensitivity_suppression')
            makeSampleOutputFileLocal(model.sensitivity_ignitions, 'sensitivity_ignitions')

        if run_Bootstraps: 
            nboots = model.burnt_area_bootstraps.shape[0]
            for i in range(nboots):
                makeSampleOutputFileLocal(model.burnt_area_bootstraps[i], 
                                          'burnt_area_full_posterior', str(paramLoc)+ '-' + str(i))
        
        open(lock_file, 'a').close()

def listFiles_recursive(PATH, fname = '', pname = '', ext = '.nc'):
    return([os.path.join(dp, f) for dp, dn, filenames in os.walk(PATH) for f in filenames if fname in f and pname in dp and os.path.splitext(f)[1] == ext])

def npLogistic(x):
    return(1/(1+np.exp(-x)))

def build_distribution_from_boots(ensembles_dir, output_dir, var, period, 
                                  output_period = "monthly", transform = None, difference = False):
    ## setup directories and path
    mkDir(output_dir)
    output_dir_var = output_dir + var
    output_dir_var = output_dir_var + '-' + output_period + '/'
    
    mkDir(output_dir_var)
    outFile = output_dir_var + period + '.nc'
    print(outFile)
    if ignore_lock_summery or not os.path.exists(outFile):
        ## open and orgaise cube
        files = listFiles_recursive(ensembles_dir, period, var)
        if difference: 
            def fileFromID(files, id):
                return([file for file in files if id in file])
            
            files = fileFromID(files, "ensemble")
            cubesC = iris.load(fileFromID(files, experiments[0]))
            cubesE = iris.load(fileFromID(files, experiments[1]))
            cubes = iris.cube.CubeList([C - E for C, E in zip(cubesC, cubesE)])
        else:
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

    output_dir_sample = output_dir_exp + '/ensembles/'
    mkDir(output_dir_sample)
    
    ngap =int(params.shape[0]/n_posterior_sample)
    sample_nos = range(0, params.shape[0], ngap)
    for i in sample_nos: runSample(i, output_dir_sample, params, input_data, period)
    
    def summery_outputs(var):
        build_distribution_from_boots(output_dir_sample, output_dir_exp, 
                                      var, period, output_period = "annual")
        build_distribution_from_boots(output_dir_sample, output_dir_exp, 
                                      var, period)
    
    summery_outputs("burnt_area_mode")
    if run_Bootstraps: summery_outputs("burnt_area_full_posterior")
    if run_standard:
        summery_outputs("standard_fuel")
        summery_outputs("standard_moisture")
        summery_outputs("standard_ignitions")
        summery_outputs("standard_suppression")
    

for period in periods:
    for experiment in experiments:
        run_for_period(period, experiment)

    def summery_difference(var):    
        build_distribution_from_boots(output_dir, output_dir + '/difference/', 
                                      var, period, output_period = "annual", difference = True)
        build_distribution_from_boots(output_dir, output_dir + '/difference/', 
                                      var, period, difference = True)

    summery_difference("burnt_area_mode")

    if run_Bootstraps: summery_difference("burnt_area_full_posterior")
    if run_standard:
        summery_difference("standard_fuel")
        summery_difference("standard_moisture")
        summery_difference("standard_ignitions")
        summery_difference("standard_suppression")
        

