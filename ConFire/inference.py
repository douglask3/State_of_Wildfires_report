import os
from   io     import StringIO
import numpy  as np
import pandas as pd
import csv
import math

import pymc  as pm
#from   pymc3.backends import SQLite
from   scipy  import optimize
from   aesara import tensor as tt

import matplotlib.pyplot as plt
import re

import sys
sys.path.append('../')
from ConFire import ConFire

from pdb import set_trace as browser


datDir       =  "/data/users/dkelley/ConFIRE_ISIMIP/isimip3_inputs/Global/inference_data/"
param_outpath = "../ConFIRE_ISIMIP/outputs/isimip3/params-for_sampling/"
param_file = "with_vpd"
sample_pc = 1
nChains = 2
nIterations = 1000
nTune = 500


#import corner


def npLogit(x):
    return np.log(x/(1-x))
    
def ttLogit(x):
    return tt.log(x/(1-x))

#########################################################
## loading data                                        ##
#########################################################

def load_with_buffer(filename, line_select, **kwargs):
    s_buf = StringIO()
    line_select = np.sort(line_select)
    with open(filename) as file:
        count = -1
        lineN = -1
        for line in file:
            lineN += 1
            if lineN == 0 or lineN == line_select[count]:
                s_buf.write(line)
                
                count += 1
                if count == len(line_select): break
            
    s_buf.seek(0)
    df = pd.read_csv(s_buf,**kwargs)
    return df

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f): pass
    return i + 1

def openDat(datPath):
    datPath = datDir + '/' + datPath
    DATAPATH = os.path.expanduser(datPath)
    
    nlines      = file_len(DATAPATH)
    npoints     = round(sample_pc * nlines / 100)
    line_select = np.random.choice(range(0, nlines), npoints, False)
    line_select = line_select[line_select > 0]
    fd          = load_with_buffer(DATAPATH, line_select)

    fd['vpd'].values[fd['vpd'] < 0.0] = 0.0

    BA = fd["fireObs"].values
    BA[BA < 10e-9] = 10e-9
    
    BA = npLogit(BA)
    fd["fireObs"].values[:] = BA[:] 
    
    vpd = fd['vpd'].values  
    vpd[vpd < 0.0] = 0.0 
    return fd

fds = [openDat(f) for f in os.listdir(datDir)]


#########################################################
## inference functions                                 ##
#########################################################

def make_zero_inflated_normal(x, mu, sigma, pz):
    '''return tt.sw1itch(
        tt.lt(x, -150),
        -p0,
        -(1.0 - p0) *(1.0/(sigma * 2.506))*tt.exp(-0.5 * ((x-mu)/sigma)**2)
    )
    '''
    return tt.switch(
        tt.lt(x, -30),
        tt.log(pz),
        tt.log(1-pz) - tt.log(sigma * tt.sqrt(2*math.pi)) - ((x-mu)**2)/(2*sigma**2)
    )

def runInference(fd, outfile):

    with pm.Model() as fire_error:
        

        params = {"fuel_x0": pm.Normal     ('fuel_x0'     , 0.5, 0.25),
                  "fuel_k": pm.Exponential('fuel_k'      , 1.0      ),
                  "c_cveg":  pm.Lognormal ('c_cveg'      , 0.0, 1.0 ),
                  "c_csoil":  pm.Lognormal ('c_csoil'      , 0.0, 1.0 ),
                  "c_vpd":  pm.Lognormal ('c_vpd'      , 0.0, 1.0 ),
                  "moisture_x0": pm.Normal     ('moisture_x0' , 0.5, 0.25),
                  "moisture_k": pm.Exponential('moisture_k'  , 1.0      ),
                  "wd_pg": pm.Exponential('wd_pg'       , 1.0      ),
                  "k_vpd1": pm.LogitNormal('k_vpd1' , 0.0, 1.0),
                  "k_vpd2": pm.LogitNormal('k_vpd2' , 0.0, 1.0),
                  "kM": pm.LogitNormal('kM' , 0.0, 1.0),
                  "pT": pm.Lognormal  ('pT' , 0.0, 1.0),
                  "c_emc": pm.Lognormal  ('c_emc'       , 0.0, 1.0 ),
                  "c_trees": pm.Lognormal  ('c_trees'     , 0.0, 1.0 ),
                  "ignition_x0": pm.Normal     ('ignition_x0', 1000.0, 50.0),
                  "ignition_k": pm.Exponential('ignition_k' , 100.0     ),
                  "bck_ignitions": pm.Lognormal  ('bck_ignitions', 0.0, 10.0 ),
                  "suppression_x0": pm.Normal ('suppression_x0'  , 0.5, 0.25),
                  "suppression_k": pm.Exponential('suppression_k', 1.0     ),
                  "max_f": pm.LogitNormal('max_f'           , 0.0, 1.0)}        
        p0 = pm.Uniform('p0', 0.0, 1.0)
        p1 = pm.Lognormal('p1', 0.0, 1.0)
        sigma = pm.Lognormal('sigma', 0.0, 1.0)
        prediction = ConFire(fd, params, True).burnt_area_mode
        
               
        pz = 1.0 - (prediction**p1) * (1.0 - p0)
        prediction = ttLogit(prediction)
        
        error = pm.DensityDist("error", prediction, sigma, pz, 
                               logp = make_zero_inflated_normal, 
                               observed = fd["fireObs"].values)
        
        # set the step-method (criteria algorithm for moving around information space)        
        istep = pm.Metropolis()
        
        # do the sampling
        idata = pm.sample(step=istep, chains = nChains) #, start=start, trace=db_save      
        posterior = idata.posterior.to_dataframe()
        browser()
        posterior.to_csv(param_outpath + '/' + param_file + '-' + 
                         outfile + str(nChains) + '.csv', index=False)

for fd, outfile in zip(fds,os.listdir(datDir)):
    runInference(fd, outfile) 
