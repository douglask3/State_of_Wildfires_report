#import multiprocessing as mp 
#mp.set_start_method('forkserver')
import sys
sys.path.append('fire_model/')
sys.path.append('libs/')

from MaxEntFire import MaxEntFire

from read_variable_from_netcdf import *
from plot_maps import *
import os
from   io     import StringIO
import numpy  as np
import math

import pymc  as pm
import pytensor
import pytensor.tensor as tt
#from   aesara import tensor as tt

import matplotlib.pyplot as plt
import re

import sys
import arviz as az


def MaxEnt_on_prob(BA, fx):
    """calculates the log-transformed continuous logit likelihood for x given mu when x 
       and mu are probabilities between 0-1. 
       Works with tensor variables.   
    Arguments:
        x -- x in P(x|mu). tensor 1-d array
	mu -- mu in P(x|mu). tensor 1-d array
    Returns:
        1-d tensor array of liklihoods.
    """
    fx = tt.switch(
        tt.lt(fx, 0.0000000000000000001),
        0.0000000000000000001, fx)
    return BA*tt.log(fx) + (1.0-BA)*tt.log((1-fx))   

def fit_MaxEnt_probs_to_data(Y, X, niterations, 
                             out_dir = 'outputs/', filename = '', grab_old_trace = True,
                             *arg, **kw):
    """ Bayesian inerence routine that fits independant variables, X, to dependant, Y.
        Based on the MaxEnt solution of probabilities. 
    Arguments:
        Y-- dependant variable as numpy 1d array
	X -- numpy 2d array of indepenant variables, each columne a different variable
	niterations -- number of iterations per chain when sampling the postior during 
                NUTS inference 
		(note default chains is normally 2 and is set by *args or **kw)
	out_dir --string of path to output location. This is where the traces netcdf file 
                will be saved.
		Defauls is 'outputs'.
	filename -- string of the start of the traces output name. Detault is blank. 
		Some metadata will be saved in the filename, so even blank will 
                save a file.
	grab_old_trace -- Boolean. If True, and a filename starting with 'filename' and 
                containing some of the same setting (saved in filename) exists,  it will open 
                and return this rather than run a new one. Not all settings are saved for 
                identifiation, so if in doubt, set to 'False'.
	*args, **kw -- arguemts passed to 'pymc.sample'

    Returns:
        pymc traces, returned and saved to [out_dir]/[filneame]-[metadata].nc
    """

    trace_file = out_dir + '/' + filename + '-nvariables_' + '-ncells_' + str(X.shape[0]) + \
                str(X.shape[1]) + '-niterations_' + str(niterations * cores) + '.nc'
    
    ## check if trace file exsts and return if wanted
    if os.path.isfile(trace_file) and grab_old_trace: 
        return az.from_netcdf(trace_file)

    with pm.Model() as max_ent_model:
        ## set priorts
        betas = pm.Normal('betas', mu = 0, sigma = 1, shape = X.shape[1], 
                          initval =np.repeat(0.5, X.shape[1]))

        powers = pm.Normal('powers', mu = 0, sigma = 1, shape = [2, X.shape[1]])
        ## build model
        
        prediction = MaxEntFire(betas, powers, inference = True).fire_model(X)  
        
        
        ## define error measurement
        error = pm.DensityDist("error", prediction, logp = MaxEnt_on_prob, observed = Y)
                
        ## sample model
        trace = pm.sample(niterations, return_inferencedata=True, *arg, **kw)
        ## save trace file
        trace.to_netcdf(trace_file)
    return trace


if __name__=="__main__":
    dir = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/"
    #dir = "/gws/nopw/j04/jules/mbarbosa/driving_and_obs_overlap/AllConFire_2000_2009/"
    
    dir_outputs = 'outputs/'
    grab_old_trace = True
    y_filen = "GFED4.1s_Burned_Fraction.nc"
    

    x_filen_list=["precip.nc", "lightn.nc", "crop.nc", "humid.nc","vpd.nc", "csoil.nc", 
                  "lightn.nc", "rhumid.nc", "cveg.nc", "pas.nc", "soilM.nc", 
                   "totalVeg.nc", "popDens.nc", "trees.nc"]


    fraction_data_for_sample = 0.1    

    niterations = 100
    sample_for_plot = 20

    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    cmap = 'OrRd'

    months_of_year = [7]

    cores = 4

    def read_data_to_cols(frac_random_sample, *args, **kw):
        return read_all_data_from_netcdf(y_filen, x_filen_list, 
                                         add_1s_columne = True, dir = dir,
                                         x_normalise01 = True, 
                                         frac_random_sample = frac_random_sample,
                                         subset_function = sub_year_months, 
                                         subset_function_args = {'months_of_year': 
                                                                 months_of_year})

    Y, X, lmask, scalers = read_data_to_cols(fraction_data_for_sample)
    
    filename = '_'.join([file[:-3] for file in x_filen_list]) + \
              '-frac_points_' + str(fraction_data_for_sample) + \
              '-Month_' +  '_'.join([str(mn) for mn in months_of_year])
    
    trace = fit_MaxEnt_probs_to_data(Y, X, out_dir = dir_outputs, filename = filename, 
                                     niterations = niterations, cores = cores,
                                     grab_old_trace = grab_old_trace)
    
    az.plot_trace(trace)
    plt.savefig('figs/' + filename + '-traces.png')

    Y, X, lmask, scalers = read_data_to_cols(1.0)
    Obs = read_variable_from_netcdf(y_filen, dir, subset_function = sub_year_months, 
                                     subset_function_args = {'months_of_year': months_of_year})

    def select_post_param(name): 
        out = trace.posterior[name].values
        A = out.shape[0]
        B = out.shape[1]
        new_shape = ((A * B), *out.shape[2:])
        return np.reshape(out, new_shape)

    def sample_model(i): 
        powers =select_post_param('powers')[i,:]
        betas =select_post_param('betas')[i,:]
        return MaxEntFire(betas, powers).fire_model(X)

    nits = np.prod(trace.posterior['betas'].values.shape[0:2])
    idx = range(0, nits, int(np.floor(nits/sample_for_plot)))

    Sim = np.array(list(map(sample_model, idx)))
    Sim = np.percentile(Sim, q = [10, 90], axis = 0)
    
    def insert_sim_into_cube(x):
        Pred = Obs.copy()
        pred = Pred.data.copy().flatten()
        
        pred[lmask] = x
        Pred.data = pred.reshape(Pred.data.shape)
        return(Pred)

    def plot_map(cube, plot_name, plot_n):
        plot_annual_mean(cube, levels, cmap, plot_name = plot_name, scale = 100*12, 
                     Nrows = 1, Ncols = 3, plot_n = plot_n)
  
    plot_map(Obs, "Observtations", 1)
    plot_map(insert_sim_into_cube(Sim[0,:]), "Simulation - 10%", 2)
    plot_map(insert_sim_into_cube(Sim[1,:]), "Simulation - 90%", 3)
    plt.gcf().set_size_inches(8, 6)
    plt.savefig('figs/' + filename + '-maps.png')
    
    '''
    #Run the model with first iteration
    simulation1 = fire_model(trace.posterior['betas'].values[0,0,:], X, False)

    #Plot against observations (Y)
    plt.plot(Y, simulation1, '.')
    plt.show()

    #when developing plots, use betas = trace.posterior['betas'].values[0,0,:]
    and run with model with fire_model(betas, X, False)
    '''
    set_trace()
