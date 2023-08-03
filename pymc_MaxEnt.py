import sys
sys.path.append('fire_model/')
sys.path.append('libs/')

from MaxEntFire import MaxEntFire
from BayesScatter import *

from read_variable_from_netcdf import *
from plot_maps import *
import os
from   io     import StringIO
import numpy  as np
import pandas as pd
import math
from scipy.special import logit, expit

import pymc  as pm
import pytensor
import pytensor.tensor as tt

import matplotlib.pyplot as plt
#import re
import arviz as az

from scipy.stats import wilcoxon

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
                             out_dir = 'outputs/', filename = '',
                             grab_old_trace = True,
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

    trace_callback = None
    try:
        if "SLURM_JOB_ID" in os.environ:
            def trace_callback(trace, draw):        
                if len(trace) % 10 == 0:
                    print('chain' + str(draw[0]))
                    print('trace' + str(len(trace)))
    except:
        pass        

    with pm.Model() as max_ent_model:
        ## set priors
        priors = {#"q":     pm.LogNormal('q', mu = 0.0, sigma = 1.0),
                  "betas": pm.Normal('betas', mu = 0, sigma = 1, shape = X.shape[1], 
                                      initval =np.repeat(0.5, X.shape[1])),
                   "powers": pm.Normal('powers', mu = 0, sigma = 1, shape = [2, X.shape[1]])#,
                   #"x2s": pm.Normal('x2s', mu = 0, sigma = 1, shape = [2, X.shape[1]])
                    # Maria: Add response curve priors
                 }
            # Will get used as follows:
            # y = beta[0] * X[:,0] + beta[1] * X[:,1] + beta[2] *X[:,2] + ....
            #           + powers[0,0] * X[:,0]^powers[0,1] + powers[1,0] * X[:,1]^powers[1,1] + powers[2,0] * X[:,2]^powers[2,1]
            #           + x2s[0,0] * (X2s[0,1] + X[:,0])^2 + x2s[1,0] * (X2s[1,1] + X[:,1])^2 + ...
            #           + Maria: response curve .....
        
        ## run model
        model = MaxEntFire(priors, inference = True)
        prediction = model.burnt_area(X)  
        
        ## define error measurement
        error = pm.DensityDist("error", prediction, logp = MaxEnt_on_prob, observed = Y)
                
        ## sample model
        attempts = 1
        while attempts <= 10:
            try:
                trace = pm.sample(niterations, return_inferencedata=True, 
                                  callback = trace_callback, *arg, **kw)
                attempts = 100
            except:
                print("sampling attempt " + str(attempts) + " failed. Trying a max of 10 times")
                attempts += 1
        ## save trace file
        trace.to_netcdf(trace_file)
    return trace


def train_MaxEnt_model(y_filen, x_filen_list, dir = '', filename_out = '',
                       dir_outputs = '',
                       frac_random_sample = 1.0,
                       subset_function = None, subset_function_args = None,
                       niterations = 100, cores = 4, model_title = '', grab_old_trace = False):
    
   
    Y, X, lmask, scalers = read_all_data_from_netcdf(y_filen, x_filen_list, 
                                                     add_1s_columne = True, dir = dir,
                                                     x_normalise01 = True, 
                                                     frac_random_sample = frac_random_sample,
                                                     subset_function = subset_function, 
                                                     subset_function_args = subset_function_args)
    
    dir_outputs = dir_outputs + '/' +  model_title
    if not os.path.exists(dir_outputs): os.makedirs(dir_outputs)

    trace = fit_MaxEnt_probs_to_data(Y, X, out_dir = dir_outputs, filename = filename, 
                                     niterations = niterations, cores = cores,
                                     grab_old_trace = grab_old_trace)
    
    return trace, scalers

def predict_MaxEnt_model(trace, y_filen, x_filen_list, scalers, dir = '', 
                         dir_outputs = '', model_title = '', filename_out = '',
                         subset_function = None, subset_function_args = None,
                         sample_for_plot = 1,
                         run_evaluation = False, run_projection = False):

    if not run_evaluation and not run_projection:
        return 

    Y, X, lmask, scalers = read_all_data_from_netcdf(y_filen, x_filen_list, 
                                                     add_1s_columne = True, dir = dir,
                                                     x_normalise01 = True, scalers = scalers,
                                                     subset_function = subset_function, 
                                                     subset_function_args = subset_function_args)
    if run_evaluation:
        Obs = read_variable_from_netcdf(y_filen, dir,
                                        subset_function = subset_function, 
                                        subset_function_args = subset_function_args)
    
    def select_post_param(name): 
        out = trace.posterior[name].values
        A = out.shape[0]
        B = out.shape[1]
        new_shape = ((A * B), *out.shape[2:])
        return np.reshape(out, new_shape)

    params = trace.to_dict()['posterior']
    params_names = params.keys()
    params = [select_post_param(var) for var in params_names]
    
    dir_outputs = dir_outputs + '/' + model_title + '/'
    if not os.path.exists(dir_outputs): os.makedirs(dir_outputs)
    dir_samples = dir_outputs + '/samples/' 
    if not os.path.exists(dir_samples): os.makedirs(dir_samples)
      
    def sample_model(i, run_name = 'control'):   
        file_sample = dir_samples + '/' + run_name + '-' + filename_out + '-' + str(i) + '.nc'
        
        if os.path.isfile(file_sample) and grab_old_trace:
            return pd.read_csv(file_sample).values.T[0]
           
        param_in = [param[i] if param.ndim == 1 else param[i,:] for param in params]
        param_in = dict(zip(params_names, param_in))
        out = MaxEntFire(param_in).burnt_area(X)
        pd.DataFrame(out).to_csv(file_sample, index = False)
        return out

    nits = np.prod(trace.posterior['betas'].values.shape[0:2])
    idx = range(0, nits, int(np.floor(nits/sample_for_plot)))
                         
    Obs = read_variable_from_netcdf(y_filen, dir,
                                    subset_function = subset_function, 
                                    subset_function_args = subset_function_args)

    Y, X, lmask, scalers = read_all_data_from_netcdf(y_filen, x_filen_list, 
                                                     add_1s_columne = True, dir = dir,
                                                     x_normalise01 = True, scalers = scalers,
                                                     subset_function = subset_function, 
                                                     subset_function_args = subset_function_args)
    x_copy = X[:, 1].copy()
    
    Sim = np.array(list(map(lambda id: sample_model(id, "control"), idx)))
    
    '''X[:, 1] = 1
    
    Sim2 = np.array(list(map(sample_model, idx)))
    
    fig, ax = plt.subplots()
    
    for rw in range(Sim.shape[0]):
        ax.plot(x_copy, (Sim[rw,:] / Sim2[rw,:]), '.')
    
    
    plt.show()
    
    
    set_trace()
    '''

    
    if run_evaluation:
        evaluate_model(filename_out, dir_outputs, Obs, Sim, lmask, levels, cmap)

    if run_projection:
        project_model(filename_out, dir_outputs, Sim, lmask, levels, cmap, eg_cube = Obs)


def plot_model_maps(Sim, lmask, levels, cmap, Obs = None, eg_cube = None, Nrows = 1, Ncols = 2):
    Sim = np.percentile(Sim, q = [10, 90], axis = 0)

    def plot_map(cube, plot_name, plot_n):
        plot_annual_mean(cube, levels, cmap, plot_name = plot_name, scale = 100*12, 
                     Nrows = Nrows, Ncols = Ncols, plot_n = plot_n)
  

    if eg_cube is None: eg_cube = Obs
    if Obs is not None: plot_map(Obs, "Observtations", 1)
    plot_map(insert_data_into_cube(Sim[0,:], eg_cube, lmask), "Simulation - 10%", Ncols - 1)
    plot_map(insert_data_into_cube(Sim[1,:], eg_cube, lmask), "Simulation - 90%", Ncols)

    

def evaluate_model(filename_out, dir_outputs, Obs, Sim, lmask, levels, cmap):
    qSim = np.percentile(Sim, q = np.arange(5, 100, 5), axis = 0)
    
    ax = plt.subplot(2, 3, 4)
    BayesScatter(Obs.data.flatten()[lmask], qSim.T, 0.000001, 0.000001, ax)
    plot_model_maps(Sim, lmask, levels, cmap, Obs, Nrows = 2, Ncols = 3)

    Y = Sim.reshape([Sim.shape[0], Obs.shape[0], int(Sim.shape[1]/Obs.shape[0])])
    X = Obs.data.flatten()[lmask].reshape([Obs.shape[0], Y.shape[2]])
    Yi = Y[:,:,0]
    Xi = X[:,0]

    pos = np.mean(X[np.newaxis, :, :] > Y, axis = 0)
    _, p_value = wilcoxon(pos - 0.5, axis = 0)
    apos = np.mean(pos, axis = 0)
    
    mask = lmask.reshape([ X.shape[0], int(lmask.shape[0]/X.shape[0])])[0]
    apos_cube = insert_data_into_cube(apos, Obs[0], mask)
    p_value_cube = insert_data_into_cube(p_value, Obs[0], mask)
    
    
    plot_annual_mean(apos_cube,[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
                     'RdYlBu_r',  plot_name = "mean bias",  Nrows = 2, Ncols = 3, plot_n = 5)

    plot_annual_mean(p_value_cube, np.array([0, 0.01, 0.05, 0.1, 0.5, 1.0]), 'copper',   
                     plot_name = "mean bias p-value",   Nrows = 2, Ncols = 3, plot_n = 6)
    
    plt.gcf().set_size_inches(8, 8)

    fig_dir = dir_outputs + '/figs/'
    if not os.path.exists(fig_dir): os.makedirs(fig_dir)

    plt.savefig(fig_dir + filename_out + '-evaluation.png')

    az.plot_trace(trace)
    plt.savefig(fig_dir + filename + '-traces.png')

def project_model(filename_out, dir_outputs, *arg, **kw):
    plot_model_maps(*arg, **kw)

    plt.gcf().set_size_inches(8*2/3, 6)
    fig_dir = dir_outputs + '/figs/'
    if not os.path.exists(fig_dir): os.makedirs(fig_dir)
    plt.savefig(fig_dir + filename_out + '-projections.png')

if __name__=="__main__":
    """ Running optimization and basic analysis. 
    Variables that need setting:
    For Optimization:
        dir_training -- The directory of the training data inclduing 
            dependant and independant variables
        y_filen -- filename of dependant variable (i.e burnt area)
        x_filen_list -- filanames of independant variables
            (ie bioclimate, landscape metrics etc)
        cores - how many chains to start (confusiong name, I know).
            When running on slurm, also number of cores
        fraction_data_for_sample -- fraction of gridcells used for optimization
        niterations -- number of iterations or samples )after warm-up) in optimixation for each
            chain. Equilivent to number of ensemble members.
        months_of_year --- which months to extact on training and projecting
        grab_old_trace -- Boolean. If True and there's an appripritate looking old trace file, 
            then  optimisation is skipped that file is loaded instead. 
            This isn't totally infalable, so if doing a final run and in doubt, set to False
    For Projection/evaluating:
        dir_outputs -- where stuff gets outputted
        dir_projecting -- The directory of the data used for prjections. 
            Should contain same files for independant varibales as dir_training 
            (though you should be able to adpated this easily if required). 
            Doesnt need dependant variable, but if there, this will (once
            we've implmented it) attempt some evaluation.
        sample_for_plot -- how many iterations (samples) from optimixation should be used 
            for plotting and evaluation.
        levels -- levels on the colourbar on observtation and prodiction maps
        cmap -- levels on teh colourbar on observtation and prodiction maps
    Returns:
        trace file, maps, etc (to be added too)
    """

    """ 
        SETPUT 
    """
    """ optimization """

    model_title = 'Example_model-X2'

    dir_training = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/"
    #dir_training = "/gws/nopw/j04/jules/mbarbosa/driving_and_obs_overlap/AllConFire_2000_2009/"

    y_filen = "GFED4.1s_Burned_Fraction.nc"
    
    x_filen_list=["trees.nc", "pr_mean.nc", "consec_dry_mean.nc", 
                  "lightn.nc", "popDens.nc",
                  "crop.nc", "pas.nc", 
                  "humid.nc", "csoil.nc", "tas_max.nc",
                  "totalVeg.nc"]
    
    x_filen_list=["trees.nc", "consec_dry_mean.nc", 
                  "lightn.nc", "popDens.nc",
                  "crop.nc", "pas.nc", 
                  "tas_max.nc",
                  "totalVeg.nc"]


    grab_old_trace = True
    cores = 2
    fraction_data_for_sample = 0.05
    niterations = 100

    months_of_year = [7]
    
    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = dir_training

    sample_for_plot = 20

    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    cmap = 'OrRd'

    run_evaluation = True
    run_projection = True
     
    """ 
        RUN optimization 
    """
    subset_function = sub_year_months
    subset_function_args = {'months_of_year': months_of_year}

    filename = '_'.join([file[:-3] for file in x_filen_list]) + \
              '-frac_points_' + str(fraction_data_for_sample) + \
              '-Month_' +  '_'.join([str(mn) for mn in months_of_year])

    #### Optimize
    trace, scalers = train_MaxEnt_model(y_filen, x_filen_list, dir_training, 
                                        filename, dir_outputs,
                                        fraction_data_for_sample,
                                        subset_function, subset_function_args,
                                        niterations, cores, model_title, grab_old_trace)


    """ 
        RUN projection 
    """
    predict_MaxEnt_model(trace, y_filen, x_filen_list, scalers, dir_projecting,
                         dir_outputs, model_title, filename,
                         subset_function, subset_function_args,
                         sample_for_plot, 
                         run_evaluation = run_evaluation, run_projection = run_projection)
    
    
