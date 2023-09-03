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

def combine_path_and_make_dir(path1, path2):
    path = path1 + '/'+ path2 + '/'
    if not os.path.exists(path): os.makedirs(path)
    return path

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


def fit_MaxEnt_probs_to_data(Y, X, niterations, *arg, **kw):
    """ Bayesian inerence routine that fits independant variables, X, to dependant, Y.
        Based on the MaxEnt solution of probabilities. 
    Arguments:
        Y-- dependant variable as numpy 1d array
	X -- numpy 2d array of indepenant variables, each columne a different variable
	niterations -- number of iterations per chain when sampling the postior during 
                NUTS inference 
		(note default chains is normally 2 and is set by *args or **kw)
		Defauls is 'outputs'.
	*args, **kw -- arguemts passed to 'pymc.sample'

    Returns:
        pymc traces, returned and saved to [out_dir]/[filneame]-[metadata].nc
    """

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
        nvars = X.shape[1]
        priors = {#"q":     pm.LogNormal('q', mu = 0.0, sigma = 1.0),
                  "lin_betas": pm.Normal('lin_betas', mu = 0, sigma = 1, shape = nvars),
                  "pow_betas": pm.Normal('pow_betas', mu = 0, sigma = 1, shape = nvars),
                  "pow_power": pm.Normal('pow_power', mu = 0, sigma = 1, shape = nvars),
                  "x2s_betas": pm.Normal('x2s_betas', mu = 0, sigma = 1, shape = nvars),
                  "x2s_X0"   : pm.Normal('x2s_X0'   , mu = 0, sigma = 1, shape = nvars),
                  "comb_betas": pm.Normal('comb_betas', mu = 0, sigma = 1, shape = nvars),
                  "comb_X0": pm.Normal('comb_X0', mu = 0.5, sigma = 1, shape = nvars),
                  "comb_p": pm.Normal('comb_p', mu = 0, sigma = 1 , shape = nvars)
                  
                   #"x2s": pm.Normal('x2s', mu = 0, sigma = 1, shape = [2, X.shape[1]])
                    # Maria: Add response curve priors
                    #"X0": pm.Normal('X0', mu = 0.5, sigma = 1, shape = X.shape[1])
                    #"p": pm.Normal('p', mu = 0, sigma = 1 , shape = X.shape[1])
                    #"gama": pm.Normal('gama', mu = 0 , sigma = 1, shape = X.shape[1])
                 }
            # Will get used as follows:
            # y = beta[0] * X[:,0] + beta[1] * X[:,1] + beta[2] *X[:,2] + ....
            #           + powers[0,0] * X[:,0]^powers[0,1] + powers[1,0] * X[:,1]^powers[1,1] + powers[2,0] * X[:,2]^powers[2,1]
            #           + x2s[0,0] * (X2s[0,1] + X[:,0])^2 + x2s[1,0] * (X2s[1,1] + X[:,1])^2 + ...
            #           + Maria: response curve .....
        
        ## run model
        model = MaxEntFire(priors, inference = True)
        prediction = model.burnt_area_uninflated(X)  
        
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
    return trace


def train_MaxEnt_model(y_filen, x_filen_list, dir = '', filename_out = '',
                       dir_outputs = '',
                       frac_random_sample = 1.0,
                       subset_function = None, subset_function_args = None,
                       niterations = 100, cores = 4, model_title = '', grab_old_trace = False):
    '''
	out_dir --string of path to output location. This is where the traces netcdf file 
                will be saved.
        
	filename_out -- string of the start of the traces output name. Detault is blank. 
		Some metadata will be saved in the filename, so even blank will 
                save a file.
	grab_old_trace -- Boolean. If True, and a filename starting with 'filename' and 
                containing some of the same setting (saved in filename) exists,  it will open 
                and return this rather than run a new one. Not all settings are saved for 
                identifiation, so if in doubt, set to 'False'.
    '''

    dir_outputs = combine_path_and_make_dir(dir_outputs, model_title)
    out_file =   filename_out + '-nvariables_' + \
                 '-frac_random_sample' + str(frac_random_sample) + \
                 '-nvars_' +  str(len(x_filen_list)) + \
                 '-niterations_' + str(niterations * cores)
    
    trace_file = dir_outputs + '/trace-'   + out_file + '.nc'
    scale_file = dir_outputs + '/scalers-' + out_file + '.csv'
    
    ## check if trace file exsts and return if wanted
    if os.path.isfile(trace_file) and os.path.isfile(scale_file) and grab_old_trace:
        return az.from_netcdf(trace_file), pd.read_csv(scale_file).values   
    print("opening data for inference")
    Y, X, lmask, scalers = read_all_data_from_netcdf(y_filen, x_filen_list, 
                                                     add_1s_columne = True, dir = dir,
                                                     x_normalise01 = True, 
                                                     frac_random_sample = frac_random_sample,
                                                     subset_function = subset_function, 
                                                     subset_function_args = subset_function_args)
    print("Running trace")
    trace = fit_MaxEnt_probs_to_data(Y, X,niterations = niterations, cores = cores)
    
    ## save trace file
    trace.to_netcdf(trace_file)
    pd.DataFrame(scalers).to_csv(scale_file, index = False)
    return trace, scalers

def predict_MaxEnt_model(trace, y_filen, x_filen_list, scalers, dir = '', 
                         dir_outputs = '', model_title = '', filename_out = '',
                         subset_function = None, subset_function_args = None,
                         sample_for_plot = 1,
                         run_evaluation = False, run_projection = False, grab_old_trace = False,
                         *args, **kw):

    if not run_evaluation and not run_projection:
        return 

    Y, X, lmask, scalers = read_all_data_from_netcdf(y_filen, x_filen_list, 
                                                     add_1s_columne = True, dir = dir,
                                                     x_normalise01 = True, scalers = scalers,
                                                     subset_function = subset_function, 
                                                     subset_function_args = subset_function_args)
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
    
    dir_outputs = combine_path_and_make_dir(dir_outputs, model_title)
    dir_samples = combine_path_and_make_dir(dir_outputs, '/samples/')     
    dir_samples = combine_path_and_make_dir(dir_samples, filename_out)
      
    def sample_model(i, run_name = 'control'):   
        dir_sample =  combine_path_and_make_dir(dir_samples, run_name)
        file_sample = dir_sample + '/sample' + str(i) + '.nc'
        
        if os.path.isfile(file_sample) and grab_old_trace:
            return iris.load_cube(file_sample)
        print("Generating Sample:" + file_sample)
        param_in = [param[i] if param.ndim == 1 else param[i,:] for param in params]
        param_in = dict(zip(params_names, param_in))
        out = MaxEntFire(param_in).burnt_area(X)
        out = insert_data_into_cube(out, Obs, lmask)
        iris.save(out, file_sample)
        
        return out
    
    nits = len(trace.posterior.chain)*len(trace.posterior.draw)
    idx = range(0, nits, int(np.floor(nits/sample_for_plot)))
         
    Sim = np.array(list(map(lambda id: sample_model(id, "control"), idx)))
    
    
    x_copy = X[:, 1].copy()
    for col in range(X.shape[1]-1):
        x_copy = X[:, col].copy()  # Copy the values of the current column
        
        variable = x_filen_list[col].replace('.nc', '')
        print(col)
        
        X[:, col] = 0.0  # Set the current column to 0
        
        Sim2 = np.array(list(map(lambda id: sample_model(id, variable + '_to_zero-'), 
                                 idx))) # Sample model for the modified column
        
        fcol = math.floor(math.sqrt(X.shape[1]))
        frw = math.ceil(X.shape[1]/fcol)
        
        ax = plt.subplot(frw,fcol, col + 1)  # Select the corresponding subplot
        
        def non_masked_data(cube):
            return cube.data[cube.data.mask == False].data

        for rw in range(len(Sim)):
            ax.plot(x_copy, non_masked_data(Sim[rw]) / non_masked_data(Sim2[rw]), '.')  # Plot the data for the current variable
        
        X[:, col] = x_copy 
        
    plt.show()
    
    if run_evaluation:
        evaluate_model(filename_out, dir_outputs, Obs, Sim, lmask, *args, **kw)

    if run_projection:
        project_model(filename_out, dir_outputs, Sim, lmask, eg_cube = Obs, *args, **kw)


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
    set_trace()
    qSim = np.percentile(Sim, q = np.arange(5, 100, 5), axis = 0)
    
    ax = plt.subplot(2, 3, 4)
    BayesScatter(Obs.data.flatten()[lmask], qSim.T, 0.000001, 0.000001, ax)
    plot_model_maps(Sim, lmask, levels, cmap, Obs, Nrows = 2, Ncols = 3)
    set_trace()
    Y = Sim.reshape([Sim.shape[0], Obs.shape[0], int(Sim.shape[1]/Obs.shape[0])])
    X = Obs.data.flatten()[lmask].reshape([Obs.shape[0], Y.shape[2]])
    Yi = Y[:,:,0]
    Xi = X[:,0]
    #set_trace()
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
    
    fig_dir = combine_path_and_make_dir(dir_outputs, '/figs/')

    plt.savefig(fig_dir + filename_out + '-evaluation.png')

    az.plot_trace(trace)
    plt.savefig(fig_dir + filename + '-traces.png')

def project_model(filename_out, dir_outputs, *arg, **kw):
    plot_model_maps(*arg, **kw)

    plt.gcf().set_size_inches(8*2/3, 6)
    fig_dir = combine_path_and_make_dir(dir_outputs, '/figs/')
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

    model_title = 'Example_model'

    dir_training = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/"
    #dir_training = "/gws/nopw/j04/jules/mbarbosa/driving_and_obs_overlap/AllConFire_2000_2009/"
    #dir_training = "D:/Doutorado/Sanduiche/research/maxent-variables/2002-2011/"

    y_filen = "GFED4.1s_Burned_Fraction.nc"
    #y_filen = "Area_burned_NAT.nc"
    #y_filen = "Area_burned_NON2.nc"

    #x_filen_list=["trees.nc", "pr_mean.nc", "consec_dry_mean.nc", 
                  #"lightn.nc", "popDens.nc",
                  #"crop.nc", "pas.nc", 
                  #"humid.nc", "csoil.nc", "tas_max.nc",
                  #"totalVeg.nc"]
    

    x_filen_list=["trees.nc", "pr_mean.nc", "consec_dry_mean.nc", 
                  "lightn.nc", "popDens.nc",
                  "crop.nc", "pas.nc", 
                  "humid.nc", "csoil.nc", "tas_max.nc",
                  "totalVeg.nc"]
    
    #x_filen_list=["trees.nc", "consec_dry_mean.nc", 
    #              "lightn.nc", "popDens.nc",
    #              "crop.nc", "pas.nc", 
    #              "tas_max.nc",
    #              "totalVeg.nc", "MPA.nc"]

    #x_filen_list=["consec_dry_mean.nc", "Savanna.nc", "cveg.nc", "rhumid.nc",
    #              "lightn.nc", "popDens.nc", "Forest.nc", "precip.nc",
    #              "crop.nc", "pas.nc", "Grassland.nc",
    #              "tas_max.nc", "tas_mean.nc",
    #              "totalVeg.nc", "vpd.nc", "csoil.nc"]



    grab_old_trace = True # set to True till you get the code running. Then set to False when you start adding in new response curves

    cores = 2
    fraction_data_for_sample = 0.05
    niterations = 100

    months_of_year = [7]
    
    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = dir_training

    sample_for_plot = 50

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
                         run_evaluation = run_evaluation, run_projection = run_projection,
                         grab_old_trace = grab_old_trace,
                         levels = levels, cmap = cmap)
    
    
