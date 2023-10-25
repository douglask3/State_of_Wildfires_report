import sys
sys.path.append('fire_model/')
sys.path.append('libs/')

from MaxEntFire import MaxEntFire
from BayesScatter import *
#from train import *

from read_variable_from_netcdf import *
from combine_path_and_make_dir import * 
from plot_maps import *
import os
from   io     import StringIO
import numpy  as np
import pandas as pd
import math
from scipy.special import logit, expit

import matplotlib.pyplot as plt
import arviz as az

from scipy.stats import wilcoxon
from sklearn.metrics import mean_squared_error


def predict_MaxEnt_model(trace, y_filen, x_filen_list, scalers, CA_filen = None, dir = '', 
                         dir_outputs = '', model_title = '', filename_out = '',
                         subset_function = None, subset_function_args = None,
                         sample_for_plot = 1,
                         run_evaluation = False, run_projection = False, grab_old_trace = False,
                         *args, **kw):

    """ Runs prediction and evalutation of the sampled model based on previously run trace.
    Arguments:
        trace - pymc traces nc file, probably from a 'train_MaxEnt_model' run
	y_filen -- filename of dependant variable (i.e burnt area)
        x_filen_list -- filanames of independant variables
	scalers -- the scalers used during the generation of the trace file that scales 
            x data between 0 and 1 (not that it might not scale the opened variables 
            here between 0 and 1 as different data maybe selected for evaluation)
        dir -- dir y_filen and  x_filen_list are stored in. default it current dir
        dir_outputs -- directiory where all outputs are stored. String
        model_title - title of model run. A str default to 'no_name'. Used to initially to name 
                the dir everythings stored in.
        filename_out -- string of the start of the traces output name. Detault is blank. 
		Some metadata will be saved in the filename, so even blank will 
                save a file.
        subset_function -- a list of constrain function useful for constraining and resticting 
                data to spatial locations and time periods/months. Default is not to 
                constrain (i.e "None" for no functions")
        subset_function_args -- list of arguements that feed into subset_function
        sample_for_plot -- fraction of gridcells used for optimization
        run_evaluation -- Logical, default False - run the standard evalaution 
                plots/metrics or not
        run_projection -- Logical, default False - run projections
        grab_old_trace -- Boolean. If True, and a filename starting with 'filename' and 
                containing some of the same setting (saved in filename) exists,  it will open 
                and return this rather than run new samples. Not all settings are saved for 
                identifiation, so if in doubt, set to 'False'.
	*args, **kw -- arguemts passed to 'evaluate_model' and 'project_model'
    
    Returns:
        look in dir_outputs + model_title, and you'll see figure and tables from evaluation, 
        projection, reponse curevs, jacknifes etc (not all implmenented yet)
    """
    if not run_evaluation and not run_projection:
        return 

    common_args = {
    'y_filename': y_filen,
    'x_filename_list': x_filen_list,
    'add_1s_columne': True,
    'dir': dir,
    'x_normalise01': True,
    'subset_function': subset_function,
    'subset_function_args': subset_function_args
}

    if CA_filen is not None:
        Y, X, lmask, scalers, CA = read_all_data_from_netcdf(CA_filename = CA_filen, **common_args)
        
    else:
        
        Y, X, lmask, scalers = read_all_data_from_netcdf(**common_args)
    
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
        coord = iris.coords.DimCoord(i, "realization")
        out.add_aux_coord(coord)
        iris.save(out, file_sample)
        
        return out
    
    nits = len(trace.posterior.chain)*len(trace.posterior.draw)
    idx = range(0, nits, int(np.floor(nits/sample_for_plot)))
    
    def runSim(id_name):  
        out = np.array(list(map(lambda id: sample_model(id, id_name), idx)))
        return iris.cube.CubeList(out).merge_cube()
               
    Sim = runSim("control") 
    '''  
    for col_to_keep in range(X.shape[1]-1):
        other_cols = np.arange(X.shape[1]-1)  # Create an array of all columns
        other_cols = other_cols[other_cols != col_to_keep]  # Exclude col_to_keep
        original_X = X[:, other_cols].copy()
        X[:, other_cols] = 0.0  
        
        Sim2 = runSim("_to_zero")
        
        fcol = math.floor(math.sqrt(X.shape[1]))
        frw = math.ceil(X.shape[1]/fcol)
        
        ax = plt.subplot(frw,fcol, col_to_keep + 1)  # Select the corresponding subplot
        
        variable_name = x_filen_list[col_to_keep].replace('.nc', '')
        ax.set_title(variable_name)
        
        def non_masked_data(cube):
            return cube.data[cube.data.mask == False].data
        
        
        num_bins = 10
        hist, bin_edges = np.histogram(X[:, col_to_keep], bins=num_bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        median_values = []
        percentile_10 = []
        percentile_90 = []
        
        for i in range(num_bins):
        
            mask = (X[:, col_to_keep] >= bin_edges[i]) & (X[:, col_to_keep] < bin_edges[i + 1])
            values_in_bin = []
        
            for rw in range(Sim.shape[0]):
                sim_final = non_masked_data(Sim[rw]) - non_masked_data(Sim2[rw])                      
                values_in_bin.append(sim_final[mask])
            values_in_bin = np.array(values_in_bin).flatten()    
            #set_trace()    
            median_values.append(np.median(values_in_bin))
            percentile_10.append(np.percentile(values_in_bin, 10))
            percentile_90.append(np.percentile(values_in_bin, 90))
      
                
        set_trace()   
        ax.plot(bin_centers, median_values, marker='.', label='Median')
        ax.fill_between(bin_centers, percentile_10, percentile_90, alpha=0.3, label='10th-90th Percentiles')                           
                
        X[:, other_cols] = original_X
    
    #plt.show()
    #set_trace()     
      
    fig_dir = combine_path_and_make_dir(dir_outputs, '/figs/')
    
    plt.savefig(fig_dir + 'control-response-curves.png')    
    '''
    
    #plot sensitivity response curves
    plt.figure(figsize=(14, 12))
    for col in range(X.shape[1]-1):
        x_copy = X[:, col].copy()  # Copy the values of the current column
        
        print(col)

        dx = 0.001
        X[:, col] -= dx/2.0  # Subtract 0.1 of all values for the current column

        Sim3 = runSim("subtract_01")    
        
        X[:, col] = x_copy #restore values
        
        X[:, col] += dx/2.0 #add 0.1 to all values for the current column
        
        Sim4 = runSim("add_01") 
        
        fcol = math.floor(math.sqrt(X.shape[1]))
        frw = math.ceil(X.shape[1]/fcol)
        
        ax = plt.subplot(frw,fcol, col + 1)
        variable_name = x_filen_list[col].replace('.nc', '')
        ax.set_title(variable_name)
        
        def non_masked_data(cube):
            return cube.data[cube.data.mask == False].data
            
        num_bins = 20
        hist, bin_edges = np.histogram(X[:, col], bins=num_bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        median_values = []
        percentile_10 = []
        percentile_90 = []
        
        for i in range(num_bins):
        
            mask = (X[:, col] >= bin_edges[i]) & (X[:, col] < bin_edges[i + 1])
            values_in_bin = []
        
            for rw in range(Sim.shape[0]):
                sim_final = (non_masked_data(Sim3[rw]) - non_masked_data(Sim4[rw]))/dx                      
                values_in_bin.append(sim_final[mask])
            values_in_bin = np.array(values_in_bin).flatten()   
            median_values.append(np.median(values_in_bin))
            try:
                percentile_10.append(np.percentile(values_in_bin, 10))
                percentile_90.append(np.percentile(values_in_bin, 90))
            except:
                percentile_10.append(np.nan)
                percentile_90.append(np.nan)
            
            
        #set_trace()    
        ax.plot(bin_centers, median_values, marker='.', label='Median')
        ax.fill_between(bin_centers, percentile_10, percentile_90, alpha=0.3, label='10th-90th Percentiles')      
        
        #for rw in range(Sim.shape[0]):
            #ax.plot(X[:, col], non_masked_data(Sim3[rw]) - non_masked_data(Sim4[rw]), '.', color = "darkred", markersize = 0.5, linewidth=0.5)
    
        X[:, col] = x_copy 
    
    plt.show()
    set_trace() 
    
    fig_dir = combine_path_and_make_dir(dir_outputs, '/figs/')
    
    plt.savefig(fig_dir + 'sensitivity-response-curves.png')  

    if run_evaluation:
        evaluate_model(filename_out, dir_outputs, Obs, Sim, lmask, *args, **kw)

    if run_projection:
        project_model(filename_out, dir_outputs, Sim, lmask, eg_cube = Obs, *args, **kw)


def plot_model_maps(Sim, lmask, levels, cmap, Obs = None, eg_cube = None, Nrows = 1, Ncols = 2):
    Sim = Sim.collapsed('realization', iris.analysis.PERCENTILE, 
                          percent=[10, 90])

    def plot_map(cube, plot_name, plot_n):
        plot_annual_mean(cube, levels, cmap, plot_name = plot_name, scale = 100*12, 
                     Nrows = Nrows, Ncols = Ncols, plot_n = plot_n)

    if eg_cube is None: eg_cube = Obs
    if Obs is not None: plot_map(Obs, "Observations", 1)
    plot_map(Sim[0,:], "Simulation - 10%", Ncols - 1)
    plot_map(Sim[1,:], "Simulation - 90%", Ncols)
    

def evaluate_model(filename_out, dir_outputs, Obs, Sim, lmask, levels, cmap):    
    ax = plt.subplot(2, 3, 4)
    BayesScatter(Obs, Sim, lmask,  0.000001, 0.000001, ax)
    plot_model_maps(Sim, lmask, levels, cmap, Obs, Nrows = 2, Ncols = 3)
    
    
    X = Obs.data.flatten()[lmask]
    ncells = int(len(X)/Obs.shape[0])
    X = X.reshape([Obs.shape[0], ncells])
    Y = [Sim[i].data.flatten()[lmask].reshape([Obs.shape[0],ncells]) \
         for i in range(Sim.shape[0])]
    Y = np.array(Y)
   
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
    
    plt.gcf().set_size_inches(12, 8)
    
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
        model_title -- name of model run. Used as directory and filename.
        trace_file -- netcdf filename containing trace (produced in pymc_MaxEnt_train.py)
        y_filen -- filename of dependant variable (i.e burnt area)
        x_filen_list -- filanames of independant variables
            (ie bioclimate, landscape metrics etc)
        months_of_year --- which months to extact on training and projecting
        dir_outputs -- where stuff gets outputted
        dir_projecting -- The directory of the data used for prjections. 
            Should contain same files for independant varibales as dir_training 
            (though you should be able to adpated this easily if required). 
            Doesnt need dependant variable, but if there, this will (once
            we've implmented it) attempt some evaluation.
        sample_for_plot -- how many iterations (samples) from optimixation should be used 
            for plotting and evaluation.
        levels -- levels on the colourbar on observtation and prodiction maps
        cmap -- levels on the colourbar on observtation and prodiction maps
    Returns:
         (to be added too)
    """

    """ 
        SETPUT 
    """
    ### input data paths and filenames
    model_title = 'simple_example_model'
    
    trace_file = "outputs//trace-trees_consec_dry_mean_crop_pas_humid_totalVeg-frac_points_0.01-Month_7-nvariables_-frac_random_sample0.01-nvars_6-niterations_100.nc"
    scaler_file = "outputs//scalers-trees_consec_dry_mean_crop_pas_humid_totalVeg-frac_points_0.01-Month_7-nvariables_-frac_random_sample0.01-nvars_6-niterations_100.csv"
   
    y_filen = "GFED4.1s_Burned_Fraction.nc"
    CA_filen = None
    x_filen_list=["trees.nc","consec_dry_mean.nc",
                  "crop.nc", "pas.nc", "humid.nc", "totalVeg.nc"] 
    
    
    months_of_year = [7]
    
    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/"

    sample_for_plot = 20

    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    cmap = 'OrRd'

    
    subset_function = sub_year_months
    subset_function_args = {'months_of_year': months_of_year}

    """ 
        RUN evaluation 
    """
    evaluate_MaxEnt_model(trace, y_filen, x_filen_list, scalers, CA_filen, dir_projecting,
                         dir_outputs, model_title, filename,
                         subset_function, subset_function_args,
                         sample_for_plot, 
                         run_evaluation = run_evaluation, run_projection = run_projection,
                         grab_old_trace = grab_old_trace,
                         levels = levels, cmap = cmap)
    
    
