import numpy as np
import math
from pymc_extras import *
from plot_maps import *

import matplotlib.pyplot as plt

from pdb import set_trace

def non_masked_data(cube):
    if type(cube) is tuple: cube = cube[0]
    out = cube.data[cube.data.mask == False].data
    return out[~np.isnan(out)]


def initial_curve_experiment(Sim, Xi, col_to_keep, name, trace, sample_for_plot, 
                              eg_cube, lmask, *args, **kw):
    Sim1, Sim2 = standard_curve_experiment(Sim, Xi, col_to_keep, name, trace, sample_for_plot, 
                              eg_cube, lmask, *args, **kw)
    X = Xi.copy()
    X[:,:] = 0.0
    Sim1 = runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, 
                             "/all_to_0", *args, **kw)
    #set_trace()
    return Sim1, Sim2


def standard_curve_experiment(Sim, Xi, col_to_keep, name, trace, sample_for_plot, 
                              eg_cube, lmask, *args, **kw):
    X = Xi.copy()
    other_cols = np.arange(X.shape[1])  # Create an array of all columns
    other_cols = other_cols[other_cols != col_to_keep]  # Exclude col_to_keep
    
    
    X[:, other_cols] = 0.0 
    
    Sim2 = runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, 
                             name + "/all_but_to_0", *args, **kw)

    
    return Sim, Sim2

def potential_curve_experiment(Sim, Xi, col_to_go, name, trace, sample_for_plot, 
                              eg_cube, lmask, *args, **kw):
    X = Xi.copy()
    X[:, col_to_go] = 0.0
        
    Sim2 = runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, 
                             name + "/to_0", *args, **kw)
    return Sim, Sim2


def sensitivity_curve_experiment(Sim, Xi, col, name, trace, sample_for_plot, 
                              eg_cube, lmask, *args, **kw):
    dx = 0.001
    X = Xi.copy()
    X[:, col] -= dx/2.0  # Subtract 0.1 of all values for the current column
    
    Sim1 = runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, 
                             name + "/subtract_" + str(dx/2.0), *args, **kw)
        
    X = Xi.copy() #restore values
    X[:, col] += dx/2.0 #add 0.1 to all values for the current column
        
    Sim2 = runSim_MaxEntFire(trace, sample_for_plot, X, eg_cube, lmask, 
                             name + "/add_" + str(dx/2.0), *args, **kw)
    return Sim1, Sim2


def response_curve(Sim, curve_type, trace, sample_for_plot, X, eg_cube, lmask, 
                   dir_samples, fig_dir, grab_old_trace, x_filen_list, 
                   levels, cmap, dlevels, dcmap, scalers = None, *args, **kw):  

    figure_filename = fig_dir + curve_type + '-response'
    figure_dir =  combine_path_and_make_dir(figure_filename + '-maps/') 
    
    print("Plotting reponse curve: " + curve_type)
    
    map_type = 1
    if curve_type == "initial":
        response_FUN = initial_curve_experiment
        map_type = 0
    elif curve_type == "standard":
        response_FUN = standard_curve_experiment
    elif curve_type == "potential":
        response_FUN = potential_curve_experiment
    elif curve_type == "sensitivity":
        response_FUN = sensitivity_curve_experiment
        map_type = 2
    else:
        set_trace()

    #fig_map, ax_map = plt.subplots(len(x_filen_list) + 1, 4* map_type)
    fig_curve = plt.figure()
    fig_map = plt.figure()

    fcol = math.floor(math.sqrt(X.shape[1]))
    frw = math.ceil(X.shape[1]/fcol)
    
    Ncol = 6
    if map_type == 0: Ncol = 3
    if map_type == 2: Ncol = 7
    
    def plotFun(cube, ylab = '', plot0 = 0, lvls = levels, cm = cmap, **kw2): 
        if map_type > -1:
            plot_BayesModel_maps(cube, lvls, cm, ylab = ylab,
                                 Nrows = len(x_filen_list) + 1, Ncols = Ncol, plot0 = plot0, 
                                 colourbar = True, fig = fig_map, **kw2)

    #def dplotFun(cube, *args, **kw): plotFun(cube, dlevels, dcmap, *args, **kw)
    
    plotFun(Sim, 'Control')
    for col in range(X.shape[1]-1):
        
        varname = x_filen_list[col]
        if varname.endswith(".nc"):
            varname = varname[:-3]
        makeDir(varname)
        Sim1, Sim2 = response_FUN(Sim, X, col, varname, trace, sample_for_plot, 
                                      eg_cube, lmask, dir_samples, grab_old_trace)
        
        plotN = Ncol * (col + 1)
        if map_type > 0: plotFun(Sim2, varname, plotN)
        if map_type == 2:
            plotFun(Sim1, '', plotN + 2, figure_filename = figure_dir + varname + '-absolute')

        if Sim1 is not None:
            diff = Sim2.copy()
            diff.data = Sim2.data - Sim1.data
        else:
            diff = Sim2
        
        diffP = diff.collapsed('time', iris.analysis.MEAN) 
        
        plotFun(diffP, '', plotN + 2 * map_type, dlevels, dcmap, 
                figure_filename = figure_dir + varname + '-difference')    

        ax = fig_curve.add_subplot(frw,fcol, col + 1)  # Select the corresponding subplot
        
        variable_name = x_filen_list[col].replace('.nc', '')
        ax.set_title(variable_name)
        
        num_bins = 10
        
        hist, bin_edges = np.histogram(X[:, col], bins=num_bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        median_values = []
        percentile_10 = []
        percentile_90 = []
        
        for i in range(num_bins):
            
            mask = (X[:, col] >= bin_edges[i]) & (X[:, col] < bin_edges[i + 1])
            if np.any(mask):
                values_in_bin = []
        
                for rw in range(Sim2.shape[0]):
                    values_in_bin.append(non_masked_data(diff[rw])[mask])
                values_in_bin = np.array(values_in_bin).flatten()    
                   
                median_values.append(np.median(values_in_bin))
                percentile_10.append(np.percentile(values_in_bin, 10))
                percentile_90.append(np.percentile(values_in_bin, 90))
            else:
                median_values.append(np.nan)
                percentile_10.append(np.nan)
                percentile_90.append(np.nan)
        
        if scalers is not None:
            bin_centers = bin_centers*(scalers[1, col] - scalers[0, col]) + scalers[0, col]
        ax.plot(bin_centers, median_values, marker='.', label='Median')
        ax.fill_between(bin_centers, percentile_10, percentile_90, alpha=0.3, 
                        label='10th-90th Percentiles')                           
    
          
    if map_type > -1:
        fig_map.set_size_inches(12, 4*X.shape[1])
        fig_map.tight_layout()
        fig_map.subplots_adjust(left=0.15)
        fig_map.savefig(figure_filename + '-maps.png')
    plt.close(fig_map)
    
    fig_curve.savefig(figure_filename + '-curves.png')   
    plt.close(fig_curve)
    plt.clf()


