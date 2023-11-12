import sys
sys.path.append('fire_model/')
sys.path.append('libs/')

from train import *
from evaluate import *

import os
from   io     import StringIO
import numpy  as np

import matplotlib.pyplot as plt

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

    person = 'Doug'

    if person == 'Maria':
        model_title = 'Example_model-biome-6'
        #dir_training = "/gws/nopw/j04/jules/mbarbosa/driving_and_obs_overlap/AllConFire_2000_2009/"
        dir_training = "D:/Doutorado/Sanduiche/research/maxent-variables/2002-2011/"

        #y_filen = "GFED4.1s_Burned_Fraction.nc"
        y_filen = "Area_burned_NAT.nc"
        #y_filen = "Area_burned_NON.nc"
        
        CA_filen = "brazil_NAT.nc"
        #CA_filen = "brazil_NON.nc"
        biome_IDs = [6]
        
        x_filen_list=["consec_dry_mean.nc", "savanna.nc", "cveg.nc", "rhumid.nc",
                      "lightn.nc", "popDens.nc", "forest.nc", "precip.nc",
                      "crop.nc", "pas.nc", "grassland.nc", "ed.nc", "np.nc",
                      "tas_max.nc", "tas_mean.nc", "tca.nc", "te.nc", "mpa.nc",
                      "totalVeg.nc", "vpd.nc", "csoil.nc", "SoilM.nc"]

        cores = 2
        fraction_data_for_sample = 0.05
    else:
        model_title = 'Example_model-q-reduced'

        dir_training = "../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/"
        y_filen = "Area_burned_NAT.nc"
        CA_filen = "brazil_NAT.nc"
        biome_IDs = [1,2]
        x_filen_list=["trees.nc", "pr_mean.nc", "consec_dry_mean.nc", 
                  #"lightn.nc", "popDens.nc",
                  #"crop.nc", "pas.nc", 
                  "humid.nc", "csoil.nc", "tas_max.nc",
                  "totalVeg.nc"] 

        cores = 1
        fraction_data_for_sample = 0.01

    grab_old_trace = True # set to True till you get the code running. Then set to False when you start adding in new response curves
    
    niterations = 100
    min_data_points_for_sample = 1000
    months_of_year = [7]
    
    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = dir_training
    #dir_projecting = "D:/Doutorado/Sanduiche/research/maxent-variables/2012-2021/"
    sample_for_plot = 20
    
    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    dlevels = [-20, -10, -5, -2, -1, -0.1, 0.1, 1, 2, 5, 10, 20]
    cmap = 'OrRd'
    dcmap = 'RdBu_r'
     
    """ 
        RUN optimization 
        
    """
    for biome_ID in biome_IDs:
    
        subset_function = [sub_year_months, constrain_BR_biomes]  
        subset_function_args = [{'months_of_year': months_of_year}, {'biome_ID': [biome_ID]}]

        filename = '_'.join([file[:-3] for file in x_filen_list]) + \
                '-frac_points_' + str(fraction_data_for_sample) + \
                '-Month_' +  '_'.join([str(mn) for mn in months_of_year]) + \
                f'-Biome_{biome_ID}'
    

        #### Optimize
        trace, scalers, variable_info_file = \
                    train_MaxEnt_model(y_filen, x_filen_list, CA_filen, None, dir_training, 
                                      filename, dir_outputs,
                                      fraction_data_for_sample,
                                      subset_function, subset_function_args,
                                      niterations, cores, model_title, grab_old_trace,
                                      min_data_points_for_sample = min_data_points_for_sample)

        """ 
            RUN projection 
        """

        evaluate_MaxEnt_model_from_namelist(variable_info_file, dir = dir_projecting,
                                            grab_old_trace = grab_old_trace,
                                            sample_for_plot = sample_for_plot,
                                            levels = levels, cmap = cmap,
                                            dlevels = dlevels, dcmap = dcmap)
    
