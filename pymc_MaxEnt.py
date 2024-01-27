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

    person = 'Maria'
    quick = False

    if person == 'Maria':
        #dir_training = "D:/Doutorado/Sanduiche/research/maxent-variables/2002-2011/"
        dir_training = "D:/Doutorado/Sanduiche/research/maxent-variables/variables_masked_na/fixed_gfed_reference/2002-2021/" 
    else:
        dir_training = "/data/users/dkelley/MaxEnt/driving_Data/2002-2021/"



    year_range = [2002, 2009]
    
    #x_filen_list= ["ed.nc", "consec_dry_mean.nc", "savanna.nc", "rhumid.nc",
    #               "lightn.nc", "popDen.nc", "forest.nc", "precip.nc",
    #               "pasture.nc", "cropland.nc", "grassland.nc", "np.nc",
    #               "tas_max.nc", "csoil.nc", "wetland.nc",
    #               "vpd.nc", "soilM.nc", "road_density.nc"]
    
    
    
    #x_filen_list= ["ed.nc", "consec_dry_mean.nc", "savanna.nc", "cveg.nc", "rhumid.nc",
    #               "lightn.nc", "popDen.nc", "forest.nc", "precip.nc",
    #               "pasture.nc", "cropland.nc", "grassland.nc", "np.nc",
    #               "tas_max.nc", "csoil.nc",
    #               "vpd.nc", "soilM.nc", "road_density.nc"]
    
    #"tas_mean.nc", "tca.nc", "mpa.nc"
                   
    #x_filen_list= ["tas_max.nc", "precip.nc", "ed.nc", "csoil.nc", 
    #                "road_density.nc", "forest.nc", "pasture.nc",
    #                "lightn.nc", "grassland.nc", "consec_dry_mean.nc", 
    #                "cropland.nc", "soilM.nc", "savanna.nc", "rhumid.nc"] 

    x_filen_list= ["ed.nc", "csoil.nc", "tas_max.nc", "precip.nc",
                    "road_density.nc", "forest.nc", "pasture.nc"]                    
                  
    

    if quick:

        model_title = 'model-test-7varlin-'
        biome_IDs = range(0,1)

        fraction_data_for_sample = 0.001
        min_data_points_for_sample = 1000 #minimum grid cells to use
        cores = 2
        months_of_year = [8, 9, 10]
        niterations = 100
    else:
        model_title = 'final-full-lin-pow-7-1000-'
        biome_IDs = range(1,7)
        fraction_data_for_sample = 0.2
        min_data_points_for_sample = 6000 #minimum grid cells to use
        cores = 5
       
        months_of_year = [8,9,10]
        niterations = 1000

    CA_filen = "brazil_NAT.nc"
    y_filen = "Area_burned_NAT"

    #CA_filen = "brazil_NON.nc"
    #y_filen = "Area_burned_NON"

    #CA_filen = None
    #y_filen = "GFED4.1s_Burned_Fraction"

    model_title = model_title + y_filen
    y_filen = y_filen + '.nc'

    grab_old_trace = True # set to True till you get the code running. Then set to False when you start adding in new response curves
    
    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = dir_training
    
    sample_for_plot = 1000
    
    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    dlevels = [-20, -10, -5, -2, -1, -0.1, 0.1, 1, 2, 5, 10, 20]
    cmap = 'OrRd'
    dcmap = 'RdBu_r'

    #response_grouping = None
    
    #response_grouping = [["tas_max.nc", "consec_dry_mean.nc"], ["savanna.nc", "grassland.nc", "pasture.nc"], ["csoil.nc", "forest.nc"],
    #                    [ "precip.nc", "rhumid.nc", "soilM.nc", "lightn.nc"],
    #                    ["ed.nc", "road_density.nc", "cropland.nc"]]
    
    #response_grouping = [["tas_max.nc", "precip.nc", "consec_dry_mean.nc", "rhumid.nc", "soilM.nc"], 
    #                    ["road_density.nc", "lightn.nc"], ["forest.nc", "pasture.nc",
    #                    "grassland.nc","cropland.nc", "savanna.nc", "ed.nc", "csoil.nc"]]
    
    #response_grouping= [["ed.nc", "tca.nc", "np.nc", "mpa.nc"], ["consec_dry_mean.nc", 
    #                    "precip.nc", "tas_max.nc", "tas_mean.nc", "vpd.nc", "rhumid.nc",
    #                    "soilM.nc"], ["savanna.nc", "forest.nc", "pasture.nc", "grassland.nc",
    #                    "cropland.nc"], ["cveg.nc", "csoil.nc"], ["lightn.nc", "popDens.nc", 
    #                    "road_density.nc"]] 
    
    response_grouping = [["tas_max.nc", "precip.nc"], ["ed.nc", 
                    "road_density.nc"], ["forest.nc", "pasture.nc", "csoil.nc"]] 
     
    """ 
        RUN optimization 
        
    """
    for biome_ID in biome_IDs:
        subset_function = [sub_year_range, 
                            sub_year_months, constrain_BR_biomes]
        subset_function_args = [{'year_range': year_range},
                            {'months_of_year': months_of_year},
                            {'biome_ID': [biome_ID]}]

        
        filename = '_' +  str(len(x_filen_list)) + \
                '-frac_points_' + str(fraction_data_for_sample) + \
                '-Month_' +  '_'.join([str(mn) for mn in months_of_year]) + \
               f'-Biome_{biome_ID}'

        #### Optimize
        trace, scalers, variable_info_file = \
                    train_MaxEnt_model(y_filen, x_filen_list, CA_filen , dir_training, 
                                      filename, dir_outputs,
                                      fraction_data_for_sample,
                                      subset_function, subset_function_args,
                                      niterations, cores, model_title,
                                      '/biome' + str(biome_ID),
                                      grab_old_trace,
                                      min_data_points_for_sample = min_data_points_for_sample)

        """ 
            RUN projection 
       
        """
        evaluate_MaxEnt_model_from_namelist(variable_info_file, subset_function_args = [{'year_range': [2010,2019]},
                                            {'months_of_year': months_of_year},
                                            {'biome_ID': [biome_ID]}], 
                                            dir = dir_projecting,
                                            grab_old_trace = grab_old_trace,
                                            sample_for_plot = sample_for_plot,
                                            levels = levels, cmap = cmap,
                                            dlevels = dlevels, dcmap = dcmap,
                                            response_grouping = response_grouping)
        
