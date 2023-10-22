from pymc_MaxEnt import *

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
        levels -- levels on teh colourbar on observtation and prodiction maps
        cmap -- levels on teh colourbar on observtation and prodiction maps
    Returns:
        trace file, maps, etc (to be added too)
    """

    """ 
        SETPUT 
    """
    """ optimization """

    model_title = 'Global_model-linRC'

    dir_training = "~/ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Global/historic_TS_2000_2019/obsclim/"
    y_filen = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/ISIMIP3a_obs/GFED5_Burned_Percentage.nc"

    x_filen_list=["crop.nc", "csoil.nc", "cveg.nc", "rhumid.nc",
                  "vpd.nc", "tas.nc", 
                  "lightn.nc", "pas.nc", "soilM.nc", "precip.nc",
                  "popDens.nc"]

    grab_old_trace = False
    cores = 1
    fraction_data_for_sample = 0.01
    niterations = 100

    year_range = [2002, 2004]
    extent = [-180, 180, -59.5, 84]
    months_of_year = range(0, 12)

    """ Projection/evaluating """
    dir_outputs = 'outputs/'

    dir_projecting = "~/ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5/Global/historic_TS_1880_1889/obsclim/"

    sample_for_plot = 4

    levels = [0, 0.1, 1, 2, 5, 10, 20, 50, 100] 
    cmap = 'OrRd'

    run_evaluation = False
    run_projection = True
    
    """ 
        RUN optimization 
    """
    subset_function = [sub_year_range, sub_year_months, contrain_coords]
    subset_function_args = [{'year_range': year_range},
                            {'months_of_year': months_of_year},
                            {'extent': extent}]

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
    months_of_year = 1
    subset_function_args = [{'months_of_year': months_of_year},
                            {'extent': extent}]
    y_filen = "vpd.nc"
    predict_MaxEnt_model(trace, y_filen, x_filen_list, scalers, dir_projecting,
                         dir_outputs, model_title, filename,
                         subset_function, subset_function_args,
                         sample_for_plot, 
                         run_evaluation = run_evaluation, run_projection = run_projection,
                         grab_old_trace = grab_old_trace,
                         levels = levels, cmap = cmap)
    
