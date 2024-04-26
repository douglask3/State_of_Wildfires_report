
from train import *
from evaluate import *


from   io     import StringIO
import numpy  as np
import cftime

import matplotlib.pyplot as plt

def call_eval(training_namelist, namelist,
              control_run_name, extra_params = None, run_only = True, *args, **kw):
    return evaluate_MaxEnt_model_from_namelist(training_namelist, namelist,
                                               run_only = run_only, 
                                               control_run_name = control_run_name,
                                               extra_params = extra_params, *args, **kw)
    
   
def Standard_limitation(training_namelist, namelist,
                        controlID, name, control_direction, *args, **kws):   
    control_Directioni = np.array(control_direction.copy())
    control_Directioni[:] = 0.0
    control_Directioni[controlID] = control_direction[controlID]
    
    extra_params = {"control_Direction": control_Directioni}
    
    return call_eval(training_namelist, namelist,
                     name + '/Standard_'+ str(controlID), extra_params, hyper = False,
                     *args, **kws)
    
def make_time_series(cube, name, figName):
    try: 
        cube.coord('latitude').guess_bounds()
    except:
        pass
    
    try:
        cube.coord('longitude').guess_bounds()
    except:
        pass
    
    cube.data = np.ma.masked_invalid(cube.data)
    grid_areas = iris.analysis.cartography.area_weights(cube)
    #set_trace()
    area_weighted_mean = cube.collapsed(['latitude', 'longitude'], iris.analysis.WPERCENTILE, percent = 90, weights=grid_areas)
    
    #area_weighted_mean = area_weighted_mean.aggregated_by('year', iris.analysis.MAX)
    out_file = figName + '/points-' + name + '.csv'
    np.savetxt(out_file, area_weighted_mean.data, delimiter=',')
    
    TS = area_weighted_mean.collapsed('realization', 
                                      iris.analysis.PERCENTILE, percent=[25, 75])
    
    time_coord = TS.coord('time')
    time_datetime = time_coord.units.num2date(time_coord.points)
    time_datetime = cftime.date2num(time_datetime, 'days since 0001-01-01 00:00:00')/365.24
    TS = np.append(time_datetime[:, None], np.transpose(TS.data), axis = 1)
    
    out_file = figName + '/time_series' + name + '.csv'
    np.savetxt(out_file, TS, delimiter=',', header = "year,p25%,p75%")
    return TS

def run_experiment(training_namelist, namelist, control_direction, control_names, 
                   output_dir, output_file, 
                   name = '', *args, **kws):
    if "baseline" in name: 
        run_only = False
    else:
        run_only = True
    print("running: " + name)
    name = name + '-'

    temp_file = 'temp/run_ConFire_lock' + (output_dir + output_file + name).replace('/', '_') + '.txt'
    #if os.path.isfile(temp_file): return None

    Control, Y, X, lmask, scalers  = call_eval(training_namelist, namelist,
                        name + '/control', run_only = run_only, return_inputs = True,
                        *args, **kws)
    
    Standard = [Standard_limitation(training_namelist, namelist,i, name, control_direction, 
                                    *args, Y = Y, X = X, lmask = lmask, scalers = scalers,
                                    **kws) \
                for i in range(len(control_direction))]

    figName = output_dir + 'figs/' + output_file + '-' + name + 'control_TS'
    makeDir(figName + '/')
    control_TS = make_time_series(Control[0], 'Control', figName)
    standard_TS = np.array([make_time_series(cube[0], name, figName) \
                           for cube, name in zip(Standard, control_names)])
    
        
    open(temp_file, 'a').close() 

def run_ConFire(namelist):   
    
    run_info = read_variables_from_namelist(namelist) 

    regions = run_info['regions']
    subset_function_args = run_info['subset_function_args']
    

    for region in regions:
        model_title = run_info['model_title'].replace('<<region>>', region)
        dir_training = run_info['dir_training'].replace('<<region>>', region)
        
        subset_function_args['months_of_year'] = run_info['region_months'][region]

        trace, scalers, training_namelist = \
                        train_MaxEnt_model_from_namelist(namelist, model_title = model_title,
                                                         dir_training = dir_training,
                                                         subset_function_args = subset_function_args)
        
        params = read_variables_from_namelist(training_namelist)
        output_dir = params['dir_outputs']
        output_file = params['filename_out']
        
        control_direction = read_variables_from_namelist(params['other_params_file'])
        control_direction = control_direction['control_Direction']
        control_names = read_variables_from_namelist(namelist)['control_names']
    
        def find_replace_period_model(exp_list):
            exp_list_all = [item.replace('<<region>>', region) for item in exp_list \
                            if "<<experiment>>" not in item and "<<model>>" not in item]
            looped_items = [item for item in exp_list \
                            if "<<experiment>>" in item and "<<model>>" in item]
            for experiment, period in zip(experiments, periods):
                for model in models:
                    dirs = [item.replace("<<period>>", period) for item in looped_items]    
                    dirs = [item.replace("<<model>>", model) for item in dirs]   
                    dirs = [item.replace("<<experiment>>", experiment) for item in dirs] 
                    dirs = [item.replace('<<region>>', region) for item in dirs] 
                    exp_list_all += dirs
             
            return exp_list_all
    
        dir_projecting = run_info['dir_projecting'].replace('<<region>>', region)
        experiment_dirs  = run_info['experiment_dir']
        experiment_names = run_info['experiment_names']
        experiments = run_info['experiment_experiment']
        periods = run_info['experiment_period']
        models = run_info['experiment_model']
        experiment_dirs = find_replace_period_model(experiment_dirs)
        experiment_names = find_replace_period_model(experiment_names)
        
        y_filen = run_info['x_filen_list'][0]
    
        run_experiment(training_namelist, namelist, control_direction, control_names,
                                  output_dir, output_file, 'baseline', 
                                  model_title = model_title, 
                                  subset_function_args = subset_function_args)
        
        [run_experiment(training_namelist, namelist, control_direction, 
                                     control_names,
                                     output_dir, output_file, name, dir = dir, 
                                     y_filen = y_filen, model_title = model_title,
                                     subset_function_args = subset_function_args) \
                          for name, dir in zip(experiment_names, experiment_dirs)]


if __name__=="__main__":
    namelist = 'namelists/ConFire_Canada.txt'
    namelist = 'namelists/Greece.txt'
    namelist = 'namelists/tuning.txt'
    #namelist = 'namelists/SOW2023.txt'
    
    run_ConFire(namelist)
    set_trace()
    #experiment_TS = np.array([make_time_series(cube[0], name)  for cube, name in zip(experiment, experiment_names)])
    
    

    
    
    colors = [(0, 1, 0, 0.5), (0, 0, 1, 0.5), (1, 0, 0, 0.5), (0, 0, 0, 0.5)]  # Green, Blue, Red, Black with 50% opacity
    
    fig, ax = plt.subplots()
    
    for i in range(standard_TS.shape[0]):
        ax.fill_between(range(standard_TS.shape[2]), 
                        standard_TS[i, 0, :], standard_TS[i, 1, :], color=colors[i])

    # Set labels and title
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Polygons with Color Vision Impaired Color Combination')
    
    
