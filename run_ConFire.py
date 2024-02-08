from train import *
from evaluate import *


from   io     import StringIO
import numpy  as np
import cftime

import matplotlib.pyplot as plt


if __name__=="__main__":
    namelist = 'namelists/ConFire_example.txt'
    control_Direction = [1, -1, 1, -1]
    
    trace, scalers, training_namelist = \
                    train_MaxEnt_model_from_namelist(namelist)

    params = read_variables_from_namelist(training_namelist)
    output_dir = params['dir_outputs']
    output_file = params['filename_out']
    control_direction = read_variables_from_namelist(params['other_params_file'])
    control_direction = control_direction['control_Direction']
    control_names = read_variables_from_namelist(namelist)['control_names']

    def call_eval(control_run_name, extra_params = None, run_only = True, *args, **kw):
        return evaluate_MaxEnt_model_from_namelist(training_namelist, namelist,
                                                    run_only = run_only, 
                                                    control_run_name = 'Standard_fuel',
                                                    extra_params = extra_params)
    def Standard_limitation(controlID):
        
        control_Directioni = control_Direction
        control_Directioni[-controlID] = 0.0
        extra_params = {"control_Direction": control_Directioni}
        
        return call_eval('Standard_'+ str(controlID), extra_params)

    Control = call_eval('control', run_only = False)
    Standard = [Standard_limitation(i) for i in range(len(control_Direction))]

    def make_time_series(cube, name):
        try: 
            cube.coord('latitude').guess_bounds()
        except:
            pass

        try:
            cube.coord('longitude').guess_bounds()
        except:
            pass
        grid_areas = iris.analysis.cartography.area_weights(cube)
        area_weighted_mean = cube.collapsed(['latitude', 'longitude'], 
                                            iris.analysis.MEAN, weights=grid_areas)
        out_file = figName + '/points-' + name + '.csv'
        np.savetxt(out_file, area_weighted_mean.data, delimiter=',')
        
        TS = area_weighted_mean.collapsed('realization', 
                                          iris.analysis.PERCENTILE, percent=[10, 90])
        time_coord = TS.coord('time')
        time_datetime = time_coord.units.num2date(time_coord.points)
        time_datetime = cftime.date2num(time_datetime, 'days since 0001-01-01 00:00:00')/365.24
        TS = np.append(time_datetime[:,None], np.transpose(TS.data), axis = 1)

        out_file = figName + '/time_series' + name + '.csv'
        np.savetxt(out_file, TS, delimiter=',', header = "year,p10%,p90%")
        return TS
    
    figName = output_dir + 'figs/' + output_file + '-control_TS'
    makeDir(figName + '/')
    control_TS = make_time_series(Control[0], 'Control')
    standard_TS = np.array([make_time_series(cube[0], name) \
                            for cube, name in zip(Standard, control_names)])
    set_trace()
    colors = [(0, 1, 0, 0.5), (0, 0, 1, 0.5), (1, 0, 0, 0.5), (0, 0, 0, 0.5)]  # Green, Blue, Red, Black with 50% opacity
    
    fig, ax = plt.subplots()
    
    for i in range(standard_TS.shape[0]):
        ax.fill_between(range(standard_TS.shape[2]), 
                        standard_TS[i, 0, :], standard_TS[i, 1, :], color=colors[i])

    # Set labels and title
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Polygons with Color Vision Impaired Color Combination')
    
    
