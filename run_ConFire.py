from train import *
from evaluate import *


from   io     import StringIO
import numpy  as np

import matplotlib.pyplot as plt


if __name__=="__main__":
    namelist = 'namelists/ConFire_example.txt'
    control_Direction = [1, -1, 1, -1]
    
    trace, scalers, training_namelist = \
                    train_MaxEnt_model_from_namelist(namelist)

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

    def make_time_series(cube):
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

        return area_weighted_mean.collapsed('realization', 
                                            iris.analysis.PERCENTILE, percent=[10, 90])
    
    
    control_TS = make_time_series(Control[0]).data
    standard_TS = [make_time_series(cubes[0]).data for cubes in Standard]

    colors = ['green', 'blue', 'red', 'black']  
    fig, ax = plt.subplots()

    for i in range(standard_TS.shape[0]):
        ax.fill_between(range(standard_TS.shape[2]), 
                        standard_TS[i, 0, :], standard_TS[i, 1, :], color=colors[i])

    # Set labels and title
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_title('Polygons with Color Vision Impaired Color Combination')
    set_trace()
    
