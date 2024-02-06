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
        set_trace()
        call_eval('Standard_'+ str(controlID), extra_params)

    Control = call_eval('control', run_only = False)
    Standard = [Standard_limitation(i) for i in range(len(control_Direction))]
    
    set_trace()
    
