from train import *
from evaluate import *


from   io     import StringIO
import numpy  as np

import matplotlib.pyplot as plt


if __name__=="__main__":
    namelist = 'namelists/ConFire_example.txt'

    trace, scalers, training_namelist = \
                    train_MaxEnt_model_from_namelist(namelist)

    def call_eval(control_run_name, extra_params = None, *args, **kw):
        return evaluate_MaxEnt_model_from_namelist(training_namelist, namelist, 
                                                    control_run_name = 'Standard_fuel',
                                                    extra_params = extra_params)
        
    Sim = call_eval('control')
    set_trace()
    extra_params = {"control_Direction": [0, -1, 1, -1]}
