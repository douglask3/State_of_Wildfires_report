import numpy as np
import matplotlib.pyplot as plt

from response_curves import *
from pymc_extras import *
from evaluate import *


import numpy as np
import matplotlib.pyplot as plt

def jackknife(x_filen_list, X, **common_args):

    
    contributions = {}
    
    Sim = runSim_MaxEntFire(X = X, **common_args, run_name = "control", test_eg_cube = True)

    for col in range(X.shape[1] - 1):
        
        varname = x_filen_list[col] 
        if varname.endswith(".nc"):
            varname = varname[:-3]
        makeDir(varname)
        
        X = X.copy()
        X[:, col] = 0.0 
        
        Sim2 = runSim_MaxEntFire(X = X, **common_args, run_name=varname + "deleted", test_eg_cube=True)
        
        #return Sim, Sim2
        
        contributions[varname] = np.mean(np.abs(non_masked_data(Sim) - non_masked_data(Sim2))) * 100

    sorted_contributions = {k: v for k, v in sorted(contributions.items(), key=lambda item: item[1], reverse=True)}
    sorted_variable_names = list(sorted_contributions.keys())
    sorted_values = list(sorted_contributions.values())

    # Create the bar plot
    fig_jackknife = plt.figure(figsize=(10, 6))
    fig_jackknife = plt.barh(sorted_variable_names, sorted_values, color='blue')
    plt.xlabel('Contribution (%)')
    plt.title('Variable Contributions')
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.gca().invert_yaxis()
    fig_jackknife.savefig('jackknife.png')
    set_trace()
    #plt.show()
    return sorted_contributions
    
    #fig_jackknife.savefig('jackknife.png')   
    #plt.close(fig_jackknife)
    #plt.clf()

