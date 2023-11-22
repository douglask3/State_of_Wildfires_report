import numpy as np
import matplotlib.pyplot as plt

from response_curves import *
from pymc_extras import *
from evaluate import *


import numpy as np
import matplotlib.pyplot as plt

def jackknife(x_filen_list, X, **common_args):
    
    contributions = {}

    for col in range(X.shape[1] - 1):
        
        varname = x_filen_list[col] 
        if varname.endswith(".nc"):
            varname = varname[:-3]
        makeDir(varname)

        X = np.delete(X, col, axis=1)
        
        Sim2 = runSim_MaxEntFire(X = X, **common_args, run_name=varname + "deleted", test_eg_cube=True)

        contributions[varname] = np.mean(np.abs(non_masked_data(Sim) - non_masked_data(Sim2))) * 100

    sorted_contributions = {k: v for k, v in sorted(contributions.items(), key=lambda item: item[1], reverse=True)}
    sorted_variable_names = list(sorted_contributions.keys())
    sorted_values = list(sorted_contributions.values())

    # Create the bar plot
    plt.figure(figsize=(10, 6))
    plt.barh(sorted_variable_names, sorted_values, color='blue')
    plt.xlabel('Contribution (%)')
    plt.title('Variable Contributions')
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.gca().invert_yaxis()
    plt.show()
    set_trace()
    return sorted_contributions

# Usage example:
# jackknife(x_filen_list, **common_args)

