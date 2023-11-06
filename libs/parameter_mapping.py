import matplotlib.pyplot as plt
import arviz as az
from pymc_extras import *

import numpy as np
import pandas as pd
import seaborn as sns

from pdb import set_trace

def plot_basic_parameter_info(traces, fig_dir):
    """ plots parameter distributions from a trace file.
    More info to follow
    """
    az.plot_trace(traces)
    plt.savefig(fig_dir + 'traces.png')



def paramter_map(traces, x_filen_list, fig_dir):
    params, parameter_names = select_post_param(traces)

    def remove_nc(varname):
        if varname.endswith(".nc"):
            varname = varname[:-3]
        return varname
    variable_names = [remove_nc(varname) for varname in x_filen_list]


    # Create an empty dictionary to store data
    data_dict = {}
    
    # Iterate through the array_list and assign column names accordingly
    for i, array in enumerate(params):
        if array.ndim == 2:
            # If it's a 2D array, use column names that combine variable and parameter names
            
            param = parameter_names[i]
            col_names = [f"{var}_{param}" for var in variable_names]
        else:
            # If it's a 1D array, use column names with just parameter names
            col_names = [parameter_names[i]]
        
        # Convert the array to a DataFrame and assign column names
        df = pd.DataFrame(array, columns=col_names)
        # Add the DataFrame to the data_dict
        data_dict[f"instance_{i + 1}"] = df


    # Concatenate the DataFrames in data_dict to create the final DataFrame
    
    final_df = pd.concat(data_dict, axis=1)

    set_trace()
