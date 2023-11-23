import numpy as np
import matplotlib.pyplot as plt

from response_curves import *
from pymc_extras import *


import numpy as np
import matplotlib.pyplot as plt

def jackknife(x_filen_list, X, fig_dir, **common_args):
    
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
        
        contributions[varname] = np.mean(np.abs(non_masked_data(Sim) - non_masked_data(Sim2))) * 100

    sorted_contributions = {k: v for k, v in sorted(contributions.items(), key=lambda item: item[1], reverse=True)}
    sorted_variable_names = list(sorted_contributions.keys())
    sorted_values = list(sorted_contributions.values())
    
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 12})

    fig_jackknife = plt.figure(figsize=(10, 6))
    ax = fig_jackknife.add_subplot(111)

    # Using terracotta color (#c86558) for the bars
    bars = ax.barh(sorted_variable_names, sorted_values, color='#c86558', edgecolor='black', linewidth=0.8)
    
    for i, v in enumerate(sorted_values):
        ax.text(v + 0.1, i, str(round(v, 2)), color='black', va='center', fontsize=9)
    
    ax.set_xlabel('Contribution (%)', fontname='Times New Roman', fontsize=14)
    #ax.set_title('Variables Contribution', fontname='Times New Roman', fontsize=16, pad=20)  # Adding a title with some padding
    ax.grid(axis='x', linestyle='-', linewidth=0.5, alpha=0.7)  # Adjusting grid appearance
    ax.invert_yaxis()
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ax.annotate('', xy=(0, 1.01), xytext=(0, -0.01), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))
            
    # Annotating arrow-like tips at the end of the bottom spine facing towards the right
    ax.annotate('', xy=(1, 0), xytext=(-0.01, 0), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))

    # Customizing tick labels and formatting contributions
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.set_xticks(np.arange(0, max(sorted_values) + 1, 5))
    
    plt.show()
    set_trace()
    
    '''
    # Create the bar plot
    print("Plotting Jackknife")
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 12})
    fig_jackknife = plt.figure(figsize=(10, 6))
    figure_filename = fig_dir +'jackknife'
    figure_dir =  combine_path_and_make_dir(figure_filename)
    plt.barh(sorted_variable_names, sorted_values, color='#c86558')
    plt.xlabel('Contribution (%)')
    #plt.title('Variables Contribution')
    plt.grid(axis='x', linestyle='--', alpha=0.6)
    plt.gca().invert_yaxis()
    
    plt.show()
    set_trace()

    # Save the figure as a PNG file
    fig_jackknife.savefig(figure_filename + '.png')
    plt.clf()
    color='#FF7F50'
    '''
