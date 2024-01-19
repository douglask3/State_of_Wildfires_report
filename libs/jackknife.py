import numpy as np
import matplotlib.pyplot as plt

from response_curves import *
from pymc_extras import *
from BayesScatter import *

import numpy as np
import matplotlib.pyplot as plt

#trying stuff out
from sklearn.metrics import roc_auc_score
from sklearn.utils import shuffle


def jackknife(x_filen_list, X, Sim, fig_dir, response_grouping = None, **common_args):
    
    diff = {}
    contributions = {}
    
    if response_grouping is not None:
        group_contributions = {}
        for group_index, group in enumerate(response_grouping):
            for col in group:
                Xi = X.copy()
                Xi[:, col] = 0.0
            varname = f"group_{group_index}"
            makeDir(varname)

            print("Plotting", varname)

            Sim2 = runSim_MaxEntFire(X=Xi, **common_args, run_name=varname + "/deleted", test_eg_cube=True)

            diff[varname] = np.abs(non_masked_data(Sim) - non_masked_data(Sim2))
            contributions[varname] = np.mean(diff[varname] / non_masked_data(Sim))

            group_contributions[varname] = contributions[varname]

    else:
        for col in range(X.shape[1]):
            Xi = X.copy()
            Xi[:, col] = 0.0
            varname = x_filen_list[col]
            if varname.endswith(".nc"):
                varname = varname[:-3]
            makeDir(varname)

            print("Plotting", varname)

            Sim2 = runSim_MaxEntFire(X=Xi, **common_args, run_name=varname + "/deleted", test_eg_cube=True)
            
            diff[varname] = np.abs(np.mean(non_masked_data(Sim)) - np.mean(non_masked_data(Sim2)))
            contributions[varname] = (diff[varname] / np.mean(non_masked_data(Sim)))
            
            #diff[varname] = np.mean(np.abs(non_masked_data(Sim[1]) - non_masked_data(Sim2[1])))
            #contributions[varname] = (diff[varname] / np.mean(non_masked_data(Sim[1])))

    #set_trace()  
            #np.mean(np.abs(non_masked_data(Sim) - non_masked_data(Sim2)))
    #set_trace()
    sorted_contributions = {k: v for k, v in sorted(contributions.items(), key=lambda item: item[1], reverse=True)}
    sorted_variable_names = list(sorted_contributions.keys())
    sorted_values = list(sorted_contributions.values())
    
    print("Plotting Jackknife")
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 12})

    fig_jackknife = plt.figure(figsize=(10, 6))
    ax = fig_jackknife.add_subplot(111)

    # Using terracotta color (#c86558) for the bars
    bars = ax.barh(sorted_variable_names, sorted_values, color='#c86558', edgecolor='black', linewidth=0.8)
    
    #adding values at the end of each bar
    for i, v in enumerate(sorted_values):
        ax.text(v + 0.1, i, str(round(v, 2)), color='black', va='center', fontsize=9)
    
    ax.set_xlabel('Contribution (%)', fontname='Times New Roman', fontsize=14)
    #ax.set_title('Variables Contribution', fontname='Times New Roman', fontsize=16, pad=20)  # Adding a title with some padding
    ax.grid(axis='x', linestyle='-', linewidth=0.5, alpha=0.7)  # Adjusting grid appearance
    ax.invert_yaxis()
    
    #removing borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Annotating arrow-like tips at the end of the left spine     
    ax.annotate('', xy=(0, 1.01), xytext=(0, -0.01), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))
            
    # Annotating arrow-like tips at the end of the bottom spine 
    ax.annotate('', xy=(1, 0), xytext=(-0.01, 0), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))

    # Customizing tick labels and formatting contributions
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
    #ax.set_xticks(np.arange(0, max(sorted_values) + 1, 5))
    
    plt.show()
    set_trace()
    
    figure_filename = fig_dir +'jackknife'
    figure_dir =  combine_path_and_make_dir(figure_filename)
    
    # Save the figure as a PNG file
    fig_jackknife.savefig(figure_filename + '.png')
    plt.clf()
    
def jackknife_plotting(Sim, jackknife_type, X, x_filen_list, fig_dir, **common_args):

    figure_filename = fig_dir + jackknife_type
    figure_dir =  combine_path_and_make_dir(figure_filename)
    
    print("Plotting : " + jackknife_type)  
    
    if jackknife_type == "full":
        jackknife_FUN = jackknife_full
    elif jackknife_type == "percentile":
        jackknife_FUN = jackknife_for_percentile
    else:
        set_trace()
    
    contributions = {}
    
    for col in range(X.shape[1] - 1):
        Xi = X.copy()
        Xi[:, col] = 0.0
        varname = x_filen_list[col] 
        if varname.endswith(".nc"):
            varname = varname[:-3]
        makeDir(varname)
        
        Sim, Sim2 = jackknife_FUN(x_filen_list, X, Sim, varname, col, **common_args)
    
        contributions[varname] = \
            np.mean(np.abs(non_masked_data(Sim) - non_masked_data(Sim2))) * 100
    
    
    sorted_contributions = {k: v for k, v in sorted(contributions.items(), key=lambda item: item[1], reverse=True)}
    sorted_variable_names = list(sorted_contributions.keys())
    sorted_values = list(sorted_contributions.values())
    
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 12})

    fig_jackknife = plt.figure(figsize=(10, 6))
    ax = fig_jackknife.add_subplot(111)

    # Using terracotta color (#c86558) for the bars
    bars = ax.barh(sorted_variable_names, sorted_values, color='#c86558', edgecolor='black', linewidth=0.8)
    
    #adding values at the end of each bar
    for i, v in enumerate(sorted_values):
        ax.text(v + 0.1, i, str(round(v, 2)), color='black', va='center', fontsize=9)
    
    ax.set_xlabel('Contribution (%)', fontname='Times New Roman', fontsize=14)
    #ax.set_title('Variables Contribution', fontname='Times New Roman', fontsize=16, pad=20)  # Adding a title with some padding
    ax.grid(axis='x', linestyle='-', linewidth=0.5, alpha=0.7)  # Adjusting grid appearance
    ax.invert_yaxis()
    
    #removing borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Annotating arrow-like tips at the end of the left spine     
    ax.annotate('', xy=(0, 1.01), xytext=(0, -0.01), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))
            
    # Annotating arrow-like tips at the end of the bottom spine 
    ax.annotate('', xy=(1, 0), xytext=(-0.01, 0), xycoords='axes fraction', textcoords='axes fraction',
            arrowprops=dict(arrowstyle='->', linewidth=1.5, color='black'))

    # Customizing tick labels and formatting contributions
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
    ax.set_xticks(np.arange(0, max(sorted_values) + 1, 5))
    
    plt.show()
    #set_trace()
    
    figure_filename = fig_dir +'jackknife'
    plt.show()
    figure_dir =  combine_path_and_make_dir(figure_filename)
    '''
    # Save the figure as a PNG file
    fig_jackknife.savefig(figure_filename + '.png')
    plt.clf()
    '''
    
            #contributions[varname] = \
        #    np.abs((non_masked_data(Sim) - non_masked_data(Sim2)))*100
    
    #total_contribution = sum(np.abs(value) for value in contributions.values())
    
    # Calculate contributions as a percentage of the total
    #contributions_percentage = {key: (np.abs(value) / total_contribution) * 100 for key, value in contributions.items()}
    
    
    # Calculate the total sum of absolute contributions
    #total_sum = sum(np.abs(value) for value in contributions.values())

    # Compute the percentage contributions
    #contributions_percentage = {
    #    varname: np.abs(value) / total_sum * 100 
    #    for varname, value in contributions.items()
    #}
    
    '''    
    # Extract values for plotting
    #data = list(contributions.values())
    data = list(contributions_percentage.values())

    # Calculate percentiles and medians for each variable
    percentiles_10 = np.percentile(data, 10, axis=1)
    percentiles_90 = np.percentile(data, 90, axis=1)
    medians = np.median(data, axis=1)
    
    #set_trace()
    # Create violin plot
    plt.figure(figsize=(10, 6))
    plt.violinplot(data, vert=False)

    # Plotting markers and error bars for percentiles
    for i, var in enumerate(contributions.keys()):
        plt.scatter(medians[i], i + 1, color='red')
        plt.errorbar(medians[i], i + 1, xerr=[[medians[i] - percentiles_10[i]], [percentiles_90[i] - medians[i]]],
                     fmt='none', ecolor='black')

    plt.xlabel('Contributions')
    plt.ylabel('Variables')
    plt.title('Contributions of Variables')
    plt.yticks(np.arange(1, len(contributions) + 1), list(contributions.keys()))
    plt.tight_layout()
    plt.show()
    set_trace()
    
    
    '''
