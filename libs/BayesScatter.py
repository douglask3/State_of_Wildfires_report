import matplotlib.pyplot as plt    
import numpy as np
from pdb import set_trace
import iris

def BayesScatter(X, Y, lmask = None, logXmin = None, logYmin = None, ax = None):
    pc_size = 10
    percentiles = np.arange(pc_size, 100, pc_size)    
    if lmask is not None:
        X = X.data.flatten()[lmask]
        Y = Y.collapsed('realization', iris.analysis.PERCENTILE, percent = percentiles)
        Y = np.array([Y[i].data.flatten()[lmask] for i in range(Y.shape[0])]).T
    else:
        Y = np.percentile(Y, q = percentiles, axis = 0).T
    

    if logXmin is not None: X = logXmin + X
    if logYmin is not None: Y = logYmin + Y
    if ax is None: fig, ax = plt.subplots()
    
    ncols = int(Y.shape[1]/2)
    line_widths = np.linspace(0.2, 2, ncols)

    for i in range(ncols):        
        ax.vlines(X, ymin = Y[:,i], ymax = Y[:, -i-1], 
                  linewidth = line_widths[i],
                  alpha = 0.01, color = 'black')
    
    vals = [min(np.min(X), np.min(Y)),  max(np.max(X), np.max(Y))]
    ax.plot(vals, vals, color = 'blue', linestyle = '--')
    ax.scatter(X, Y[:, ncols + 1], s = 0.2, color = 'red')
    if logYmin is not None: ax.set_yscale('log')
    if logXmin is not None: ax.set_xscale('log')

    plt.xlabel("Observation")
    plt.ylabel("Simulation")
