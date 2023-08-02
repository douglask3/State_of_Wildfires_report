import matplotlib.pyplot as plt    
import numpy as np

def BayesScatter(X, Y, logXmin = None, logYmin = None, ax = None):
    
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
