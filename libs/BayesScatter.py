import matplotlib.pyplot as plt    
import numpy as np
from pdb import set_trace
import iris

def sqidge_01s(X, logXmin):
    if logXmin is not None: X = logXmin + X*(1 - 2*logXmin)
    return(X)

def Bayes_benchmark(Y, X, lmask, logXmin = None, logYmin = None, ax = None):
    X = X.data.flatten()[lmask]
    Y = np.array([Y[i].data.flatten()[lmask] for i in range(Y.shape[0])]).T
    set_trace()

def BayesScatter(X, Y, lmask = None, logXmin = None, logYmin = None, ax = None):
    pc_size = 10
    percentiles = np.arange(pc_size, 100, pc_size)    
    if lmask is  None: 
        Y = np.percentile(Y, q = percentiles, axis = 0).T       
    else:
        X = X.data.flatten()[lmask]
        Y = Y.collapsed('realization', iris.analysis.PERCENTILE, percent = percentiles)
        Y = np.array([Y[i].data.flatten()[lmask] for i in range(Y.shape[0])]).T
    
    if len(X) > 1000:
        select = np.random.choice(len(X), 1000)
        X = X[select]
        Y = Y[select, :]
    
    X = sqidge_01s(X, logXmin)
    Y = sqidge_01s(Y, logYmin)
    if ax is None: fig, ax = plt.subplots()
    
    ncols = int(Y.shape[1]/2)
    line_widths = np.linspace(0.2, 2, ncols)
    alpha = 1.0/Y.shape[1]

    for i in range(ncols):        
        ax.vlines(X, ymin = Y[:,i], ymax = Y[:, -i-1],
                  linewidth = line_widths[i],
                  alpha = alpha, color = 'black')
    
    vals = [min(np.min(X), np.min(Y)),  max(np.max(X), np.max(Y))]
    ax.plot(vals, vals, color = 'blue', linestyle = '--')
    ax.scatter(X, Y[:, ncols + 1], s = 0.2, color = 'red')
    if logYmin is not None: ax.set_yscale('logit')
    if logXmin is not None: ax.set_xscale('logit')

    plt.xlabel("Observation")
    plt.ylabel("Simulation")
