import numpy as np
import iris

import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from pdb import set_trace


class ConFire(object):
    def __init__(self, params, inference = False):
        """
        Initalise parameters and calculates the key variables needed to calculate burnt area.
        """
        self.inference = inference
        if self.inference:
            self.numPCK =  __import__('pytensor').tensor
        else:
            self.numPCK =  __import__('numpy')
        
        self.params = params
        
    def burnt_area(self, X, return_controls = False, return_limitations = False):
        ## finds controls        
        def cal_control(cid = 0):
            cid = 0
            ids = self.params['controlID'][cid]
            betas = self.params['betas'][cid] * self.params['driver_Direction'][cid]

            return self.numPCK.dot(X[:,ids], betas)
            

        controls = [cal_control(i) for i in range(len(self.params['controlID']))]
        if return_controls: return controls

        cdir = self.params['control_Direction']

        def sigmoid(y, k):
            return 1.0/(1.0 + self.numPCK.exp(-y * k))

        limitations = [sigmoid(y, k) for y, k in zip(controls, cdir)]
        if return_limitations: return limitations
    
        return self.numPCK.prod(limitations, axis = 0)
    
    
    def emc_weighted(self, emc, precip, wd_pg):
        
        try:
            wet_days = 1.0 - self.numPCK.exp(-wd_pg * precip)
            emcw = (1.0 - wet_days) * emc + wet_days
        except:
            emcw = emc.copy()
            emcw.data  = 1.0 - self.numPCK.exp(-wd_pg * precip.data)
            emcw.data = emcw.data + (1.0 - emcw.data) * emc.data
        return(emcw)


      

