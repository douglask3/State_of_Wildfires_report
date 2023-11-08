import numpy as np
import matplotlib.pyplot as plt
import os
from   io     import StringIO
import numpy  as np
import math

import pymc  as pm
import pytensor
import pytensor.tensor as tt

from pdb import set_trace

def select_key_or_defualt(dirc, key, default):
    if key in dirc:
        return dirc[key]
    else:
        return default

class MaxEntFire(object):
    """
    Maximum Entropy fire model which takes independent variables and coefficients. 
    At the moment, just a linear model fed through a logistic function to convert to 
    burnt area/fire probablity. But we'll adapt that.  
    """ 
    def __init__(self, params, inference = False):
        """
        Sets up the model based on betas and repsonse curve pararameters (response curve 
            not yet implmented
        Arguments:
	    X -- numpy or tensor 2d array of indepenant variables, each columne a different 
                    variable, no. columns (no. variables) is same as length of betas.
	    inference -- boolean.   If True, then used in bayesian inference and uses 
                                        tensor maths. 
			            If False, used in normal mode or prior/posterior sampling 
                                        and uses numpy.
        Returns:
            callable functions. 'fire_model' (below) is the main one to use.
        """
        self.inference = inference
        if self.inference:
            self.numPCK =  __import__('pytensor').tensor
        else:
            self.numPCK =  __import__('numpy')
       
        
        self.lin_betas = params['lin_betas']
        self.lin_beta_constant = select_key_or_defualt(params, 'lin_beta_constant', 0.0)
        self.pow_betas = select_key_or_defualt(params, 'pow_betas', None)
        self.pow_power = select_key_or_defualt(params, 'pow_power', None)
        self.x2s_betas = select_key_or_defualt(params, 'x2s_betas', None)
        self.x2s_X0    = select_key_or_defualt(params, 'x2s_X0'   , 0.0 )
        self.q = select_key_or_defualt(params, 'q', 0.0)
        self.comb_betas = select_key_or_defualt(params, 'comb_betas', None)   
        self.comb_X0 = select_key_or_defualt(params, 'comb_X0', None) 
        self.comb_p = select_key_or_defualt(params, 'comb_p', None)
        #Maria: add your response curve parameter selection thing
        

    def burnt_area(self, X):
        """calculated predicted burnt area based on indepedant variables. 
            At the moment, just a linear model fed through a logistic function to convert to 
            burnt area/fire probablity. But we'll adapt that.   
        Arguments:
	    X -- numpy or tensor 2d array of indepenant variables, each columne a different 
                    variable, no. columns (no. variables) is same as length of betas.
        Returns:
            numpy or tensor (depending on 'inference' option) 1 d array of length equal to 
	    no. rows in X of burnt area/fire probabilities.
        """
        
        y = self.numPCK.dot(X, self.lin_betas)

        def add_response_curve(Rbetas, FUN, y):
            if Rbetas is not None:
                XR = FUN(X)
                y = y + self.numPCK.dot(XR, Rbetas) 
            return(y)

        y = add_response_curve(self.pow_betas, self.power_response_curve, y)
        y = add_response_curve(self.x2s_betas, self.X2_response_curve   , y)
        y = add_response_curve(self.comb_betas, self.linear_combined_response_curve , y)
        
        # y = add_response_curve(paramers, function, y)
        # Maria: add yours here 

        BA = 1.0/(1.0 + self.numPCK.exp(-y))
        
        return BA
    
    def burnt_area_spread(self, X):
        BA = self.burnt_area(X)
        if self.q == 0.0: return BA
        return BA / (1 + self.q * (1 - BA))
     
    def power_response_curve(self, X):  
        return self.pow_power**X


    def X2_response_curve(self, X):  
        return (X - self.x2s_X0)**2.0

    def linear_combined_response_curve(self, X):
        return np.log(((np.exp(X - self.comb_X0)) ** self.comb_p) + 1)
