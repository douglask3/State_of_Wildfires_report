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
    Maximum Entropy fire model which takes indepedant variables and coefficants. 
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
       
        self.q = select_key_or_defualt(params, 'q', 0.0)
        self.betas = params['betas']
        self.powers = select_key_or_defualt(params, 'powers', None)
        self.x2s = select_key_or_defualt(params, 'x2s', None)
        

    def burnt_area(self, X):
        """calculated predicted burnt area based on indepedant variables. 
            At the moment, just a linear model fed through a logistic function to convert to 
            burnt area/fire probablity. But we'll adapt that.   
        Arguments:
	    X -- numpy or tensor 2d array of indepenant variables, each columne a different 
                    variable, no. columns (no. variables) is same as length of betas.
        Returns:
            numpy or tensor (depdaning on 'inference' option) 1 d array of length equal to 
	    no. rows in X of burnt area/fire probabilities.
        """
        def dot_fun(x, p): return self.numPCK.dot(x, p)
        
        y = dot_fun(X, self.betas)

        def add_response_curve(params, FUN, y):
            if params is not None:
                XR = FUN(X)
                y = y + dot_fun(XR, params[0,:]) 
            return(y)

        y = add_response_curve(self.powers, self.power_response_curve, y)
        y = add_response_curve(self.x2s, self.X2_response_curve, y)

        BA = 1.0/(1.0 + self.numPCK.exp(-y))
        
        return BA
    
    def burnt_area_uninflated(self, X):
        BA = self.burnt_area(X)
        if self.q == 0.0: return BA
        return BA / (1 + self.q * (1 - BA))
     
    def hinge_1(x0, y0, a, b):
        """ fits a hinge curve function
        x -- numpy array 
        -- hinge point
        """
    
        if np.all(x1 > x0):
            y = a*x1 + b   
        else:
            y = y0
    
        return y
   

    def hinge_2(self, x, a1, b1, a2, b2):
  
        x0 = (b2-b1)/(a1-a2)
        print("x0 = ", x0)
        y=[]
    
        for xi in x:
            if xi > x0:
                yi = a2*xi + b2   
            else:
                yi = a1*xi + b1
    
            y.append(yi)
        return y

#
    def power_response_curve(self, X):  
        return X**self.powers[1,:]

    def X2_response_curve(self, X):  
        return (X - self.x2s[1,:])**2.0



