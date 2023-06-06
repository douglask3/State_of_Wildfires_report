import numpy as np
import matplotlib.pyplot as plt
import os
from   io     import StringIO
import numpy  as np
import math

import pymc  as pm
import pytensor
import pytensor.tensor as tt

from pdb import set_trace as browser
class MaxEntFire(object):
    """
    Maximum Entropy fire model which takes indepedant variables and coefficants. 
    At the moment, just a linear model fed through a logistic function to convert to 
    burnt area/fire probablity. But we'll adapt that.  
    """ 
    def __init__(self, betas, inference = False):
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
        
        self.betas = betas

    def fire_model(self, X):
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
        y = self.numPCK.dot(X, self.betas)

        BA = 1.0/(1.0 + self.numPCK.exp(-y))
    
        return BA
'''     
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
   

#calling the function

x1 = np.round(np.linspace(0, 1, num=10), decimals=1)

hin1 = hinge_1(0.5, 0.7, 2, 3)

print(hin1)
 
#notes
 
# np.linspace is a function from NumPy that generates a sequence of evenly spaced numbers within a specified range
'''
'''
def hinge_2(x, a1, b1, a2, b2):

    
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

#calling the function

x = np.linspace(0, 1, num=100)
print(x)
hin2 = hinge_2(x, -1, 1, 2, 0.5)

print("y = ", hin2)
plt.plot(x, hin2)
plt.show()
'''
'''
def exp(x, p):
    
    y = x**p
    
    return y
 

#calling the function

x = np.linspace(0, 1, num=100)
e = exp(x,2)

print(e)
plt.plot(x, e)
plt.show()
'''



