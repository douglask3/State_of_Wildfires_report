import numpy as np
import matplotlib.pyplot as plt


plt.rcParams.update({'font.size': 14, 'font.family': 'Times New Roman'})

'''
#def f(x, x0, p):
    """
    Calculate f(x, x0, p) = ln|e^(x-x0)^p + 1|.
    
    Arguments:
        x -- The input variable.
        x0 -- The parameter x0.
        p -- The parameter p.

    """
    
 #   y = (np.exp(x - x0)) ** p
 #   y = np.log(y + 1)
 #   return y

#def linear_combined(x, beta, gama, x0, p):
    """
    Calculate f(x, x0, p) = ln|e^(x-x0)^p + 1|.
    
    Arguments:
        x -- The input variable.
        x0 -- The parameter x0.
        p -- The parameter p.

    """
    #y = (beta_lin * x) + beta_comb * (np.log(((np.exp(x - x0)) ** p) + 1))
    #return y

# Generate x values
x_values = np.linspace(-10, 10, 1000)

# Define different combinations x0 and p
parameters = [(1, 5)]

# Plot the curves for different combinations of x0 and p
#plt.figure(figsize=(10, 6))
for x0, p in parameters:
    y_values = f(x_values,x0, p)
    #plt.plot(x_values, y_values, label=f"x0={x0}, p={p}")
    plt.plot(x_values, y_values)

#plt.rcParams.update({'font.size': 14}) 
plt.xlabel("Explanatory variables")
plt.ylabel("Target variable")
plt.legend()
plt.grid(True)
plt.show() 
    
''' 


#def linear_curve(a, b, x): 
#    return a*x + b
    
def power_response_curve(X):  
    return 2**X


#def X2_response_curve(self, X):  
#    return (X - self)**2.0

x_values = np.linspace(0, 10, 1000)

# Define different combinations self
#parameters = [(0.1, 0.5)] #linear

#parameters = [(0.5)] 

# Plot the curves for different combinations of x0 and p
#plt.figure(figsize=(10, 6))
#for self in parameters:
y_values = power_response_curve(x_values)
    #plt.plot(x_values, y_values, label=f"x0={x0}, p={p}")
plt.plot(x_values, y_values)
    

plt.xlabel("Explanatory variables")
plt.ylabel("Target variable")
plt.legend()
plt.grid(True)
plt.show()

