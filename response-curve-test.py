import numpy as np
import matplotlib.pyplot as plt

'''def f(x, x0, p):
    """
    Calculate f(x, x0, p) = ln|e^(x-x0)^p + 1|.
    
    Arguments:
        x -- The input variable.
        x0 -- The parameter x0.
        p -- The parameter p.

    """
    
    y = (np.exp(x - x0)) ** p
    y = np.log(y + 1)
    return y
'''
def linear_combined(x, beta, gama, x0, p):
    """
    Calculate f(x, x0, p) = ln|e^(x-x0)^p + 1|.
    
    Arguments:
        x -- The input variable.
        x0 -- The parameter x0.
        p -- The parameter p.

    """
    y = (beta * x) + gama * (np.log(((np.exp(x - x0)) ** p) + 1))
    return y

# Generate x values
x_values = np.linspace(0, 1, 1000)

# Define different combinations of beta, gama, x0 and p
parameters = [(10, 50, 0.5, 20), (10, -2, 0.5, 20), (10,5,0.5, 20)]

# Plot the curves for different combinations of x0 and p
plt.figure(figsize=(10, 6))
for beta, gama, x0, p in parameters:
    y_values = linear_combined(x_values, beta, gama,x0, p)
    plt.plot(x_values, y_values, label=f"x0={x0}, p={p}, beta={beta}, gama={gama}")

plt.xlabel("x")
plt.ylabel("f(x, beta, gama, x0, p)")
plt.legend()
plt.grid(True)
plt.show()



