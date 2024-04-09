import numpy as np
import pandas as pd
import statsmodels.api as sm
from pdb import set_trace
import matplotlib.pyplot as plt

# Step 0; Load fwi data from CSV using pandas
df_sim = pd.read_csv('synthetic_fwi_data-sim.csv')
df_obs = pd.read_csv('synthetic_fwi_data-obs.csv')

# Extract years and FWI values
years = df_sim['Year'].values
fwi_sim = df_sim['FWI'].values
fwi_obs = df_obs['FWI'].values

# Step 1a: Fit a linear regression model to obs and sim
t = years - 2022  # shift years to be relative to 2022
X = sm.add_constant(t)  # add a constant term for intercept
def find_regression_parameters(fwi):
    model = sm.OLS(fwi, X)
    results = model.fit()

    # Step 1b: Get the coefficients (slope and intercept)
    fwi0, delta = results.params
    
    return fwi0, delta, np.std(fwi - delta * t) 

fwi0_sim, delta_sim, std_sim =  find_regression_parameters(fwi_sim)
fwi0_obs, delta_obs, std_obs =  find_regression_parameters(fwi_obs)

# Step 2: Detrend the sim and scale to obs
fwi_detrended_sim = fwi0_obs + (fwi_sim - delta_sim * t - fwi0_sim) * (std_obs/ std_sim)

# Step 3: Plotting
plt.figure(figsize=(10, 6))
plt.plot(years, fwi_obs, label='Obs', color='blue')
plt.plot(years, fwi_sim, label='Sim', color='red')
plt.plot(years, fwi_detrended_sim, label='Detrended & Shifted FWI', color='purple')
plt.xlabel('Year')
plt.ylabel('FWI')
plt.title('Time Series Plot of FWI and Detrended & Shifted FWI')
plt.legend()
plt.grid(True)
plt.show()
set_trace()
# Now you have the detrended data shifted to 2022 in fwi_detrended_shifted
