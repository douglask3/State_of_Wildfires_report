import numpy as np

# Generate synthetic data for FWI
np.random.seed(0)  # For reproducibility

# Define the parameters for generating synthetic data
slope = 0.2
intercept = 25
noise_level = 2

# Generate years from 1960 to 2013
years = np.arange(1960, 2014)

# Generate synthetic FWI data with linear trend and added noise
t = years - 2022  # shift years to be relative to 2022
fwi_trend = slope * t + intercept
fwi_noise = np.random.normal(loc=0, scale=noise_level, size=len(years))
fwi = fwi_trend + fwi_noise

data = np.column_stack((years, fwi))
np.savetxt('synthetic_fwi_data.csv', data, delimiter=',', header='Year,FWI', comments='')

# Now you have synthetic FWI data in the array 'fwi'
