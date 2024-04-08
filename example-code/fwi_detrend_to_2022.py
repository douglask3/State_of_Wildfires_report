import numpy as np
import pandas as pd
import statsmodels.api as sm
from pdb import set_trace
import matplotlib.pyplot as plt

# Step 0; Load fwi data from CSV using pandas
df = pd.read_csv('synthetic_fwi_data.csv')

# Extract years and FWI values
years = df['Year'].values
fwi = df['FWI'].values

# Step 1: Fit a linear regression model
t = years - 2022  # shift years to be relative to 2022
X = sm.add_constant(t)  # add a constant term for intercept
model = sm.OLS(fwi, X)
results = model.fit()

# Step 2: Get the coefficients (slope and intercept)
fwi0, delta = results.params

# Step 3: Calculate the trend
trend = delta * t

# Step 4: Detrend the data
fwi_detrended = fwi - trend

# Step 5: Plotting
plt.figure(figsize=(10, 6))
plt.plot(years, fwi, label='FWI', color='blue')
plt.plot(years, fwi_detrended, label='Detrended & Shifted FWI', color='red')
plt.xlabel('Year')
plt.ylabel('FWI')
plt.title('Time Series Plot of FWI and Detrended & Shifted FWI')
plt.legend()
plt.grid(True)
plt.show()
set_trace()
# Now you have the detrended data shifted to 2022 in fwi_detrended_shifted
