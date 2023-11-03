import numpy as np
import iris
import iris.analysis
import iris.coord_categorisation
from pdb import set_trace


annual_cube = iris.load_cube("../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/Forest.nc")

print(annual_cube)

monthly_cube = iris.load_cube("../ConFIRE_attribute/isimip3a/driving_data/GSWP3-W5E5-20yrs/Brazil/AllConFire_2000_2009/GFED4.1s_Burned_Fraction.nc")


numbers = monthly_cube.coord('time').points
annual_cube = annual_cube.interpolate([('time', numbers)], iris.analysis.Linear())
    
print(annual_cube)

