import numpy as np
import iris
import iris.analysis
import iris.coord_categorisation
from pdb import set_trace


annual_cube = iris.load_cube("D:/Doutorado/Sanduiche/research/maxent-variables/2002-2011/Forest.nc")
iris.coord_categorisation.add_year(annual_cube, 'time')
iris.coord_categorisation.add_month_number(annual_cube, 'time', name='month')
print(annual_cube)

#monthly_cube = iris.load_cube("D:/Doutorado/Sanduiche/research/maxent-test/driving_and_obs_overlap/AllConFire_2000_2009/GFED4.1s_Burned_Fraction.nc")
#iris.coord_categorisation.add_year(monthly_cube, 'time')
#iris.coord_categorisation.add_month_number(monthly_cube, 'time', name='month')
#print(monthly_cube)



if np.all(annual_cube.coord("month").points == 1):
# Regrid the monthly data to the same grid as the monthly example cube
    numbers = np.tile(np.arange(1, 13), len(annual_cube.coord("time").points))
    annual_cube = annual_cube.interpolate([('time', numbers)], iris.analysis.Linear())
    
print(annual_cube)
#print(monthly_cube.coord('month').points)
#set_trace()




#print(numbers)

#print(dataset_monthly_regridded)
