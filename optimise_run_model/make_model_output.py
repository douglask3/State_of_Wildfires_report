import sys
sys.path.append('../')

import warnings
warnings.filterwarnings('ignore')

import os
from   io     import StringIO
import numpy  as np
import pandas as pd
import csv

import iris
import matplotlib.pyplot as plt
import numpy.ma as ma

import cartopy.crs as ccrs
from   libs.plot_maps    import *
from scipy.stats import norm
from scipy.stats import skewnorm
from scipy.stats import lognorm
import random
from pdb import set_trace as browser


def mkDir(dir):
    try: os.mkdir(dir)  
    except: pass

def npLogit(x):
    return np.log(x/(1.0-x))

# We're using the row which has minimum sigma
def newCubes3D(variable, step, eg_cube_in, dimname = 'model_level_number', minV = -10, maxV = None):
    
    
    #maxV = np.max(eg_cube_in.data)    
    #maxV = np.min([1.0, np.ceil(maxV) * 1.5])
    
    def newCube(i):
        coord = iris.coords.AuxCoord(i, dimname)
        eg_cube = eg_cube_in.copy()
        eg_cube.data[eg_cube.data > 0.0] = 0.0
        try:
            eg_cube.remove_coord(dimname)
        except:
            pass
        eg_cube.add_aux_coord(coord)
        return(eg_cube)
    
    if not hasattr(step, '__len__'):
        
        if minV is None: minV = np.round(npLogit(np.min(eg_cube_in.data[eg_cube_in.data>0.0]))) #* 2.0
        if maxV is None: maxV = -minV#np.round(npLogit(np.max(eg_cube_in.data[eg_cube_in.data<1.0])))
        
        minV = minV - 1
        step = np.arange(minV ,maxV, step)     
    
    eg_cubes = iris.cube.CubeList([newCube(i) for i in step])
    eg_cubes = eg_cubes.merge()[0]
    
    return(eg_cubes)

################
## input info ##
################

dir = "../ConFIRE_ISIMIP/inputs2/"
model_names = os.listdir(dir)

model_names = [model_names[3]]
experiments = os.listdir(dir + model_names[0])
browser()
experiments = [experiments[0]]
files = {'soilwMax'           : 'soil12.nc',
         'shallow_soilw'      : 'soilM_top.nc',
         'deep_soilw'         : 'soilM_bottom.nc',
         'precip'             : 'precip.nc',
         'emc'                : 'humid.nc',
         'treeCover'          : 'trees.nc',
         'pasture'            : 'pas.nc',
         'cropland'           : 'crop.nc',
         'vegcover'           : 'totalVeg.nc',
         'burnt_area'         : '../../../burnt_area_GFED4sObs.nc'}

param_file = '../ConFIRE_ISIMIP/outputs/params-for_sampling/'
title_output = 'attempt4-full'

###############
## open data ##
##############
def loadInputsParams(model, experiment):
    print(model)
    print(experiment)
    input_data = {}
    for key, file in files.items():
        try:
            data = iris.load_cube(dir + model + '/' + experiment + '/' + file)
        except:
            browser()
        input_data[key] = data

    params = pd.read_csv(param_file + model + '.csv')
    return input_data, params
    
inputs = [[loadInputsParams(model, experiment) for experiment in experiments] for model in model_names]

###########
## model ##
###########

class ConFIRE(object):
    def __init__(self, data, params):
        """
        Initalise parameters and calculates the key variables needed to calculate burnt area.
        """
        self.params = params


        ## finds controls
        self.fuel = self.control_fuel(data['vegcover'], data['soilwMax'], self.params['fuel_pw'],
                                      self.params['fuel_pg'])
        
        self.emcw = self.emc_weighted(data['emc'], data["precip"], self.params['wd_pg'])
        
        self.moisture = self.control_moisture(data['shallow_soilw'], data['deep_soilw'],
                                              self.emcw, data['treeCover'],
                                              self.params['cMs'], self.params['cM'], self.params['cMT'], 
                                              self.params['kM'], self.params['pT'])

        self.ignitions = self.control_ignitions(data['pasture'])

        self.suppression = self.control_suppression(data['cropland'])

        ## calculates limiting factor of each control.
        self.standard_fuel        = self.sigmoid(self.fuel       ,
                                            self.params[       'fuel_x0'], self.params[       'fuel_k'])  
        
        self.standard_moisture    = self.sigmoid(self.moisture   ,
                                            self.params[   'moisture_x0'], -self.params[   'moisture_k'])
        self.standard_ignitions   = self.sigmoid(self.ignitions   ,
                                            self.params[  'ignition_x0'], self.params[  'ignition_k'])
        self.standard_suppression = self.sigmoid(self.suppression,
                                            self.params['suppression_x0'], -self.params['suppression_k'])
        
        
        self.error = self.params['sigma']
        ## burnt area us just limitation of each control muliplied together.
        self.burnt_area_mode = self.standard_fuel * self.standard_moisture * self.standard_ignitions * \
            self.standard_suppression * self.params['max_f']
        
        ## find the mean burnt area
        self.burnt_area_calPDF(data, self.params['p0'], self.params['pp'])
        
        self.burnt_area = self.burnt_area_mean.copy() #* (1.0-self.p0)
        
        #browser()
        
        self.standard_moisture    = self.standard_moisture    / self.sigmoid(0.0, self.params['moisture_x0'],
                                                 -self.params['moisture_k'])
        self.standard_suppression = self.standard_suppression / self.sigmoid(0.0, self.params['suppression_x0'],
                                                 -self.params['suppression_k'])

        self.potential_fuel = self.potential(self.standard_fuel, "potential_fuel")
        self.potential_moisture = self.potential(self.standard_moisture, "potential_moisture")
        self.potential_ignitions = self.potential(self.standard_ignitions, "potential_ignitions")
        self.potential_suppression = self.potential(self.standard_suppression, "potential_suppression")

        self.sensitivity_fuel = self.sensitivity(self.fuel, self.params['fuel_x0'], self.params['fuel_k'],
                                    self.standard_fuel, "sensitivity_fuel")

        self.sensitivity_moisture = self.sensitivity(self.moisture, self.params['moisture_x0'], -self.params['moisture_k'],
                                    self.standard_moisture, "sensitivity_moisture")

        self.sensitivity_ignitions = self.sensitivity(self.ignitions, self.params['ignition_x0'], self.params['ignition_k'],
                                    self.standard_ignitions, "sensitivity_ignitions")


        self.sensitivity_suppression = self.sensitivity(self.suppression, self.params['suppression_x0'], -self.params['suppression_k'] ,
                                    self.standard_suppression, "sensitivity_suppression")


        ## if the inputs are iris cubes, we can add some useful metadata
        try:
            self.burnt_area.long_name = "burnt_area"
            self.burnt_area_mode.long_name = "burnt_area_mode"
            #self.burnt_area_median.long_name = "burnt_area_median"
            self.burnt_area_mean.long_name = "burnt_area_mean"
            self.burnt_area_pdf.long_name = "burnt_area_pdf"        

            self.fuel.long_name = "fuel continuity"
            self.fuel.units = '1'

            self.moisture.long_name = "moisture content"
            self.moisture.units = '1'

            self.ignitions.long_name = "ignitions"
            self.ignitions.units = 'km-2'

            self.suppression.long_name = "suppression"
            self.suppression.units = '1'

            self.standard_fuel.long_name = "standard_fuel"
            self.standard_moisture.long_name = "standard_moisture"
            self.standard_ignitions.long_name = "standard_ignitions"
            self.standard_suppression.long_name = "standard_suppression"

            self.standard_fuel.units = '1'
            self.standard_moisture.units = '1'
            self.standard_ignitions.units = '1'
            self.standard_suppression.units = '1'
        except:
            pass        
    
        
    def control_fuel(self, vegcover, alphaMax, fuel_pw, fuel_pg):
        """
        Definition to describe fuel load: while return the input; capability to be modified later.
        """
        return (vegcover**(fuel_pw+1)) * (fuel_pg * (alphaMax-1) + 1) / (1 + fuel_pg)
    
    def emc_weighted(self, emc, precip, wd_pg):
        try:
            wet_days = 1 - np.exp(-wd_pg * precip)
            emcw = wet_days + (1-wet_days) * emc
        except:
            emcw = emc.copy()
            emcw.data  = 1 - np.exp(-wd_pg * precip.data)
            emcw.data = emcw.data + (1-emcw.data) * emc.data
        return(emcw)

    def control_moisture(self, shallow_soil, deeps_soil, emc, treeCover, cMs, cM, cMT, kM, pT):
        """
        Definition to describe moisture
        """
        moist = (shallow_soil/100.0 + cMs * deeps_soil/100.0 + cM*emc + cMT * (treeCover**pT)) / (1 + cMs + cM + cMT)
        moist.data = 1 - np.log(1 - moist.data*kM)
        return moist


    def control_ignitions(self,pasture):
        """
        Definition for the measure of ignition
        """
        return pasture


    def control_suppression(self, cropland):
        """
        Definition for the measure of fire supression
        """
        return cropland

        """
        Defines potential limitation for each control in turn
        """
   
    def sensitivity(self, x, x0, k, fi, long_name = None):

        gradient = self.gradient(x, x0, k)
        sens = gradient * self.control_removal(fi)

        try: sens.units = '1'
        except: pass

        if long_name is not None:
            try: sens.long_name = long_name
            except: browser()
        return sens


    def control_removal(self, fi):
        return self.burnt_area_mode/fi


    def potential(self, fi, long_name = None):
        out = fi.copy()
        out.data = self.burnt_area_mode.data * ((1/out.data) - 1)

        try: out.units = '1'
        except: pass

        if long_name is not None:
            try: out.long_name = long_name
            except: pass
        return out


    def gradient(self, x, x0, k, dx = 0.0001):

        f1 = self.sigmoid(x + dx, x0, k)
        f2 = self.sigmoid(x - dx, x0, k)

        n1 = self.sigmoid(x0 + dx, x0, k)
        n2 = self.sigmoid(x0 - dx, x0, k)

        try:
            f1 = (f1 - f2)/(n1 - n2)
        except:
            f1.data = (f1.data - f2.data) / (n1.data - n2.data)

        return f1


    def sigmoid(self, x, x0,k):
        """
        Sigmoid function to describe limitation using tensor
        """
        try:
            out = x.copy()
            out.data = -k*(x.data - x0)
            out.data = 1.0/(1.0 + np.exp(out.data))
            x = out
        except:
            x = -k * (x - x0)
            x = 1.0/(1.0 + np.exp(x))
        return x
    
    def burnt_area_calPDF(self, data, p0, pp):
        
        mask = np.logical_not(self.burnt_area_mode.data.mask)
        self.burnt_area_pdf = newCubes3D('burnt_area', 0.5, data['burnt_area'])
        
        self.burnt_area_mean = self.burnt_area_mode.copy()
        self.burnt_area_mean.data[mask] = 0.0
        
        level_no = self.burnt_area_pdf.coord('model_level_number').points
        
        dist = norm(npLogit(self.burnt_area_mode.data[mask]), self.error)
        
        self.pz = 1.0 - (self.burnt_area_mode.data[mask]**pp) * (1.0 - p0)
        
        self.burnt_area_pdf.data[0][mask] = self.pz 
        x = self.burnt_area_pdf.coord('model_level_number').points
        for k in range(1, self.burnt_area_pdf.shape[0]):       
            self.burnt_area_pdf.data[k][mask] = dist.pdf(x[k]) * (1.0 - self.pz)
            self.burnt_area_mean.data[mask]  = self.burnt_area_mean.data[mask] +  dist.pdf(x[k]) * (1.0 - self.pz)  *(1/(1+np.exp(-x[k])))
                                     
        PDFtot = self.burnt_area_pdf.collapsed(['model_level_number'], iris.analysis.SUM)
        
        self.burnt_area_mean.data = self.burnt_area_mean.data/ PDFtot.data
       
        
###################
## seupup inputs ##
###################

input_data = inputs[0][0][0].copy()
for key in input_data.keys(): input_data[key] = inputs[0][0][0][key][0]
inputs[0][0][0]["burnt_area"].shape


###############
## bootstrap ##
##############


def bootSamples(x, samples, nsample = 2):
    def bootSample(x):
        #browser()or np.ma.is_masked(test)
        test = np.ma.sum(x)
        out = np.repeat(0, nsample)
        if test==0 : return out
        try:
            out = np.random.choice(samples, nsample, p =x/np.sum(x))
        except:
            out = np.repeat(0, nsample)
        return out
   
    x = np.apply_along_axis(bootSample, 0, x)
    
    return(x)



def runModelExperiment(input, paramLoc = None, nmns = None, nyr = None):
    input_data = input[0].copy()
    params = input[1]
    if paramLoc is None: paramLoc = params["sigma"].idxmin()
    print("running")
    
    if nyr is None:
        mns = range(0, input[0]["precip"].shape[0])
    else:
        mns = np.random.choice(np.int(input[0]["precip"].shape[0]/12), 12*nyr)*12
        mns += (np.arange(0, 12, 1/nyr)*nyr).astype(int)
        mns = np.sort(mns)
    
    if nmns is not None and nmns < len(mns): mns = mns[0:nmns]
    
    nmns = len(mns)
    
    for mn in mns:
        print(mn)
        print(mn/nmns) 
        for key in input_data.keys():
            try:
                input_data[key] = input[0][key][mn]
            except:
                #print(key)
                try:
                    input_data[key] = input[0][key][mn-12]
                except:
                    browser()
        model = ConFIRE(input_data, params.loc[paramLoc])
        
        PDF = model.burnt_area_pdf.data
        if mn == mns[0]:
            mask = PDF[0,:,:].mask == False
            BAs = model.burnt_area_pdf.coord("model_level_number").points
            BAs = 1/(1+np.exp(-BAs))
            BAs[0] = 0.0
            
        PDF = bootSamples(PDF[:,mask], BAs, 10)
        
        if mn == mns[0]:
            model_out = model
            PDF_out = PDF
        else:
            model_out.burnt_area_mean.data += model.burnt_area_mean.data
            model_out.fuel.data += model.fuel.data
            model_out.moisture.data += model.moisture.data
            model_out.ignitions.data += model.ignitions.data
            model_out.suppression.data += model.suppression.data
            
            model_out.standard_fuel.data += model.standard_fuel.data 
            model_out.standard_moisture.data += model.standard_moisture.data 
            model_out.standard_ignitions.data += model.standard_ignitions.data 
            model_out.standard_suppression.data += model.standard_suppression.data 
            
            model_out.potential_fuel.data += model.potential_fuel.data 
            model_out.potential_moisture.data += model.potential_moisture.data 
            model_out.potential_ignitions.data += model.potential_ignitions.data 
            model_out.potential_suppression.data += model.potential_suppression.data 
            
            model_out.sensitivity_fuel.data += model.sensitivity_fuel.data 
            model_out.sensitivity_moisture.data += model.sensitivity_moisture.data 
            model_out.sensitivity_ignitions.data += model.sensitivity_ignitions.data 
            model_out.sensitivity_suppression.data += model.sensitivity_suppression.data 
            PDF_out += PDF
        #browser()
    nyrs = nmns / 12
    model_out.burnt_area_mean.data /= nyrs
    
    model_out.fuel.data /= nmns
    model_out.moisture.data /= nmns
    model_out.ignitions.data /= nyrs
    model_out.suppression.data /= nmns 
            
    model_out.standard_fuel.data /= nmns
    model_out.standard_moisture.data /= nmns
    model_out.standard_ignitions.data /= nmns
    model_out.standard_suppression.data /= nmns
            
    model_out.potential_fuel.data /= nyrs
    model_out.potential_moisture.data /= nyrs
    model_out.potential_ignitions.data /= nyrs
    model_out.potential_suppression.data /= nyrs

    model_out.sensitivity_fuel.data /= nyrs
    model_out.sensitivity_moisture.data /= nyrs
    model_out.sensitivity_ignitions.data /= nyrs
    model_out.sensitivity_suppression.data  /= nyrs
    
    PDF_out /= nmns
    model_out.burnt_area_pdf.data[:] = 0.0
    
    def posterize(i): return np.argmin(np.abs(i-BAs))
    vpost = np.vectorize(posterize)
    PDF_out0 = PDF_out.copy()
    PDF_out = vpost(PDF_out)
   
    dummy = PDF_out[0].copy()
    
    for PDF in PDF_out:
        for i in np.unique(PDF):
            dummy[:] = 0.0
            dummy[PDF == i] = 1.0
            try:
                model_out.burnt_area_pdf.data[i][mask] += dummy 
            except:
                browser()  
    
    return model_out

output_dir = '../ConFIRE_ISIMIP/outputs/sampled_posterior_ConFire_ISIMIP_solutions/'
mkDir(output_dir)

#models = [[runModelExperiment(input) for input in i] for i in inputs]
#model_example_file = output_dir + 'example_Confire.nc'
model = runModelExperiment(inputs[0][0], nmns = 2)
#browser()
#iris.save(model, model_example_file)
     

##########################
## perform boottrapping ##
##########################
output_controls = True
output_standard_limitation = True
output_potential_limitation = True
output_sensitivity = True
output_fullPost = True

n_posterior_sample = 50
qs = np.arange(1, 100, 1)

BAs = model.burnt_area_pdf.coord("model_level_number").points
BAs = 1/(1+np.exp(-BAs))


def weighted_percentile(data, percents, weights=None):
    ''' percents in units of 1%
        weights specifies the frequency (count) of data.
    '''
    if weights is None:
        return np.percentile(data, percents)
    ind=np.argsort(data)
    d=data[ind]
    w=weights[ind]
    p=1.*w.cumsum()/w.sum()*100
    y=np.interp(percents, p, d)
    return y

output_dir = output_dir + title_output +'/'
mkDir(output_dir)

def bootModel(input, model, experiment):
    output_diri = output_dir + '/' + model
    mkDir(output_diri)
    output_diri = output_diri + '/' + experiment + '/'
    
    outFileSumm = output_diri + 'model_summary-' + str(n_posterior_sample) + '.nc'
    outFilePost = output_diri + 'fullPost-' + str(n_posterior_sample) + '.nc'

    if os.path.isfile(outFilePost): return()
    mkDir(output_diri)
    fire_outi = []
    n_posterior = input[1].shape[0]
    ngap = int(n_posterior/n_posterior_sample)
    if ngap == 0: ngap = 1
    
    for i in range(0, n_posterior, ngap):
        outFile = output_diri + 'sample_no_' + str(i) +'.nc'
        if os.path.isfile(outFile):
            cubes = iris.load(outFile)
        else:          
            model = runModelExperiment(input, i, nyr = 1)        
            
            #model = ConFIRE(input_data, params.iloc[i])
            burnt_area = model.burnt_area_mean
            
            cubes = [model.burnt_area_mean]
            if output_fullPost:
                cubes = cubes + [model.burnt_area_pdf]
            
            if output_controls:
                cubes = cubes + [model.fuel, model.moisture, 
                                 model.ignitions, model.suppression]
            
            if output_standard_limitation:
                cubes = cubes + [model.standard_fuel, model.standard_moisture, 
                                 model.standard_ignitions, model.standard_suppression]
            
            if output_potential_limitation:
                cubes = cubes + [model.potential_fuel, model.potential_moisture,                
                                 model.potential_ignitions, model.potential_suppression]

            if output_sensitivity:
                cubes = cubes + [model.sensitivity_fuel, model.sensitivity_moisture,            
                                 model.sensitivity_ignitions, model.sensitivity_suppression]
            
            cubes = iris.cube.CubeList(cubes)             
            
            iris.save(cubes, outFile)
        nmes = [cube.name() for cube in cubes]   
        pdfID = np.where(np.array(nmes) == "burnt_area_pdf")[0][0]     
        MPDF = cubes[pdfID].copy()
        del cubes[pdfID]
        nmes = [cube.name() for cube in cubes]   
        cubes = iris.cube.CubeList([x for y, x in sorted(zip(nmes, cubes))])
        
        if output_fullPost:
            if i == 0:
                fullPost = MPDF.copy()
            else:
                fullPost.data += MPDF.data
        
        fire_outi = fire_outi + [cubes]
        
        print(outFile)
    
    fire_out = []
    for i in range(len(fire_outi[0])):
        #print(i)
        outi = []
        for out in fire_outi:
            outi = outi + [out[i].data]

        percentile_cube = newCubes3D("burnt_area", qs, out[i]) 
        percentile_cube.data = np.percentile(np.array(outi), qs, 0)

        fire_out = fire_out + [percentile_cube]

    fire_out = iris.cube.CubeList(fire_out) 
    print(outFileSumm)
    iris.save(fire_out, outFileSumm)    
    
    print(outFilePost)
    iris.save(fullPost, outFilePost) 
        
for model, i in zip(model_names, inputs):
#model = model_names[3]
#i = inputs[3]
    for experiment, input in zip(experiments, i):
        print(model)
        bootModel(input, model, experiment)
