class ConFIRE(object):
    def __init__(self, data, params, inference = False):
        """
        Initalise parameters and calculates the key variables needed to calculate burnt area.
        """
        self.inference = inference
        if self.inference:
            self.numPCK =  __import__('numpy')
        else:
            self.numPCK =  __import__('theano').tensor

        self.params = params
        

        ## finds controls
        self.fuel = self.control_fuel(data['cveg'], data['csoil'], self.params['c_csoil'])
        
        self.emcw = self.emc_weighted(data['humid'], data["precip"], self.params['wd_pg'])
        
        self.moisture = self.control_moisture(data['soilM'],
                                              self.emcw, data['trees'],
                                              self.params['c_emc2'], self.params['c_trees'], 
                                              self.params['kM'], self.params['pT'])

        self.ignitions = self.control_ignitions(data['pas'])

        self.suppression = self.control_suppression(data['crop'])

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
        if (not inference):
            try:
                self.burnt_area.long_name = "burnt_area"
                self.burnt_area_mode.long_name = "burnt_area_mode"
                #self.burnt_area_median.long_name = "burnt_area_median"
                self.burnt_area_mean.long_name = "burnt_area_mean"

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
            wet_days = 1 - self.numPCK.exp(-wd_pg * precip)
            emcw = wet_days + (1-wet_days) * emc
        except:
            emcw = emc.copy()
            emcw.data  = 1 - self.numPCK.exp(-wd_pg * precip.data)
            emcw.data = emcw.data + (1-emcw.data) * emc.data
        return(emcw)

    def control_moisture(self, shallow_soil, deeps_soil, emc, treeCover, cMs, cM, cMT, kM, pT):
        """
        Definition to describe moisture
        """
        moist = (shallow_soil/100.0 + cMs * deeps_soil/100.0 + cM*emc + cMT * (treeCover**pT)) / (1 + cMs + cM + cMT)
        moist.data = 1 - self.numPCK.log(1 - moist.data*kM)
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
        if self.inference:
            return sens
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
            out.data = 1.0/(1.0 + self.numPCK.exp(out.data))
            x = out
        except:
            x = -k * (x - x0)
            x = 1.0/(1.0 + self.numPCK.exp(x))
        return x
    
    def burnt_area_calPDF(self, data, p0, pp):
        
        mask = self.numPCK.logical_not(self.burnt_area_mode.data.mask)
        self.burnt_area_pdf = newCubes3D('burnt_area', 0.5, data['fireObs'])
        
        self.burnt_area_mean = self.burnt_area_mode.copy()
        self.burnt_area_mean.data[mask] = 0.0
        
        level_no = self.burnt_area_pdf.coord('model_level_number').points
        
        dist = norm(npLogit(self.burnt_area_mode.data[mask]), self.error)
        
        self.pz = 1.0 - (self.burnt_area_mode.data[mask]**pp) * (1.0 - p0)
        
        self.burnt_area_pdf.data[0][mask] = self.pz 
        x = self.burnt_area_pdf.coord('model_level_number').points
        for k in range(1, self.burnt_area_pdf.shape[0]):       
            self.burnt_area_pdf.data[k][mask] = dist.pdf(x[k]) * (1.0 - self.pz)
            self.burnt_area_mean.data[mask]  = self.burnt_area_mean.data[mask] +  dist.pdf(x[k]) * (1.0 - self.pz)  *(1/(1+self.numPCK.exp(-x[k])))
                                     
        PDFtot = self.burnt_area_pdf.collapsed(['model_level_number'], iris.analysis.SUM)
        
        self.burnt_area_mean.data = self.burnt_area_mean.data/ PDFtot.data
       
        
        
