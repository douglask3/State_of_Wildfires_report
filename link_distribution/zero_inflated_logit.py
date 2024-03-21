import pytensor
import pytensor.tensor as tt
import math     

class zero_inflated_logit(object):
    def __init__(self):
        pass

    def obs_given_model(Y, fx, sigma, p0, p1):
        '''return tt.sw1itch(
            tt.lt(Y, -150),
            -p0,
            -(1.0 - p0) *(1.0/(sigma * 2.506))*tt.exp(-0.5 * ((Y-fx)/sigma)**2)
        )
        '''
        pz = 1.0 - (fx**p1) * (1.0 - p0)
        return tt.switch(
            tt.lt(Y, -30),
            tt.log(pz),
            tt.log(1-pz) - tt.log(sigma * tt.sqrt(2*math.pi)) - ((Y-fx)**2)/(2*sigma**2)
        )


    
    def model_given_obs(self, Y, X, sigma, p0, p1):
        browser()
        return prob
