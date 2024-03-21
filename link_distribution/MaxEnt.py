import pytensor
import pytensor.tensor as tt

class MaxEnt(object):

    def obs_given_model(Y, fx, qSpread = None, CA = None):
        """calculates the log-transformed continuous logit likelihood for Y given fx when Y
            and fx are probabilities between 0-1 with relative areas, CA
            Works with tensor variables.   
        Arguments:
            Y  -- Y  in P(Y|fx). numpy 1-d array
	    fx -- fx in P(Y|fx). tensor 1-d array, length of Y
            qSpread -- parameter that inflates target Y, to account for potential 
                    multi-fire overlap
            CA -- Area for the cover type (cover area). numpy 1-d array, length of Y. 
            Default of None means everything is considered equal area.
        Returns:
            1-d tensor array of liklihoods.
        """
    
        fx = tt.switch( tt.lt(fx, 0.0000000001), 0.0000000001, fx)
        fx = tt.switch( tt.gt(fx, 0.9999999999), 0.9999999999, fx)
        
        if qSpread is not None:
            Y = Y *(1 + qSpread) / (Y * qSpread + 1)
      
        if CA is not None: 
            prob =  Y*CA*tt.log(fx) + (1.0-Y)*CA*tt.log((1-fx))
        else:
            prob = Y*tt.log(fx) + (1.0-Y)*tt.log((1-fx))
        return prob
    
    def model_given_obs(Y, X, *args, **kw):
        X1 = 1 - X
        def prob_fun(y):
            return (y**X) * ((1-y)**X1)
        prob = prob_fun(Y)/prob_fun(X)
        return prob
    

