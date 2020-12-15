from __future__ import print_function

import unittest
import numpy as np
from linearfit import linearfit
# import linearfit

class LinearFitTestCase(unittest.TestCase):
    def setUp(self):
        #   Pearson's data retrieved from http://www.astro.rug.nl/software/kapteyn/kmpfittutorial.html#fitting-data-when-both-variables-have-uncertainties
        x = np.array([0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4])
        y = np.array([5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5])
        wy = np.array([1,1.8,4,8,20,20,70,70,100,500])
        erry = 1/np.sqrt(wy)
        N = len(x)

        #       weights in x are ignored in linearfit       
        #       wx = np.array([1000.0,1000,500,800,200,80,60,20,1.8,1.0])
        #       errx = 1/np.sqrt(wx)  

        # dependent variables
        M = np.mat(y);

        #       diagonal covariance matrix of dependent variables
        S = np.mat(np.diag(wy));        

        # matrix of independent variables, here only ones
        C = np.mat(np.vstack( (np.ones(len(x)) , x)))    
        
        # initialise object
        res = linearfit.LinearFit(M,S,C);
        
        # do the fit 
        res.fit()        
        
        print("\n\n======== Results linearfit =========")            
        res.display_results(precision=10)
        res.display_correlations()  

        self.result = res

        
        #         optional comparison with kapteyn (nonlinear fitting) package
        if 0:                    
            from kapteyn import kmpfit
            def model(p, x):
               a, b = p
               return a + b*x

            def residuals2(p, data):
               # Residuals function for data with errors in y only
               a, b = p
               x, y, ey = data
               wi = np.where(ey==0.0, 0.0, 1.0/ey)
               d = wi*(y-model(p,x))
               return d
        
            # Prepare fit routine
            beta0 = [5.0, 1.0]         # Initial estimates
        
            # Prepare fit routine
            fitobj2 = kmpfit.Fitter(residuals=residuals2, data=(x, y, erry))
            fitobj2.fit(params0=beta0)
            print("\n\n======== Results kmpfit =========")
            print("Fitted parameters:      ", fitobj2.params)
            print("Covariance errors:      ", fitobj2.xerror)
            print("Standard errors         ", fitobj2.stderr)
            print("Reduced Chi^2:          ", fitobj2.rchi2_min)
        
        
        

    def test_linearfit_coefficients(self):
        self.assertTrue(np.max(self.result.p - np.array([6.1001092933631265,-0.6108129531048464])) < 1e-7, 'incorrect solution coefficients')

if __name__ == '__main__':
    unittest.main()