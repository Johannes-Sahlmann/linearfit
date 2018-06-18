from __future__ import print_function

import numpy as np

class LinearFit(object):
    """
    Structure Class that performs a general least-squares fit of a linear model using numpy matrix inversion.
    Uncertainties in the dependent variables (but not in the independent variables) can be taken into account.
    All inputs have to be numpy matrices.
    
    Math is based on Press' 
    Numerical Receipes p661 : Section 15.2 Fitting Data to a Straight Line
    Numerical Receipes p671 : Section 15.4 General Linear Least Squares

    Code is based on an early yorick implementation by Damien Segransan (University of Geneva)
    Python implementation and tools by Johannes Sahlmann 2009-2017 (University of Geneva, European Space Agency, STScI/AURA)

    Attributes
    ----------
    dependent_variable : numpy matrix (1xN)
        dependent_variables of the linear equation system (N equations, M unknown coefficients)
    inverse_covariance_matrix : numpy matrix (NxN)
        inverse covariance matrix corresponding to the dependent_variable.
        i.e. data weights proportional to 1/sigma**2 where sigma=uncertainty
    independent_variable : numpy matrix (MxN)   
        the independent_variables that are multiplied by the unknown coefficients
    
    
    Calculated Attributes
    ----------
    p: numpy array
        coefficients of the solution
    p_formal_uncertainty: numpy array           
        formal uncertainty  of the coefficients
    p_formal_covariance_matrix: numpy matrix    
        formal covariance matrix of the coefficients (not rescaled)
    p_normalised_uncertainty: numpy array    
        normalised uncertainty (chi2 = 1) of the coefficients
    p_normalised_covariance_matrix: numpy matrix    
        normalised covariance matrix of the coefficients (rescaled to yield chi2=1)
    p_correlation_matrix: numpy matrix 
        coefficient correlation matrix
    fit: numpy array         
        values of the best-fit model        
    residuals:  numpy array    
        Observed - Calculated (O-C) residuals         
    chi2: float
        chi-square value of the best fit
    
    
    Methods
    ----------
    fit(): 
        perform the fit
    
    display_results(par_names=None,format=None,precision = 3, scaleFactor = 1., nformat='f'): 
        show results on stdout
    
    display_correlations(par_names=None, format=None): 
        show correlation matrix on stdout
        
    """    

    def __init__(self, dependent_variable, inverse_covariance_matrix, independent_variable):    
        ''' Initialize object '''    
        self.y_i  = dependent_variable
        self.X_ij = independent_variable
        self.inverse_covariance_matrix = inverse_covariance_matrix


    def fit(self):
        '''
        Perform the linear fit. This populates the object attributes.
        '''
        
        # number of measurements/constraints, i.e. number of equations
        self.n_measurement   = self.X_ij.shape[1];
        # number of coefficients, i.e. free parameters
        self.n_param = self.X_ij.shape[0];
        # number of degrees of freedom
        self.n_freedom  = self.n_measurement - self.n_param;
  
        
        A_ij = (self.X_ij) * self.inverse_covariance_matrix;
        alpha_kj = A_ij * (self.X_ij.T);
        beta_k = A_ij * self.y_i.T;
    
        a_j = np.linalg.solve(alpha_kj,beta_k);
    
        yfit_i = (self.X_ij.T * a_j).T;
        chi2 = ((yfit_i-self.y_i)*self.inverse_covariance_matrix)*(yfit_i-self.y_i).T;    
        C_jk = np.linalg.inv(alpha_kj);  
        var_aj = np.mat( np.diag(C_jk) );
        
        #       the variance on the fitted coefficients are the diagonal terms of the coefficient covariance matrix
        #       if the error bars are not well estimated, it is useful to renormalize the variance by the measured chi2
        #       divided by the expected chi2. That means the normalised variances take into account the residual dispersion in  the data
        normalized_variance_aj =  var_aj.T * chi2/self.n_freedom;
        stdev_aj = np.sqrt(normalized_variance_aj);
        
        # coefficients of the solution
        self.p = np.array(a_j).flatten()
        
        # normalised uncertainty (chi2 = 1) of the coefficients
        self.p_normalised_uncertainty = np.array(stdev_aj).flatten()
        
        # values of the best-fit model
        self.fit_values = np.array(yfit_i).flatten()
        
        # Observed - Calculated (O-C) residuals 
        self.residuals = np.array(self.y_i-yfit_i).flatten()
        
        # chi-square value of the best fit
        self.chi2 = np.array(chi2)[0][0]
        
        # formal uncertainty  of the coefficients
        self.p_formal_uncertainty = np.array(np.sqrt(var_aj)).flatten()
        
        # formal covariance matrix of the coefficients (not rescaled)
        self.p_formal_covariance_matrix = C_jk

        # normalised covariance matrix of the coefficients (rescaled to yield chi2=1)
        self.p_normalised_covariance_matrix = C_jk * self.chi2/self.n_freedom


        #     compute correlation Matrix
        tmp_v = (1./self.p_formal_uncertainty);
        tmp = np.vstack( ( tmp_v , np.tile( np.zeros(len(tmp_v)) ,  (len(tmp_v) - 1,1) )) );
        M = np.mat(tmp.T)  
        err_matrix = M.dot(M.T)
        correlation_matrix =  np.multiply(err_matrix,self.p_formal_covariance_matrix.T);
        self.p_correlation_matrix = correlation_matrix;
     
             
    def display_results(self,par_names=None,format=None,precision = 3, scale_factor = 1., nformat='f'):
        '''
        Display results on stdout. Optional arguments include parameter names, decimal precision and number format.        
        '''

        if format is None:
            pm_string = '+/-'
        elif format == 'latex': 
            pm_string = '\pm'
    
        if par_names is None:
            par_names = np.array(['Param%d' % (i+1) for i in range(len(self.p))])            
        for i in range(len(self.p)):
            if nformat=='f':
                print('%s\t = %.*f %s %.*f (normalised for chi2=1)  +/- %.*f (formal uncertainty)' % (par_names[i], precision,self.p[i]*scale_factor,pm_string,precision,self.p_normalised_uncertainty[i]*scale_factor, precision,self.p_formal_uncertainty[i]*scale_factor ));
            elif nformat=='e':
                print('%s\t = %.*e %s %.*e (normalised for chi2=1)  +/- %.*e (formal uncertainty)' % (par_names[i], precision,self.p[i]*scale_factor,pm_string,precision,self.p_normalised_uncertainty[i]*scale_factor, precision,self.p_formal_uncertainty[i]*scale_factor ));
        
        if nformat=='f':
            print('reduced chi2 = %3.3f\t' % (self.chi2/self.n_freedom), end="")
            print('O-C RMS        %3.3f'   % (np.std(self.residuals)*scale_factor))
        elif nformat=='e':
            print('reduced chi2 = %3.3e\t' % (self.chi2/self.n_freedom), end="")
            print('O-C RMS        %3.3e'   % (np.std(self.residuals)*scale_factor))
                
                
    def display_correlations(self,par_names=None, format=None):        
        '''
        Display correlation matrix on stdout.
        '''
        
        matrix = self.p_correlation_matrix;        
        # display it
        display_square_matrix(matrix,par_names=par_names, format=format)

    

def display_square_matrix(matrix,par_names=None, format=None):        
        '''
        Display a square matrix on stdout. Optional output in latex format.
        '''
    
        print('*'*10 +' Correlation Matrix '+'*'*10);
        Nmatrix = matrix.shape[0];
        Mmatrix = matrix.shape[1];        
        if par_names is None:
            par_names = np.array(['P%d' % (i+1) for i in range(Nmatrix)])  
        
        if format is None:
            print('\t\t', end="")
            for i in range(Nmatrix ):   
                print('%s\t' % par_names[i], end="")
            print('');   
            for i in range(Nmatrix ):
                print('%s \t\t' % par_names[i], end="")
                for j in range(Mmatrix ):
                    if i>=j:
                        print('%+1.2f \t' % matrix[i,j], end="")
                print('');   
                
        elif format == 'latex':             
            print('\\begin{table} \n \\caption{Correlation matrix for }');
            print('\\label{tab:XXXcorr}  \\centering ');
            print('\\begin{tabular}{l |' + 'r@{\\;\\;} '* (Nmatrix+1) + '} \n \\hline\\hline'); 
            print('& '),
            for i in range(Nmatrix-1):
                print('%s &' % par_names[i]), 
            print('%s \\\\ \\hline' % par_names[-1]);
            for i in range(Nmatrix):
                print('%s &' % par_names[i], end="")
                for j in range(Mmatrix):
                    if j == Mmatrix & i == Nmatrix:
                        print('$%+1.2f$  ' % matrix[i,j], end="")
                    elif i>=j:
                        print('$%+1.2f$ &' % matrix[i,j], end="")
                print(' \\\\');
            print('\\hline \n \\end{tabular} \n\\end{table} \n');        
    
        print('*'*50)



         

