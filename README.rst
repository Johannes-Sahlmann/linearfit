linearfit
======

**python class that implements a general least-squares fit of a linear model using numpy matrix inversion.**

Uncertainties in the dependent variables (but not in the independent
variables) can be taken into account. All inputs have to be numpy matrices.

Math is based on Press'  
Numerical Receipes p661 : Section 15.2 Fitting Data to a Straight Line  
Numerical Receipes p671 : Section 15.4 General Linear Least Squares  

Code is based on an early yorick implementation by Damien Segransan
(University of Geneva)

Python implementation and tools by Johannes Sahlmann 2009-2017 (University of Geneva, European Space Agency, STScI/AURA)

Please see the test in test_linearfit.py for an example usage.


Basic example usage (Fitting a straight line to data with uncertainties in y)
-------------

```
# independent variable
x = np.array([0.0, 0.9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5,7.4])

# dependent variable	
y = np.array([5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5])

# weights of dependent variable	
wy = np.array([1,1.8,4,8,20,20,70,70,100,500])

# prepare matrices
M = np.mat(y);
#       diagonal covariance matrix of dependent variables
S = np.mat(np.diag(wy));
# matrix of independent variables, here only ones
C = np.mat(np.vstack( (np.ones(len(x)) , x)))

# initialise object
res = linearfit.LinearFit(M,S,C);

# do the fit
res.fit()
res.display_results()
```


Documentation
-------------

All classes and methods/functions include basic documentation. 


Installation notes
------------------

This package was developed in a python 2.7 environment, but was also
successfully tested using python 3.5.

The following python packages are required:

* numpy


How to run the example script
-----------

You may use pip for installation:  
```pip install linearfit```

Or get the source files, e.g.:   
```git clone https://github.com/johannes-sahlmann/linearfit```

Install pygacs:  
```cd linearfit```  
```python setup.py install --user```

To run the test, do:  
```python test_linearfit.py```


License
-------

Copyright (c) 2017 Johannes Sahlmann, STScI/AURA

linearfit is open source and free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
