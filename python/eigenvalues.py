#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2011-11-08 18:05:24.776750
#----------------------------------------------------------------------------#

import numpy as np
import scipy.linalg
import pylab

filename = 'matrix.npz'
data = np.load(filename)

condition = np.linalg.cond(data['matrix'])
print("The condition number is "+str(condition))
eigenvalues = scipy.linalg.eigvals(data['matrix'])
print("The eigenvalue with the largest real part is "+str(eigenvalues.max()))
print("The eigenvalue with the smallest real part is "+str(eigenvalues.min()))

pylab.plot(eigenvalues.real,eigenvalues.imag,'d')
pylab.show()
