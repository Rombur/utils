#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2010-12-01 17:29:36.177796
#----------------------------------------------------------------------------#

import pickle
import numpy as np
import pylab

def flux(position,delta_z,n_division_z,group,moment,moment_flux) :
    pos = 2*int(np.floor(position/delta_z))
    offset = 2*moment*n_division_z
    value = moment_flux[pos+offset,group]+(position-pos*delta_z/2.0)/delta_z*\
            (moment_flux[pos+1+offset,group] -moment_flux[pos+offset,group])
    return value

parameters = pickle.load(open("parameters.txt"))
n_division_z = parameters[0]
delta_z = parameters[1]
n_groups = parameters[2]
n_moments = parameters[3]
groups = parameters[4]
moments = parameters[5]

moment_flux = np.loadtxt("moments.txt.gz")
x_max = n_division_z * delta_z
x = np.arange(0.0,x_max,x_max/100.0)
y = np.zeros_like(x)
for g in xrange(0,n_groups) :
    for m in xrange(0,n_moments) :
        for i in xrange(0,100) :    
            y[i]=flux(x[i],delta_z,n_division_z,g,m,moment_flux)
        pylab.figure()
        pylab.plot(x,y)
        pylab.savefig("group_{0}_moment_{1}.eps".format(groups[g],moments[m]))
        pylab.savefig("group_{0}_moment_{1}".format(groups[g],moments[m]))
        # pylab.show()
