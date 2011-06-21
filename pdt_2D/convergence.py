# Python code
# Author: Bruno Turcksin
# Date: 2011-06-20 14:19:45.851845

#----------------------------------------------------------------------------#
## Class convergence                                                        ##
#----------------------------------------------------------------------------#

import cPickle as pickle
import gzip
import numpy as np
import pylab
import utils

"""Compute the convergence order of a method"""

class convergence(object) :
  """This class computes the order of convergence of a method."""

  def __init__(self,list_filename,reference) :
    super(convergence,self).__init__()
    self.list_filename = list_filename
    self.reference = reference

#----------------------------------------------------------------------------#

  def compute_order(self) :
    """Compute the order of convergence."""

    error = []
    dofs = []
# read the refence solution
    x_ref,y_ref,val_ref = self.read_flux(self.reference)

# loop over the files
    for filename in self.list_filename :
# read the solution
      x,y,val = self.read_flux(filename)
      tmp_error = 0.
      dofs.append(len(x))

# compute the double integral
      cpt = 0  
      for i in xrange(0,len(x)/4) :
        offset = 4*i
        for j in xrange(0,len(x_ref)/4) :
          offset_ref = 4*j
          if x_ref[offset_ref]>=x[offset] and\
              x_ref[offset_ref+1]<=x[offset+1] and\
              y_ref[offset_ref]>=y[offset] and\
              y_ref[offset_ref+2]<=y[offset+2] :
                tmp_error += 0.25*(val_ref[offset_ref]-self.compute_values(\
                    x,y,val,offset,x_ref[offset_ref],y_ref[offset_ref]))**2+\
                    0.25*(val_ref[offset_ref+1]-self.compute_values(x,y,val,\
                    offset,x_ref[offset_ref+1],y_ref[offset_ref+1]))**2+\
                    0.25*(val_ref[offset_ref+2]-self.compute_values(x,y,val,\
                    offset,x_ref[offset_ref+2],y_ref[offset_ref+2]))**2+\
                    0.25*(val_ref[offset_ref+3]-self.compute_values(x,y,val,\
                    offset,x_ref[offset_ref+3],y_ref[offset_ref+3]))**2
                cpt +=1 
       
      error.append(tmp_error)

    np_error = np.array(error) 
    n_dofs = np.array(dofs)
        
    print n_dofs
    print np_error
    a,b = np.polyfit(np.log10(n_dofs),np.log10(np_error),1)
    print 'parameters ax+b : a='+str(a)+' and b='+str(b)

    pylab.plot(np.log10(n_dofs),np.log10(np_error))
    pylab.show()

#----------------------------------------------------------------------------#

  def read_flux(self,filename) :
    """Read the values saves and put them in x, y and val."""

    input_file = gzip.open(filename,'r')
    data = pickle.load(input_file)
    x = data[0]
    y = data[1]
    val = data[2]

    input_file.close()

    return x,y,val

#----------------------------------------------------------------------------#

  def compute_values(self,x,y,val,offset,x_ref,y_ref) :
    """Compute the value of the flux at the point (x_ref,y_ref)."""

    a = x[offset+1]-x[offset]
    b = y[offset+2]-y[offset]
    pos_x = x_ref-x[offset]
    pos_y = y_ref-y[offset]
    
    flux = val[offset]*((a-pos_x)/a)*((b-pos_y)/b)+val[offset+1]*\
        (pos_x/a)*((b-pos_y)/b)+val[offset+2]*(pos_x/a)*(pos_y/b)+val[offset+3]*\
        ((a-pos_x)/a)*(pos_y/b)

    return flux
