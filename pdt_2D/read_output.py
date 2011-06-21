# Python code
# Author: Bruno Turcksin
# Date: 2011-06-20 12:35:05.775594

#----------------------------------------------------------------------------#
## Class read_output                                                        ##
#----------------------------------------------------------------------------#

import cPickle as pickle
import gzip
import utils

"""This module reads the pdt output"""

class read_output(object) :
  """This class reads the pdt output and creates a file with the corner of
  each cell (2D) and the flux at these points"""

  def __init__(self,filename_input,filename_output,n_proc=1) :
    
   super(read_output,self).__init__()
   self.file_input = filename_input
   self.file_output = filename_output
   self.n_processors = n_proc;
   self.x = []
   self.y = []
   self.values = []

#----------------------------------------------------------------------------#

  def read_file(self) :
    """This function reads the input files and store the data."""

    for i in xrange(0,self.n_processors) :
      print "reading file %i"%i
      file_obj = open(self.file_input+str(i),'r')
      eof = False
      counter = 0
      while eof == False :
        line = file_obj.readline()
        counter += 1
        if counter%1000000 == 0 :
          print "Reading line %i"%counter
        eof = utils.search_in_line(line,"TIMING")
        if counter > 100000000 :
          eof = True
        cell_id = utils.search_in_line(line,"element id,")
        if cell_id :
          for j in xrange(0,4) :
            line = file_obj.readline()
            counter += 1
            value,read = utils.read_float(line,4)
            if read == False :
                utils.abort("Problem while reading the corners of the cell.")
            self.x.append(value[1])
            self.y.append(value[2])

# skipe some lines
          for j in xrange(0,2) :
            line = file_obj.readline()
            counter += 1

          for j in xrange(0,4) :
            line = file_obj.readline()
            counter += 1
            line = file_obj.readline()
            counter += 1 
            value,read = utils.read_float(line,3)
            if read == False :
                utils.abort(
                "Problem while reading the value at the corner of the cell.")
            self.values.append(value[2]) 

      file_obj.close()

#----------------------------------------------------------------------------#

  def print_flux(self) :
    """Print the values of the flux in a file."""

    outfile = gzip.open(self.file_output,'w')
    out = [self.x,self.y,self.values]
    pickle.dump(out,outfile)
    outfile.close()
