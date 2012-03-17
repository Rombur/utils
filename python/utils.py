# Python code
# Author: Bruno Turcksin
# Date: 2010-10-08 15:51:19.904088

#----------------------------------------------------------------------------#
## Module utils                                                              ##
#----------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot

"""Module containing some miscellenous function."""

def Abort(message) :
  print (message)
  print ("The program aborted.")
  exit()

#----------------------------------------------------------------------------#

def Search_in_line(line,keyword) :
  """Search for the given keyword in a given line. If the keyword is
  found, return True else return False."""

  found = False;
  position = line.find(keyword)
  if position != -1 :
      found = True

  return found

#----------------------------------------------------------------------------#

def Read_float(line,n_floats) :
  """Read the first n_floats in the given line. The floats have to be 
  separate by blanks. The values are returned in a list."""

  vector = []
  pos = 0
  read = True
  for i in xrange(0,n_floats) :
# We read the next element
      value = ' '
# We skip the blanks
      while line[pos] == " " :
          pos += 1
# We read the value
      while line[pos] != " " :
          if line[pos] != "\n" :
              value += line[pos]
              pos += 1
          else :
              break
      if value != " " :
          value = float(value)
      else :
          read = False
          break
      vector.append(value)

  return vector,read

#----------------------------------------------------------------------------#

def Plot_pattern(matrix) :
  """Plot the pattern of the given matrix, i.e. plot the non-zeros elements of
  the matrix."""

  x = []
  y = []
  for j in xrange(0,matrix.shape[1]) :
    for i in xrange(0,matrix.shape[0]) :
      if np.abs(matrix[i,j]) > 1e-12 :
        x.append(matrix.shape[0]-i)
        y.append(j+1)
  matplotlib.pyplot.plot(y,x,'x')
  matplotlib.pyplot.xlim(0,matrix.shape[1]+1)
  matplotlib.pyplot.ylim(0,matrix.shape[0]+1)
  matplotlib.pyplot.show()
