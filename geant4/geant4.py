# Python code
# Author: Bruno Turcksin
# Date: 2010-10-08 15:39:26.070620

#----------------------------------------------------------------------------#
## Class geant4                                                             ##
#----------------------------------------------------------------------------#

import utils

"""Module for Geant4"""

class geant4  :
    """Class containing functions to manipulate Geant4 inputs/outputs."""

    def __init__(self,input_file,output_file) :

# Open the input file
        self.file = open(input_file,'r')
# Store the name of the output file
        self.out = output_file

#----------------------------------------------------------------------------#

    def read_file(self) :
        """Read the file and store the position of the energy deposition.""" 

        self.x = []
        self.y = []
        self.z = []
        self.energy_deposition = []

        line = self.file.readline()
        utility = utils.utils()
        end = utility.search_in_line(line,"Quit")
        while end==False :
            step = utility.search_in_line(line,"Step#")
            if step :
                line = self.file.readline()
                Track = utility.search_in_line(line,"Track")
                while Track==False :
                    value,read = utility.read_float(line,6)
                    if read :
                        if value[5] != 0 :
                            self.x.append(value[1])
                            self.y.append(value[2])
                            self.z.append(value[3])
                            self.energy_deposition.append(value[5])
                    else :
                        utility.abort("Problem while reading the input file.")
                    line = self.file.readline()
                    Track = utility.search_in_line(line,"Track")

            else :
                line = self.file.readline()
            end = utility.search_in_line(line,"Quit")
        
        self.file.close()

#----------------------------------------------------------------------------#

    def write_output(self) :
        """Write the output file."""

# Open the output file 
        output_file = open(self.out,'w')
# Write in the file
        i_max = len(self.x)
        for i in xrange(0,i_max) :
            line = str(self.x[i])+"\t"+str(self.y[i])+"\t"+str(self.z[i])+\
                    "\t"+str(self.energy_deposition[i])+"\n"
            output_file.write(line)
        output_file.close()
