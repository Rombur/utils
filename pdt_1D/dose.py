# Python code
# Author: Bruno Turcksin
# Date: 2010-11-20 14:56:31.885611

#----------------------------------------------------------------------------#
## Class dose                                                               ##
#----------------------------------------------------------------------------#

import numpy as np
import gzip
import pickle
import os.path
import utils

"""This module computes the dose for pdt."""

class dose :
    """This class computes the dose through the flux given an pdt output file. 
    The beam has to be along the z direction"""

    def __init__(self, filename, xs_filename, n_processors, n_division_z, 
            delta_z, n_g_groups, n_e_groups, n_p_groups,x_position, y_position, 
            n_moments,moments,groups) :
        self.filenames = filename+".output_"
        self.xs_filename = xs_filename
        self.n_processors = n_processors
        self.n_division_z = n_division_z
        self.delta_z = delta_z
        self.n_g_groups = n_g_groups
        self.n_e_groups = n_e_groups
        self.n_p_groups = n_p_groups
        self.n_groups = n_g_groups+n_e_groups+n_p_groups
        self.x_position = x_position
        self.y_position = y_position
        self.n_moments = n_moments-1
        self.moments = moments
        self.groups = groups
        self.n_cells = len(x_position)*len(y_position)
        self.flux = np.zeros((2*n_division_z,self.n_groups))
        self.moment_flux = np.zeros((2*n_division_z,len(moments),len(groups)))
        self.utils = utils.utils()

#----------------------------------------------------------------------------#

    def read_moments(self,line,moment,group,element,pos) :
        """This function read the pdt input files and computes some chosen
        moments."""

        value,read = self.utils.read_float(line,3)
        if read == False :
            self.utils.abort("Problem while reading the moments of the flux.")
        else :
            mom = self.moments.index(moment)
            grp = self.groups.index(group)
            if element < 4 :
                self.moment_flux[pos,mom,grp] += value[2]/(4.0*self.n_cells)
            else :
                self.moment_flux[pos+1,mom,grp] += value[2]/(4.0*self.n_cells)

#----------------------------------------------------------------------------#

    def read_flux(self) :
        """This function reads the pdt input files and compute the flux the 1D
        flux. For a given Z, the flux is averaged by group.""" 

        if not os.path.exists("flux.txt.gz") :
            for i in xrange(0,self.n_processors) :
                print "reading file %i"%i
                file_obj = open(self.filenames+str(i),'r')
                eof = False
                counter = 0
                while eof == False :
                    line = file_obj.readline()
                    counter += 1
                    if counter%1000000 == 0 :
                        print "Reading line %i"%counter
                    eof = self.utils.search_in_line(line,"TIMING")
                    if counter > 100000000 :
                        eof = True
                    cell_id = self.utils.search_in_line(line,"cell id,")
                    if cell_id :
                        line = file_obj.readline()
                        counter += 1
                        value,read = self.utils.read_float(line,4)
                        if read == False :
                            self.utils.abort(
                                    "Problem while reading the position of the cell.")
                        if value[1] in self.x_position and\
                                value[2] in self.y_position :
                            flux_pos = 2.0*np.floor(value[3]/self.delta_z)
                            element = 0
                            while element < 8 :
                                line = file_obj.readline()
                                counter += 1
                                partial_flux = self.utils.search_in_line(line,
                                        "flux for element")
                                if partial_flux == True :
                                    for j in xrange(0,self.n_groups) :
                                        line = file_obj.readline()
                                        counter += 1
                                        value,read = self.utils.read_float(line,3)
                                        if read == False :
                                            self.utils.abort(
                                                    "Problem while reading the flux")
                                        if element < 4 :
                                            self.flux[flux_pos,j] += value[2]/(4.0*
                                                    self.n_cells)
                                        else :
                                            self.flux[flux_pos+1,j] += value[2]/(4.0*
                                                    self.n_cells)
                                        if 0 in self.moments and j in self.groups :
                                            self.read_moments(line,0,j,element,
                                                    flux_pos)
                                        for k in xrange(0,self.n_moments) :
                                            line = file_obj.readline();
                                            counter += 1
                                            if k+1 in self.moments and\
                                                    j in self.groups :
                                                        self.read_moments(line,
                                                                k+1,j,element,
                                                                flux_pos)
                                    element += 1
                file_obj.close()
        else :
            self.flux = np.loadtxt("flux.txt.gz")
            self.moment_flux = np.loadtxt("moments.txt.gz")
            self.moment_flux = self.moment_flux.reshape((2*self.n_division_z,
                len(self.moments),len(self.groups)))

#----------------------------------------------------------------------------#

    def print_flux(self) :
        """Print the value of the flux to be read faster if we need to read
        them again."""
        
        np.savetxt("flux.txt.gz",self.flux)

        outfile = gzip.open("moments.txt.gz",'w')
# any line starting with # will be ignored by numpy.loadtxt
# The 0 is because we put there the 0th element of format
        outfile.write('# array shape : {0}\n'.format(self.moment_flux.shape))
        for i in xrange(0,len(self.moments)) :
            np.savetxt(outfile,self.moment_flux[:,i,:])
# writing out a break to indicate a new moment
            outfile.write("# New moment\n")
        outfile.close()

        paramfile = open("parameters.txt","w")
        out =[self.n_division_z,self.delta_z,len(self.groups),len(self.moments), 
                self.groups,self.moments]
        pickle.dump(out,paramfile)
        paramfile.close()

#----------------------------------------------------------------------------#

    def read_xs(self) :
        """This function reads the energy deposition cross section."""

        file_obj = open(self.xs_filename,'r')

        self.g_energy_xs = np.zeros(self.n_g_groups)
        self.e_energy_xs = np.zeros(self.n_e_groups)
        self.p_energy_xs = np.zeros(self.n_p_groups)

        eof = False
        while eof == False :
            line = file_obj.readline()
            eof = self.utils.search_in_line(line,"MT 2501")
            g_xs = self.utils.search_in_line(line,"MT 542")
            e_xs = self.utils.search_in_line(line,"MT 603")
            p_xs = self.utils.search_in_line(line,"MT 703")

            if g_xs :
               n_lines = int(np.floor(self.n_g_groups/5))
               for i in xrange(0,n_lines) :
                   line = file_obj.readline()
                   values,read = self.utils.read_float(line,5)
                   if read == False :
                       self.utils.abort(
                               "Problem while reading the energy xs for the gammas")
                   for j in xrange(0,5) :
                       self.g_energy_xs[i*5+j] = values[j]
               line = file_obj.readline()
               to_read = self.n_g_groups-5*n_lines
               values,read = self.utils.read_float(line,to_read)
               if read == False :
                   self.utils.abort(
                           "Problem while reading the energy xs for the gammas")
               for i in xrange(0,to_read) :
                   self.g_energy_xs[n_lines*5+i] = values[i]

            if e_xs :
                n_lines = int(np.floor(self.n_e_groups/5))
                for i in xrange(0,n_lines) :
                    line = file_obj.readline()
                    values,read = self.utils.read_float(line,5)
                    if read == False :
                        self.utils.abort(
                                "Problem while reading the energy xs for the electrons")
                    for j in xrange(0,5) :
                        self.e_energy_xs[i*5+j] = values[j]
                line = file_obj.readline()
                to_read = self.n_e_groups-5*n_lines
                values,read = self.utils.read_float(line,to_read)
                if read == False :
                    self.utils.abort(
                            "Problem while reading the energy xs for the electrons")
                for i in xrange(0,to_read) :
                    self.e_energy_xs[n_lines*5+i] = values[i]

            if p_xs :
                n_lines = int(np.floor(self.n_p_groups/5))
                for i in xrange(0,n_lines) :
                    line = file_obj.readline()
                    values,read = self.utils.read_float(line,5)
                    if read == False :
                        self.utils.abort(
                                "Problem while reading the energy xs for the positrons")
                    for j in xrange(0,5) :
                        self.p_energy_xs[i*5+j] = values[j]
                line = file_obj.readline()
                to_read = self.n_p_groups-5*n_lines
                values,read = self.utils.read_float(line,to_read)
                if read == False :
                    self.utils.abort(
                            "Problem while reading the energy xs for the positrons")
                for i in xrange(0,to_read) :
                    self.p_energy_xs[n_lines*5+i] = values[i]

#----------------------------------------------------------------------------#

    def compute_dose(self,position) :
        """Compute the dose at a given position."""

        value = 0.0
        
        pos = 2*int(np.floor(position/self.delta_z))
        for i in xrange(0,self.n_groups) :
            flux = self.flux[pos,i] + (position-pos*self.delta_z/2.0)/self.delta_z *\
                    (self.flux[pos+1,i] - self.flux[pos,i])
            if i < self.n_g_groups :
                value += self.g_energy_xs[i]*flux
            elif i < (self.n_g_groups+self.n_e_groups) :
                value += self.e_energy_xs[i-self.n_g_groups]*flux
            else :
                value += self.p_energy_xs[i-self.n_g_groups-self.n_e_groups]*\
                        flux
        print value
