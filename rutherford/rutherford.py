# Python code
# Author: Bruno Turcksin
# Date: 2010-11-03 17:21:49.271082

#----------------------------------------------------------------------------#
## Class rutherford                                                         ##
#----------------------------------------------------------------------------#

import math
import numpy

"""Module which compute the data to run the 9 tests of Morel."""

class rutherford  :
    """Compute the cross section and the width for the tests"""

    def __init__(self,test, tc) :
        """Constructor needs the number of the test that we want to use [1,9]
        and if we want to use the transport correction."""

        self.tc = tc;
        self.L_max = 11
        self.Z = 13.
        self.A = 26.98
        self.test = test
        if test == 1 :
            self.xs_tot = 6.5/5.
            self.E = 0.01/0.511
        elif test == 2 :
            self.xs_tot = 51.5/50.
            self.E = 0.1/0.511
        elif test == 3 :
            self.xs_tot = 501.5/50.
            self.E = 1/0.511
        elif test == 4 :
            self.xs_tot = 13./10.
            self.E = 0.01/0.511
        elif test == 5 :
            self.xs_tot = 103./100.
            self.E = 0.1/0.511
        elif test == 6 :
            self.xs_tot = 1003./1000.
            self.E = 1/0.511
        elif test == 7 :
            self.xs_tot = 16./10.
            self.E = 0.01/0.511
        elif test == 8 :
            self.xs_tot = 106./100.
            self.E = 0.1/0.511
        elif test == 9 :
            self.xs_tot = 1006./1000.
            self.E = 1/0.511

#----------------------------------------------------------------------------#

    def compute_xs(self) :
        """Compute the value of the cross sections."""

# define some parameters
        N_avo = 6.022214179E23
        r_0 = 2.8179402894E-13
        eta = 0.25*(self.Z**(1./3.)/(0.885*137.))**2/self.E/(self.E+2)*\
                (1.13+3.76*(self.Z/137.)**2*(self.E+1)**2/self.E/(self.E+2))
        cst = 2*math.pi*(self.Z+1)*N_avo*r_0**2*self.Z/self.A*(self.E+
                1)**2/(self.E**2*(self.E+2)**2)
        sca_elec = cst/(2*eta*(eta+1))

        print self.xs_tot
# compute the total cross section
        self.xs_tot = self.xs_tot * sca_elec;

        print self.xs_tot
# compute the moment of the scattering cross section
        coef = numpy.zeros(self.L_max+2)
        self.xs_sca = numpy.zeros(self.L_max+1)
        for i in xrange(1,self.L_max+2) :
            if i == 1 :
                coef[i] = math.log(1+1./eta)-1./(1+eta)
            else :
                j = float(i-1)
                coef[i] = (2+1./j)*(1+2.*eta)*coef[i-1]-(1+1./j)*coef[i-2]-\
                        (2+1./j)/(1+eta);
        for i in xrange(0,self.L_max+1) :
            self.xs_sca[i] = sca_elec - cst*coef[i]

# transport correction
        if self.tc :
            correction = sca_elec - cst*coef[-1]
            self.xs_tot = self.xs_tot - correction
            for i in xrange(0,self.L_max+1) :
                self.xs_sca[i] = self.xs_sca[i] - correction
        print self.xs_tot

#----------------------------------------------------------------------------#

    def compute_width(self) :
        """Compute width of the medium."""

        if self.test == 1 :
            self.width = 5*7.16127187911E-6
        elif self.test == 2 :
            self.width = 50*3.466E-5
        elif self.test == 3 :
            self.width = 500*9.5988E-5
        elif self.test == 4 :
            self.width = 7.1613E-5
        elif self.test == 5 :
            self.width = 3.466E-3
        elif self.test == 6 :
            self.width = 9.5988E-2
        elif self.test == 7 :
            self.width = 7.1613E-5
        elif self.test == 8 :
            self.width = 3.466E-3
        elif self.test == 9 :
            self.width = 9.5988E-2

#----------------------------------------------------------------------------#

    def write(self,filename) :
        """Append in the given file the cross section and on the screen the
        width of the domain."""

        blank = "              "
        MT_2501 = "MT 2501, Moment "
        sink = "  Sink, first, last:     0   0   0\n"

        file = open(filename,'a')
        file.write("MT 1\n")
        line = blank + str(self.xs_tot)+'\n'
        file.write(line)
        file.write("MT 1452\n")
        line = blank + '0' +'\n'
        file.write(line)
        file.write("MT 1018\n")
        line = blank + '0' +'\n'
        file.write(line)
        for i in xrange(0,self.L_max+1) :
            line = MT_2501 + str(i) + '\n'
            file.write(line)
            file.write(sink)
            line = blank + str(self.xs_sca[i]) + '\n'
            file.write(line)
        file.close()

        print "The width of the medium is %f" %self.width

