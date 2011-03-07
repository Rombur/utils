#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2010-11-19 16:03:01.097206
#----------------------------------------------------------------------------#

import dose

filename = "input_cepxs"
xs_filename = "sig_cepxs.data"
n_processors = 4
n_division = 100
delta = 0.05
n_g_groups = 15
n_e_groups = 25
n_p_groups = 0
x_position = [51.25, 53.75] 
y_position = [51.25, 53.75]
z_position = []
n_moments = 168
moments = [0,30,132]
groups = [0,n_g_groups-1,n_g_groups,n_g_groups+n_e_groups-1]
method = "old"

pdt = dose.dose(filename, xs_filename, n_processors, n_division, delta, 
        n_g_groups, n_e_groups, n_p_groups, x_position, y_position,
        z_position, n_moments, moments, groups)
              
if method == "old" :
    print "reading the flux"
    pdt.read_flux()
    print "writing the flux"
    pdt.print_flux()
    print "reading the cross section"
    pdt.read_xs()
    
    out = " "
    while out not in ["quit","q","out"] :
        out = raw_input("Give a position : ")
        if out not in ["quit","q","out"] :
            pos = float(out)
            pdt.compute_dose(pos)
elif method == "new":
    print "reading the output file"
    pdt.read_dose()
    print "printing the dose"
    pdt.print_dose()
    print "plotting the dose"
    pdt.plot_dose()
else :
    print "wrong method"
