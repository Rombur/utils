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
n_division_z = 100
delta_z = 0.05
n_g_groups = 15
n_e_groups = 25
n_p_groups = 0
x_position = [51.25, 53.75] 
y_position = [51.25, 53.75]
n_moments = 64

pdt = dose.dose(filename, xs_filename, n_processors, n_division_z, delta_z, 
        n_g_groups, n_e_groups, n_p_groups, x_position, y_position, n_moments)

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
