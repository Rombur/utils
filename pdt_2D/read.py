#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2011-06-20 16:10:24.833443
#----------------------------------------------------------------------------#

"""Read the inputs files"""

import read_output
import convergence

a = 'BiLD.output_1'
b = 'BiLD.output_2'
c = 'BiLD.output_4'
d = 'BiLD.output_8'
ref = 'BiLD.output_32'
a_out = 'out_10.pckl'
b_out = 'out_20.pckl'
c_out = 'out_40.pckl'
d_out = 'out_80.pckl'
ref_out = 'out_320.pckl'

a_ro = read_output.read_output(a,a_out)
b_ro = read_output.read_output(b,b_out)
c_ro = read_output.read_output(c,c_out)
d_ro = read_output.read_output(d,d_out)
ref_ro = read_output.read_output(ref,ref_out)

a_ro.read_file()
a_ro.print_flux()

b_ro.read_file()
b_ro.print_flux()

c_ro.read_file()
c_ro.print_flux()

d_ro.read_file()
d_ro.print_flux()

ref_ro.read_file()
ref_ro.print_flux()

list_file=[a_out,b_out,c_out,d_out]
conv = convergence.convergence(list_file,ref_out)
conv.compute_order()
