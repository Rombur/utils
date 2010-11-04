#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2010-11-04 09:27:35.481599
#----------------------------------------------------------------------------#

import sys
sys.path.append("/home/turcksin/python")
import utils
import rutherford

utils = utils.utils()

print "Script to compute the cross sections and the width for the Rutherford test problem."

test = raw_input("Choose a test (between 1 and 9) : ")
test = int(test)
if test < 1 or test > 9 :
    utils.abort("The test should be between 1 and 9.")
else :
    tc = raw_input("Do want to use the transport correction ? ")
    if tc in ["true", "True", "1", "y", "yes", "yeah"] :
        tc = True
    else :
        tc = False
    ruth = rutherford.rutherford(test,tc)
    ruth.compute_xs()
    ruth.compute_width()
    filename = raw_input("Give the name of the file where you want to append the cross section : ")
    ruth.write(filename)
