#!/usr/bin/env python

####################################################
#                                                  #
# This utility prints the total mass of a system   #
# described by the prmtop.                         #
#                                                  #
#    Written by Jason Swails, 08/02/2010           #
#                                                  #
####################################################


import sys
from readparm import amberParm

if len(sys.argv) != 2 or '-help' in sys.argv[1] or sys.argv[1] == '-h':
   print 'mass.py <prmtop>'
   sys.exit()

prm = amberParm(sys.argv[1])

if prm.prm_exists:
   print >> sys.stdout, 'The mass of ' + sys.argv[1] + ' is ' + str(prm.mass()) + ' g/mol.'
