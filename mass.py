#!/usr/bin/env python

####################################################
#                                                  #
# This utility prints the total mass of a em   #
# described by the prmtop.                         #
#                                                  #
#    Written by Jason Swails, 08/02/2010           #
#                                                  #
####################################################

from sys import stderr, stdout, exit, argv
from chemistry.amber.readparm import amberParm

if len(argv) != 2 or '-help' in argv[1] or argv[1] == '-h':
   print 'mass.py <prmtop>'
   exit()

prm = amberParm(argv[1])

if prm.exists:
   print >> stdout, 'The mass of ' + argv[1] + ' is ' + str(prm.totMass()) + ' g/mol.'
else:
   print >> stderr, 'Prmtop {0} does not exist!'.format(argv[1])
