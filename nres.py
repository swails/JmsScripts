#!/usr/bin/env python

# Find the number of residues in a given prmtop

from sys import stderr, stdout, exit, argv
from readparm import amberParm

def usage():
   print >> stderr, "{0} <prmtop>".format(argv[0])
   exit()

if len(argv) != 2:
   usage()

x = amberParm(argv[1])

if x.exists:
   print x.ptr(argv[0].strip('.py'))

else:
   print >> stderr, '{0} does not exist!'.format(argv[1])
