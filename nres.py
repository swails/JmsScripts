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

if argv[0].endswith('natom.py'):
   ptr = "NATOM"
elif argv[0].endswith('nres.py'):
   ptr = "NRES"

if x.exists:
   print x.ptr(ptr)

else:
   print >> stderr, '{0} does not exist!'.format(argv[1])
