#!/usr/bin/env python

# Find the number of residues in a given prmtop

from sys import stderr, stdout, exit, argv
from readparm import amberParm
from os import path

def usage():
   print >> stderr, "{0} <prmtop>".format(path.split(argv[0])[1])
   exit()

if len(argv) != 2:
   usage()

x = amberParm(argv[1])

if x.exists:
   print x.ptr(path.split(argv[0])[1].strip('.py'))

else:
   print >> stderr, '{0} does not exist!'.format(argv[1])
