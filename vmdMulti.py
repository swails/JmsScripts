#!/usr/bin/env python

# This script will load all of the specified files as new molecules

import sys

def printusage():
   print 'vmdMulti.py pdb1 pdb2 pdb3 ... '
   sys.exit()

if len(sys.argv) < 2:
   printusage()

string = 'vmd '

for x in range(1,len(sys.argv)):
   string += '-f %s' % (sys.argv[x])

print '%s\n' % (string)
