#!/usr/bin/env python

from readparm import amberParm
import sys

def printusage():
   """ Prints the usage statement """
   print >> sys.stdout, 'Usage: getResidue.py <prmtop> <atom #>'
   sys.exit()


if len(sys.argv) != 3:
   printusage()

parm = amberParm(sys.argv[1])

if not parm.valid:
   print >> sys.stderr, 'Error: Invalid topology file (%s)' % str(parm)
   printusage()

atom = int(sys.argv[2])

for i in range(parm.ptr('nres')):
   if parm.parm_data['RESIDUE_POINTER'][i] > atom:
      print 'Atom %d (%s) in prmtop %s is in residue %s-%s' % (atom, 
           parm.parm_data['ATOM_NAME'][atom-1], str(parm), i, parm.parm_data['RESIDUE_LABEL'][i])
      break


