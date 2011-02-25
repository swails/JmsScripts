#!/usr/bin/env python

""" Changes the protonation state of a given residue in a given prmtop """

from cpin_data import getData
from chemistry.amber.readparm import amberParm
import sys, os

def printusage():
   """ Prints the usage statement and quits """
   print "ChangeState.py <prmtop> <residue> <state> <new_prmtop>"
   sys.exit(0)

if len(sys.argv) != 5 or sys.argv[1] in ['--help','-help','-h','--h']:
   printusage()

prmtop = amberParm(sys.argv[1])

if not prmtop.valid:
   print >> sys.stderr, 'Error: %s is not a valid prmtop!' % sys.argv[1]
   printusage()

res = int(sys.argv[2]) - 1
data = getData(prmtop.parm_data['RESIDUE_LABEL'][res],5)
starting_point = prmtop.parm_data['RESIDUE_POINTER'][res] - 1
state = int(sys.argv[3])

for i in range(2,len(data[state])):
   atnum = starting_point + i - 2
   print "Changing %-4s from charge %8.4f to %8.4f" % \
      ( prmtop.parm_data['ATOM_NAME'][atnum], prmtop.parm_data['CHARGE'][atnum],
        data[state][i] )
   prmtop.parm_data['CHARGE'][atnum] = data[state][i]

prmtop.writeParm(sys.argv[4])

print 'Done!'
