#!/usr/bin/env python

from readparm import amberParm
import sys

def usage():
   print 'Usage: gimmeFrcmod.py <prmtop> <frcmod>'
   sys.exit()

if len(sys.argv) != 3:
   usage()

print 'Reading the prmtop...'

prmtop = amberParm(sys.argv[1])

if not prmtop.exists:
   print 'Error: %s does not exist!' % sys.argv[1]
   usage()

print 'gimmeing your frcmod...'
prmtop.frcmod(sys.argv[2])

print 'Fine. take it. (%s)' % sys.argv[2]
