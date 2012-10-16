#!/usr/bin/env python
""" 
This script will take 2 prmtop files that have different charge sets and create
a set of interpolated prmtops that vary in the charges.  It does not yet
interpolate VDW radii, so be careful.

By Jason Swails, Updated 04/18/2011
"""

from chemistry.amber.readparm import amberParm
from MMPBSA_mods.commandline_parser import OptionParser
import sys
import os

debug_printlevel = 0

def excepthook(extype, exval, tb):
   import traceback
   global debug_printlevel
   if debug_printlevel > 0: traceback.print_tb(tb)
   print >> sys.stderr, '%s: %s' % (extype.__name__, exval)
   sys.exit(1)

sys.excepthook = excepthook

# Set up the option parser
parser = OptionParser()

parser.addOption('-p1', 'prmtop1', help='Prmtop at Lambda = 0')
parser.addOption('-p2', 'prmtop2', help='Prmtop at Lambda = 1')
parser.addOption('-prefix', 'prefix', help='Prefix for new topology file names')
parser.addOption('-n', 'number', help='Number of TOTAL prmtops you want, ' +
                 'including endpoints', default=5)
parser.addOption('-groupfile', 'groupfile', help='Optionally print out a ' +
                 'skeleton groupfile for running H-REMD simulations with ' +
                 'sander.MPI or pmemd.MPI')
parser.addOption('-debug', 'debug', help='0: Suppress tracebacks, ' + 
                 '1: print tracebacks', default=0)

opt = parser.Parse()

debug_printlevel = opt.debug
opt.number = int(opt.number)

if not opt.prmtop1 or not opt.prmtop2:
   parser.print_help()

# Set up the prmtop objects
print >> sys.stdout, 'Reading topology files'
parm1 = amberParm(opt.prmtop1)
parm2 = amberParm(opt.prmtop2)

print >> sys.stdout, 'Validating topology files'
if not (parm1.valid and parm2.valid):
   print >> sys.stderr, 'Error: Problems with input prmtop files!'
   parser.print_help()

if parm1.parm_data['RESIDUE_LABEL'] != parm2.parm_data['RESIDUE_LABEL']:
   print >> sys.stderr, 'Error: You must have the same number of residues ' + \
                        ' in each prmtop!'
   sys.exit(1)
if parm1.parm_data['ATOM_NAME'] != parm2.parm_data['ATOM_NAME']:
   print >> sys.stderr, 'Error: Atom sequences must be the same in each prmtop!'
   sys.exit(1)

# Make copies so we don't change either set when we change our parm1 charges

print >> sys.stdout, 'Copying charges...'
charges_1 = parm1.parm_data['CHARGE'][:]
charges_2 = parm2.parm_data['CHARGE'][:]

print >> sys.stdout, '\nCreating %d new topology files\n' % opt.number
for i in range(opt.number):
   lmda = 1.0 / (opt.number - 1) * i
   for j in range(parm1.ptr('natom')):
      parm1.parm_data['CHARGE'][j] = charges_2[j]*lmda + (1-lmda)*charges_1[j]
   
   print >> sys.stdout, '\tWriting %s' % (opt.prefix + '.%.2f' % lmda)
   parm1.writeParm(opt.prefix + '.%.2f' % lmda)

print >> sys.stdout, '\nDone writing topology files!'

if opt.groupfile:
   print >> sys.stdout, '\nWriting groupfile %s for amber simulations' % \
            opt.groupfile

   grpfl = open(opt.groupfile, 'w')
   for i in range(opt.number):
      lmda = 1.0 / (opt.number - 1) * i
      grpfl.write('# Replica at lambda = %6.4f\n' % lmda)
      grpfl.write(('-O -i <MDIN> -o <MDOUT> -p %s.%.2f -c <INPCRD> -r <RESTRT> '
                  + '-x <MDCRD>\n\n') % (opt.prefix, lmda))
   grpfl.close()

   print >> sys.stdout, 'Done writing %s' % opt.groupfile

print >> sys.stdout, '%s execution complete!' % (os.path.split(sys.argv[0])[1])
