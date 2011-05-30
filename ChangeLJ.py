#!/usr/bin/env python

from chemistry.amber.readparm import amberParm
from optparse import OptionParser
import sys

parser = OptionParser()

parser.add_option('-p', '--prmtop', type='str', dest='prmtop',
                  help='Input topology file')

parser.add_option('-n', '--new-prmtop', type='str', 
                  dest='new_prmtop', help='New prmtop name. ' +
                  'Defaults to original prmtop if omitted')

parser.add_option('-a', '--atom1', type='str', dest='atom_1',
                  help='Atom Type 1')

parser.add_option('-b', '--atom2', type='str', dest='atom_2',
                  help='Atom Type 2')

parser.add_option('-e', '--epsilon', type='float', dest='eps',
                  help='Combined well depth (epsilon)')

parser.add_option('-s', '--sigma', type='float', dest='sig',
                  help='Combined well radius (sigma)')

(opt, args) = parser.parse_args()

# Sanity checks
if not opt.prmtop or not opt.atom_1 or not opt.atom_2 or \
   not opt.eps or not opt.sig:
   parser.print_help()
   sys.exit()

if not opt.new_prmtop: opt.new_prmtop = opt.prmtop

parm = amberParm(opt.prmtop)
parm.overwrite = True

if not opt.atom_1 in parm.LJ_types or not opt.atom_2 in parm.LJ_types:
   print >> sys.stderr, 'Error: Atom type(s) not found. Types must be in'
   print >> sys.stderr, parm.LJ_types.keys()
   sys.exit()

# Make sure that atom type 1 comes first
if parm.LJ_types[opt.atom_1] > parm.LJ_types[opt.atom_2]:
   holder = opt.atom_1
   opt.atom_1 = opt.atom_2
   opt.atom_2 = holder

# Find the index of where the pairs of the first atom type starts
atm_1_idx = parm.LJ_types[opt.atom_1] - 1
atm_2_idx = parm.LJ_types[opt.atom_2] - 1

# Find the atom1 - atom1 interaction (adjusting for indexing from 0)
term_idx = parm.parm_data['NONBONDED_PARM_INDEX'][
                parm.parm_data['NTYPES']*(atm_1_idx - 1) + atm_2_idx] - 1

# Now change the ACOEF and BCOEF arrays, assuming the proper combining
# rules
parm.parm_data['LENNARD_JONES_ACOEF'][term_idx] = \
   opt.eps * opt.sig ** 12

parm.parm_data['LENNARD_JONES_BCOEF'][term_idx] = \
   2 * opt.eps * opt.sig ** 6

parm.writeParm(opt.new_prmtop)

print >> sys.stdout, 'Your topology file is complete (%s).' % opt.new_prmtop
