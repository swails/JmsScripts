#!/usr/bin/env python

import sys, os
from chemistry.amber.readparm import AmberParm
from MMPBSA_mods.alamdcrd import GlyMutantMdcrd
from optparse import OptionParser, OptionGroup

epilog = 'This will mutate a trajectory with a glycine mutant'

parser = OptionParser(epilog=epilog)
group = OptionGroup(parser, 'Input Files', 'Files you must provide')
group.add_option('-p', '--prmtop', dest='norm_prm', metavar='FILE',
                  default=None, help='Original prmtop file')
group.add_option('-m', '--mutant-prmtop', dest='mut_prm', metavar='FILE',
                  default=None, help='Mutant prmtop file')
group.add_option('-y', '--inptraj', dest='inptraj', metavar='FILE',
                  default=None, help='Input trajectory file you wish to '
                  'mutate.')
parser.add_option_group(group)
group = OptionGroup(parser,'Output Files', 'Files that are made by this script')
group.add_option('-o', '--outtraj', dest='outtraj', metavar='FILE',
                  default='mutant.mdcrd', help='Name of mutant trajectory '
                  'file to generate. Default %default')
parser.add_option_group(group)

opt, arg = parser.parse_args()

if opt.norm_prm is None or opt.mut_prm is None:
   sys.stderr.write('Please provide a normal and mutated topology file\n')
   parser.print_help()
   sys.exit(1)

if opt.inptraj is None:
   sys.stderr.write('Please provide a trajectory to mutate!\n')
   parser.print_help()
   sys.exit(1)

norm_prm = AmberParm(opt.norm_prm)
mut_prm = AmberParm(opt.mut_prm)

mut_traj = GlyMutantMdcrd(opt.inptraj, norm_prm, mut_prm)
mut_traj.MutateTraj(opt.outtraj)

print 'FINALLY DONE!'
