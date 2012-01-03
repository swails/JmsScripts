#!/usr/bin/env python

from mdcrd import AmberTraj, RmsdData
from optparse import OptionParser
from utilities import which
from subprocess import Popen, PIPE
import sys, os

parser = OptionParser(usage='%prog [Options] mdcrd1 [mdcrd2 [... ] ]')
parser.add_option('-p', '--prmtop', dest='prmtop', metavar='FILE',
                  default='prmtop', help='Prmtop associated with trajectories')
parser.add_option('-r', '--rmsd-file', dest='rmsfile', metavar='FILE',
                  default=None, help='Data file to dump RMSDs to')
parser.add_option('-b', '--rmsd-bin-file', dest='rmsbin', metavar='FILE',
                  default=None, help='Data file to dump RMSD bins to')
parser.add_option('-f', '--reference', dest='reffile', metavar='FILE',
                  default=None, help='Reference file for RMSDs')
parser.add_option('-m', '--mask', dest='mask', default='@CA',
                  help='Mask to calculate RMSd for. Default [%default]')
opt, args = parser.parse_args()

if not args or not opt.rmsfile or not os.path.exists(opt.prmtop):
   print >> sys.stderr, 'Error: Missing key CL arguments or missing topology!'
   parser.print_help()
   sys.exit(1)

if opt.rmsbin:
   binner = which('1Dbinning.py')
   if not binner:
      print >> sys.stderr, 'Error: Binning requires 1Dbinning.py!'
      sys.exit(1)

trajs = AmberTraj(opt.prmtop, args)

trajs.rmsd(mask=opt.mask, outfile=opt.rmsfile, ref=opt.reffile)
trajs.run()

if opt.rmsbin:
   process = Popen([binner, '-f', opt.rmsfile, '-o', opt.rmsbin, '-n', '-c', 2])
   if process.wait():
      print >> sys.stderr, 'Error: Binning program (%s) failed!' % binner
      sys.exit(1)
   else:
      print 'RMS calculation finished'
sys.exit(0)
