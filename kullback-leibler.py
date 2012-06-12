#!/usr/bin/env python

""" 
This program performs Kullback-Leibler analysis on a time series of RMSD data
"""

from __future__ import division

import dataset as ds
from optparse import OptionParser, OptionGroup
import sys
import signal
from os import path

parser = OptionParser()
group = OptionGroup(parser, "Input Options", 
                    "Options pertaining to the input data")
group.add_option('-i', '--input-data', metavar='FILE', dest='infile',
                 default=None, help='Input file with RMSD data')
group.add_option('-c', '--column', metavar='INT', type='int', dest='col',
                 help='Column with RMSD data (Default=%default)', default=2)
parser.add_option_group(group)
group = OptionGroup(parser, 'Output Options',
                    'Options pertaining to the output data')
group.add_option('-o', '--output', metavar='FILE', dest='outfile', 
                 default=None, help='Output file')
group.add_option('--interval', type='int', default=500, metavar='INT', 
                 help='Number of snapshots between adjacent histogram samples' +
                 ' (Default=%default)', dest='interval')
parser.add_option_group(group)
group = OptionGroup(parser, "Histogram Options", "Options pertaining to " +
                    "the construction of the RMSD histograms")
group.add_option('--min', dest='hmin', default=None, type='float',
                 metavar='FLOAT', help='Minimum of the histograms. Default ' +
                 'is minimum of data set.')
group.add_option('--max', dest='hmax', default=None, type='float',
                 metavar='FLOAT', help='Maximum of the histograms. Default ' +
                 'is maximum of data set.')
group.add_option('--bins', dest='bins', default=None, type='int',
                 metavar='INT', help='Number of bins in histogram. Cannot ' +
                 'be specified with --spacing')
group.add_option('--spacing', dest='spacing', default=None, type='int',
                 metavar='FLOAT', help='Spacing between bins in histogram. ' +
                 'Cannot be specified with --bins')
parser.add_option_group(group)
opt, arg = parser.parse_args()

# Set up signal processor
def interrupt_handler(*args, **kwargs):
   print '%s interrupted.' % (path.split(sys.argv[0])[1])
   parser.print_help()
   sys.exit(1)

signal.signal(signal.SIGINT, interrupt_handler)

if opt.infile is None:
   infile = sys.stdin
else:
   try:
      infile = open(opt.infile, 'r')
   except IOError:
      print >> sys.stderr, 'Could not open %s for reading!' % opt.infile
      sys.exit(1)

if opt.outfile is None:
   outfile = sys.stdout
else:
   try:
      outfile = open(opt.outfile, 'w')
   except IOError:
      print >> sys.stderr, 'Could not open %s for writing!' % opt.outfile
      sys.exit(1)

data = ds.load_from_file(infile, opt.col)

if opt.spacing is not None and opt.bins is not None:
   print >> sys.stderr, 'Error: Cannot specify --bins and --spacing!'
   sys.exit(1)

data.set_hist_params(hmin=opt.hmin, hmax=opt.hmax, nbins=opt.bins,
                     spacing=opt.spacing, norm=True)

data.KullbackLeibler(opt.interval, outfile)
