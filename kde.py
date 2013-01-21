#!/usr/bin/env python
from __future__ import division

from array import array
from argparse import ArgumentParser
import numpy as np
from scipy import stats
from signal import signal, SIGINT
import sys
import warnings

parser = ArgumentParser(version='0.1a')

group = parser.add_argument_group('Files')
group.add_argument('-i', '--input', dest='input', metavar='FILE', help='''Input
                   file with all data. Default reads from stdin.''',
                   default=None)
group.add_argument('-o', '--output', dest='output', metavar='FILE',
                   help='Output file with all data. Defaults to stdout.',
                   default=None)
group = parser.add_argument_group('Input Options')
group.add_argument('-c' ,'--column', type=int, metavar='INT', dest='column',
                   default=0, help='''Which column to extract data from.
                   Default is first column''')
group = parser.add_argument_group('Output Options')
group.add_argument('-r', '--range', metavar='FLOAT', nargs=2, type=float,
                   dest='range', default=None, help='''Range of output data to
                   print. If you are printing a torsion, use --torsion below to
                   set up a periodic range from -180 to 180''')
group.add_argument('-t', '--torsion', default=False, dest='torsion',
                   action='store_true', help='''The data is from a torsion, so
                   plot from -180 to 180 degrees''')
group.add_argument('-res', '--resolution', metavar='INT', type=int, dest='res',
                   default=100, help='Number of points to output. Default 100.')
group.add_argument('-b', '--bandwidth', default=None, type=float,
                   metavar='FLOAT', dest='bandwidth', help='''Kernel bandwidth
                   to use. Defaults to Scott's choice.''')
group.add_argument('--plot', dest='plot', action='store_true', default=False,
                   help='''Show surface plot using matplotlib. Default is not
                   to.''')


# Set up signal handler to print help
signal(SIGINT, lambda *args, **kwargs: parser.print_help())

def buffered_range(n, x, buffer=0.1):
   """ Returns a buffered range based on input miN and maX """
   r = x - n
   return n - buffer * r, x + buffer * r

opt = parser.parse_args()

if opt.plot:
   warnings.warn('matplotlib plotting is not yet supported.')

if opt.range is None and not opt.torsion:
   print('Using range of data as default range (with 10% buffer)')
   xmin, xmax = None, None
elif opt.torsion:
   print('Data is from a torsion. Using -180 to 180 range.')
   xmin, xmax = -180.0, 180.0
else:
   print('X-min: %g; X-max: %g' % opt.xrange)
   xmin, xmax = opt.xrange

if opt.input is None:
   infile = sys.stdin
else:
   infile = open(opt.input, 'r')

data = array('d')
for line in infile:
   try:
      words = line.split()
      data.append(float(line.split()[opt.column]))
   except ValueError:
      continue
   except IndexError:
      continue

# Determine ranges if we asked for default
if xmin is None:
   xmin, xmax = buffered_range(min(data), max(data), 0.05)

# Now convert to numpy array
data = np.asarray(data)

# Now pass it in and get a KDE
try:
   kernel = stats.gaussian_kde(data, bw_method=opt.bandwidth)
except TypeError:
   kernel = stats.gaussian_kde(data)
   if opt.bandwidth is not None:
      kernel.factor = opt.bandwidth

# Output the results in a gnuplot-readable way
if opt.output is None:
   outfile = sys.stdout
else:
   outfile = open(opt.output, 'w')

print ('The bandwidth is %s' % kernel.factor)
spacing = (xmax - xmin) / opt.res

outfile.write('#             X           KDE\n')
for i in range(opt.res):
   xval = xmin + spacing * i
   outfile.write('%13.7E %13.7E\n' % (xval, kernel.evaluate((xval))))

outfile.close()
