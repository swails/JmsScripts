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
group.add_argument('-c' ,'--columns', nargs=2, type=int, metavar='INT',
                   dest='columns', default=[0,1], help='''Which columns to
                   extract data from.  Default is first and second columns''')
group = parser.add_argument_group('Output Options')
group.add_argument('-xr', '--xrange', metavar='FLOAT', nargs=2, type=float,
                   dest='xrange', default=None, help='''Range of output data to
                   print in X- and Y- dimensions. If you are printing a torsion,
                   use --torsion below to set up a periodic range from -180 to
                   180''')
group.add_argument('-yr', '--yrange', metavar='FLOAT', nargs=2, type=float,
                   dest='yrange', default=None, help='''Same as -xr/--xrange,
                   except in the Y-dimension.''')
group.add_argument('-xt', '--xtorsion', default=False, dest='xtorsion',
                   action='store_true', help='''The X-dimension is a torsion, so
                   plot from -180 to 180 degrees''')
group.add_argument('-yt', '--ytorsion', default=False, dest='ytorsion',
                   action='store_true', help='''The Y-dimension is a torsion, so
                   plot from -180 to 180 degrees''')
group.add_argument('-res', '--resolution', metavar='INT', nargs=2, type=int,
                   dest='res', default=[100,100], help='''Number of points to
                   output in the X and Y dimensions, respectively. Defaults to
                   100 x 100''')
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
   return x + buffer * r, n - buffer * r

opt = parser.parse_args()

if opt.plot:
   warnings.warn('matplotlib plotting is not yet supported.')

if opt.xrange is None and not opt.xtorsion:
   print('Using range of data as default X-range (with 10% buffer)')
   xmin, xmax = None, None
elif opt.xtorsion:
   print('X-dimension is torsion. Using -180 to 180 range.')
   xmin, xmax = -180.0, 180.0
else:
   print('X-min: %g; X-max: %g' % opt.xrange)
   xmin, xmax = opt.xrange

if opt.yrange is None and not opt.ytorsion:
   print('Using range of data as default Y-range (with 10% buffer)')
   ymin, ymax = None, None
elif opt.ytorsion:
   print('Y-dimension is torsion. Using -180 to 180 range.')
   ymin, ymax = -180.0, 180.0
else:
   print('Y-min: %g; Y-max: %g' % opt.yrange)

if opt.input is None:
   infile = sys.stdin
else:
   infile = open(opt.input, 'r')

xvals = array('d')
yvals = array('d')
for line in infile:
   try:
      words = line.split()
      x = float(words[opt.columns[0]])
      y = float(words[opt.columns[1]])
   except ValueError:
      continue
   except IndexError:
      continue

   xvals.append(x)
   yvals.append(y)

# Determine ranges if we asked for default
if xmin is None:
   xmin, xmax = buffered_range(xmin, xmax, 0.05)
if ymin is None:
   ymin, ymax = buffered_range(ymin, ymax, 0.05)

# Now convert to numpy arrays
data = np.zeros(shape=(2,len(xvals)), order='F')
for i in range(len(xvals)):
   data[0][i] = xvals[i]
   data[1][i] = yvals[i]

# Get rid of the old data
del xvals, yvals

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
spacing = [(xmax - xmin) / opt.res[0], (ymax - ymin) / opt.res[1]]

outfile.write('#             X             Y           KDE\n')
for i in range(opt.res[0]):
   for j in range(opt.res[1]):
      xval = xmin + spacing[0] * i
      yval = ymin + spacing[1] * j
      outfile.write('%13.7E %13.7E %13.7E\n' % (xval, yval,
                              kernel.evaluate((xval, yval))
                                             )
                   )
   outfile.write('\n')

outfile.close()
