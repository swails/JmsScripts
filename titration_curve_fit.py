#!/usr/bin/env python

from __future__ import division

""" 
This program will calculate the chi-squared of a best-fit titration curve to
the Hill equation that it was fitted to
"""

from optparse import OptionParser
from scipy.optimize import curve_fit
import numpy as np
import sys, math

def f(ph, pka, n):
   return 1.0/(1.0+10.0**(n*(pka-ph)))

if __name__ == '__main__':
   parser = OptionParser()
   parser.add_option('-i', '--input', dest='input_file', metavar='FILE',
                     help='Input file with titration data. Default stdin.',
                     default=None)
   parser.add_option('-d', '--deprotonated', dest='protonated', default=True,
                     action='store_false', help='The input data is fraction ' +
                     'deprotonated. NOT default behavior.')
   opt, arg = parser.parse_args()

   if not opt.input_file:
      infile = sys.stdin
   else:
      try:
         infile = open(opt.input_file, 'r')
      except IOError:
         print >> sys.stderr, 'Could not open %s for reading' % opt.input_file
         sys.exit(1)
   
   lines = infile.readlines()

   xdata = np.zeros(len(lines))
   ydata = np.zeros(len(lines))

   for i, line in enumerate(lines):
      xdata[i] = float(line.split()[0])  # pH
      ydata[i] = float(line.split()[1])  # Fraction deprotonated
      if opt.protonated: ydata[i] = 1.0 - ydata[i]

   # Generate an initial guess, which is just the average of the HH pKas
   avg = 0.0
   for i, fd in enumerate(ydata):
      if fd > 0 and fd < 1: 
         avg += xdata[i] - math.log10(fd / (1-fd))
   avg /= len(lines)

   params, cov = curve_fit(f, xdata, ydata, p0=(avg,1))

   print params

   sum2 = 0.0
   for i in range(len(xdata)):
      sum2 += (ydata[i] - f(xdata[i], params[0], params[1])) ** 2

   print sum2
