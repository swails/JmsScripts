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
try:
   import matplotlib.pyplot as plt
   hasplt = True
except:
   hasplt = False

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
   parser.add_option('-t', '--title', dest='title', default='',
                     metavar='STRING',help='Title for the plot (appended with '
                     '"Titration Curve")')
   parser.add_option('--no-plot', dest='plot', default=True,
                     action='store_false', help='Do not plot the titration '
                     'curve.')
   parser.add_option('-c', '--column', dest='col', default=1, type='int',
                     metavar='INT', help='Column that titration data is in. '
                     'First column, column 0, must have the pH values. Default '
                     '(%default)')
   opt, arg = parser.parse_args()

   if not opt.input_file:
      infile = sys.stdin
   else:
      try:
         infile = open(opt.input_file, 'r')
      except IOError:
         print >> sys.stderr, 'Could not open %s for reading' % opt.input_file
         sys.exit(1)
   
   hasplt = hasplt and opt.plot

   lines = infile.readlines()

   xdata = np.zeros(len(lines))
   ydata = np.zeros(len(lines))

   for i, line in enumerate(lines):
      if line.strip().startswith('#'): continue
      xdata[i] = float(line.split()[0])  # pH
      ydata[i] = float(line.split()[opt.col])  # Fraction deprotonated
      if opt.protonated: ydata[i] = 1.0 - ydata[i]

   # Generate an initial guess, which is just the average of the HH pKas
   avg = 0.0
   ncalls = 0
   for i, fd in enumerate(ydata):
      if fd > 0 and fd < 1: 
         avg += xdata[i] - math.log10(fd / (1-fd))
         ncalls += 1
   avg /= ncalls

   params, cov = curve_fit(f, xdata, ydata, p0=(avg,1))

   sum2 = 0.0
   for i in range(len(xdata)):
      sum2 += (ydata[i] - f(xdata[i], params[0], params[1])) ** 2


   if not hasplt:
      print 'pKa = %f' % params[0]
      print 'n   = %f' % params[1]
      print 'RSS = %f' % sum2
      sys.exit(0)
   
   if opt.title:
      opt.title = opt.title.strip() + ' '

   # Plot me
   fcn_x = np.linspace(min(xdata)-1, max(xdata)+1, 1000)
   fcn_y = [f(i, params[0], params[1]) for i in fcn_x]
   plt.axis([min(xdata)-1, max(xdata)+1, 0, 1])
   plt.xlabel('pH')
   plt.ylabel('Fraction Deprotonated')
   plt.title('%sTitration Curve' % opt.title)
   plt.text(min(xdata)-0.9, 0.4, 'pKa = %.2f\nn    = %.2f\nRSS = %.4e' %
            (params[0], params[1], sum2))
   plt.grid(True)
   myplot = plt.plot(xdata, ydata, 'bo', fcn_x, fcn_y, 'r-')
   plt.show()
