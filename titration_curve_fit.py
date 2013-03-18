#!/usr/bin/env python
""" 
This program will calculate the chi-squared of a best-fit titration curve to
the Hill equation that it was fitted to
"""

from __future__ import division

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
   
   if opt.protonated:
      transform = lambda x: x
   else:
      transform = lambda x: 1 - x
   hasplt = hasplt and opt.plot

   ress = [x for i, x in enumerate(infile.readline().split()) if i > 0]
   data = np.loadtxt(infile).transpose()
   avgs = np.zeros(data.shape[0])
   ncalls = np.zeros(data.shape[0], dtype=np.int8)
   sum2s = np.zeros(data.shape[0])
   paras = np.zeros(shape=(data.shape[0], 2))

   # Generate initial guesses, which is just the average of the HH pKas
   for j in range(1, data.shape[0]):
      for i, x in enumerate(data[j]):
         fd = transform(x)
         if fd > 0 and fd < 1: 
            avgs[j] += data[0][i] - math.log10(fd / (1-fd))
            ncalls[j] += 1
      avgs[j] /= ncalls[j]

      params, cov = curve_fit(f, data[0], data[j], p0=(avgs[j],1))

      for i in range(len(data[0])):
         sum2s[j] += (data[j][i] - f(data[0][i], params[0], params[1])) ** 2

      paras[j][0], paras[j][1] = params[0], params[1]
      if not hasplt:
         print '%-20s: pKa = %8.3f n = %8.3f RSS = %8.3f' % (ress[j], params[0],
                     params[1], sum2s[j])
   
   if not hasplt:
      sys.exit(0)

   if opt.title:
      opt.title = opt.title.strip() + ' '

   # Plot me
   fig = plt.figure(1)
   ncols = max(data.shape[0] - 1, 3)
   nrows = max(1, int(np.ceil(data.shape[0] / 3)))
   for j in range(1, data.shape[0]):
      ax = fig.add_subplot(nrows, ncols, j)
      ax.axis([np.min(data[0])-1, np.max(data[0])+1, 0, 1])
      ax.xlabel('pH')
      ax.ylabel('Fraction Deprotonated')
      ax.title('%s Titration Curve' % ress[j])
      fcn_x = np.linspace(np.min(data[0])-1, np.max(data[0])+1, 1000)
      fcn_y = [f(i, paras[j][0], paras[j][1]) for i in fcn_x]
      ax.text(np.min(data[0])-0.9, 0.4, 'pKa = %.2f\nn    = %.2f\nRSS = %.4e' %
               (paras[j][0], paras[j][1], sum2s[j]))
      ax.grid(True)
      ax.plot(data[0], data[j], 'bo', fcn_x, fcn_y, 'r-')
   plt.show()
