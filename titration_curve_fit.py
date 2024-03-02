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
   import matplotlib
   matplotlib.use('TkAgg')
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
                     action='store_false', help='The input data is fraction '
                     'deprotonated. NOT default behavior.')
   parser.add_option('-t', '--title', dest='title', default='',
                     metavar='STRING',help='Title for the plot (appended with '
                     '"Titration Curve")')
   parser.add_option('--no-plot', dest='plot', default=True,
                     action='store_false', help='Do not plot the titration '
                     'curve.')
   opt, arg = parser.parse_args()

   if not opt.input_file:
      infile = sys.stdin
   else:
      try:
         infile = open(opt.input_file, 'r')
      except IOError:
         sys.exit('Could not open %s for reading' % opt.input_file)
   
   if opt.protonated:
      transform = lambda x: 1 - x
   else:
      transform = lambda x: x
   hasplt = hasplt and opt.plot

   ress = [x for i, x in enumerate(infile.readline().split()) if i > 0]
   data = np.loadtxt(infile).transpose()
   avgs = np.zeros(data.shape[0])
   ncalls = np.zeros(data.shape[0], dtype=np.int8)
   sum2s = np.zeros(data.shape[0])
   paras = np.zeros(shape=(data.shape[0], 2))

   # Generate initial guesses, which is just the average of the HH pKas
   all0, all1 = True, True
   for j in range(1, data.shape[0]):
      for i, x in enumerate(data[j]):
         fd = transform(x)
         if fd > 0 and fd < 1: 
            all0, all1 = False, False
            avgs[j] += data[0][i] - math.log10(fd / (1-fd))
            ncalls[j] += 1
         elif fd == 0:
            all1 = False
         elif fd == 1:
            all0 = False

      # If we have no initial guesses
      if ncalls[j] > 0:
         avgs[j] /= ncalls[j]
         params, cov = curve_fit(f, data[0], data[j], p0=(avgs[j],1))

      for i in range(len(data[0])):
         sum2s[j] += (data[j][i] - f(data[0][i], params[0], params[1])) ** 2

      if ncalls[j] == 0:
         params = [0,0]
         sum2s[j] = 0

      paras[j][0], paras[j][1] = params[0], params[1]
      if not hasplt and ncalls[j] > 0:
         print('%-20s: pKa = %8.3f n = %8.3f RSS = %10.3e' % (ress[j], params[0],
                     params[1], sum2s[j]))
      elif not hasplt and all0:
         print('%-20s: pKa = %8s n = %8s RSS = %10s' % (ress[j], 'Inf', 'N/A',
                  'N/A'))
   
   if not hasplt:
      sys.exit(0)

   if opt.title:
      opt.title = opt.title.strip() + ' '

   # Plot me
   fig = plt.figure(1)
   ncols = min(data.shape[0] - 1, 3)
   nrows = max(1, int(np.ceil(data.shape[0] / 3)))
   for j in range(1, data.shape[0]):
      ax = fig.add_subplot(nrows, ncols, j)
      ax.axis([np.min(data[0])-1, np.max(data[0])+1, 0, 1])
      ax.set_xlabel('pH')
      ax.set_ylabel('Fraction Deprotonated')
      ax.set_title('%s Titration Curve' % ress[j])
      fcn_x = np.linspace(np.min(data[0])-1, np.max(data[0])+1, 1000)
      fcn_y = [f(i, paras[j][0], paras[j][1]) for i in fcn_x]
      ax.text(np.min(data[0])-0.9, 0.4, 'pKa = %.2f\nn    = %.2f\nRSS = %.4e' %
               (paras[j][0], paras[j][1], sum2s[j]))
      ax.grid(True)
      ax.plot(data[0], data[j], 'bo', fcn_x, fcn_y, 'r-')
   plt.show()
