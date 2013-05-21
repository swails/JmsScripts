#!/usr/bin/env python

"""
This program plots convergence of REFEP calculations
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from remd import HRemLog
from argparse import ArgumentParser

if __name__ == '__main__':
   parser = ArgumentParser()
   group = parser.add_argument_group('Files')
   group.add_argument('-i', '--input', dest='remlog', metavar='FILE',
            default=None, help='Input REMD log file', required=True)
   group.add_argument('-o', '--output', dest='output', metavar='IMAGE_FILE',
            default=None, help='''Image with graphed convergence. By default,
            just show the graph on the screen.''')
   group = parser.add_argument_group('Simulation Properties', '''These are
            properties of the simulation you ran that are necessary to set up
            labels and axes on the plot.''')
   group.add_argument('-e', '--eaf', '--exchange-attempt-frequency', dest='eaf',
            type=float, help='Exchange attempt frequency in ps^-1',
            required=True)
   group = parser.add_argument_group('Plot Options', '''These options change how
            the plot appears.''')
   group.add_argument('-t', '--title', default='REFEP Free Energies',
            dest='title', help='Title of the plot.')
   
   opt = parser.parse_args()

   remlog = HRemLog(opt.remlog)

   left_fe = np.zeros(remlog.numexchg)
   right_fe = np.zeros(remlog.numexchg)
   time = np.arange(0, remlog.numexchg*opt.eaf, opt.eaf)

   for i in range(remlog.numexchg):
      lfe = 0.0
      rfe = 0.0
      for j, rep in enumerate(remlog.reps):
         if j == 0:
            rfe += rep.right_fe[i]
         elif j == len(remlog.reps) - 1:
            lfe += rep.left_fe[i]
         else:
            rfe += rep.right_fe[i]
            lfe += rep.left_fe[i]
      left_fe[i] = -lfe
      right_fe[i] = rfe

   fig = plt.figure(1, figsize=(8,5))
   ax = fig.add_subplot(111)

   fontdict = dict(size=20, family='sans-serif')
   fontdict2 = dict(size=18, family='sans-serif')
   ax.set_title(opt.title, fontdict=fontdict)
   ax.set_xlabel('Time (ps)', fontdict=fontdict2)
   ax.set_ylabel('Free Energy (kcal mol$^{-1}$)', fontdict=fontdict2)
   pl1, = ax.plot(time, left_fe, 'k-', lw=2)
   pl2, = ax.plot(time, right_fe, 'r-', lw=2)
   ax.grid(lw=1)
   ax.legend((pl1, pl2), ('Forward FEP', '-Reverse FEP'), loc='best')

   if opt.output:
      fig.savefig(opt.output)
      print ('Saved %s. Done!' % opt.output)
   else:
      plt.show()
