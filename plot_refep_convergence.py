#!/usr/bin/env python

"""
This program plots convergence of REFEP calculations
"""

from __future__ import division
import numpy as np
import math
from remd import HRemLog
from argparse import ArgumentParser

if __name__ == '__main__':
   parser = ArgumentParser()
   group = parser.add_argument_group('Files')
   group.add_argument('-i', '--input', dest='remlog', metavar='FILE',
            default=[], help='Input REMD log file', required=True,
            nargs='*')
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
            dest='title', help='Title of the plot. Default [[ %(default)s ]]')
   group.add_argument('-d', '--plot-diff', dest='diff', default=False,
            action='store_true', help='''Plot the difference between the forward
            and reverse free energies to show convergence.''')
   group.add_argument('-l', '--lag', default=0, type=int, metavar='INT',
            help='''The number of points to omit from plotting the differences.
            Default is %(default)s''')
   group.add_argument('--data-file', dest='datafile', metavar='FILE',
            default=None, help='Save the data to a file instead.')
   
   opt = parser.parse_args()

   if opt.datafile is None:
      import matplotlib
      if opt.output is not None:
         matplotlib.use('Agg')
      import matplotlib.pyplot as plt

   remlogs = [HRemLog(remlog) for remlog in opt.remlog]

   nreps = len(remlogs[0].reps)
   tot_exchg = sum([remlog.numexchg for remlog in remlogs])
   time = np.arange(0, tot_exchg*opt.eaf, opt.eaf*2)

   KBT = 0.00199 * 300
   ONEKT = 1 / KBT
   # Now compute the left- and right- free energies
   rdiffs = [[] for i in range(nreps)]
   ldiffs = [[] for i in range(nreps)]
   for remlog in remlogs:
      for i, rep in enumerate(remlog.reps):
         for j, oidx in enumerate(rep.neighbor_index):
            orep = remlog.reps[oidx]
            if rep.index[j] == oidx+1 or (rep.index[j] == 0 and oidx == nreps-1):
               ldiffs[i].append(math.exp((rep.potene1[j]-orep.potene2[j])*ONEKT))
            else:
               rdiffs[i].append(math.exp((rep.potene1[j]-orep.potene2[j])*ONEKT))

   left_fe = np.zeros(len(ldiffs[0]))
   right_fe = np.zeros(len(rdiffs[0]))
   for i in range(len(ldiffs[0])):
      for j in range(1, len(ldiffs)):
         left_fe[i] -= KBT * math.log(sum(ldiffs[j][:i+1]) / (i+1))
   for i in range(len(rdiffs[0])):
      for j in range(len(rdiffs)-1):
         right_fe[i] += KBT * math.log(sum(rdiffs[j][:i+1]) / (i+1))

   middle = (left_fe[i] + right_fe[i]) * 0.5
   print 'The final, average free energy is %.4f' % middle
   if opt.datafile is not None:
      with open(opt.datafile, 'w') as f:
         f.write('%10s %15s %15s\n' % ('Frame', 'Left FE', 'Right FE'))
         f.write('-'*(42))
         f.write('\n')
         for i in range(min(len(left_fe), len(right_fe))):
            f.write('%10d %15g %15g\n' % i, left_fe[i], right_fe[i])
      sys.exit(0)

   fig = plt.figure(1, figsize=(8,5))
   ax = fig.add_subplot(111)
   low = min(np.min(left_fe), np.min(right_fe))
   high = max(np.max(left_fe), np.max(right_fe))
   tdiff = max(middle - low, high - middle)
   ax.set_ylim(middle-tdiff, middle+tdiff)

   fontdict = dict(size=20, family='sans-serif')
   fontdict2 = dict(size=18, family='sans-serif')
   ax.set_title(opt.title, fontdict=fontdict)
   ax.set_xlabel('Time (ps)', fontdict=fontdict2)
   ax.set_ylabel('Free Energy (kcal mol$^{-1}$)', fontdict=fontdict2)
   pl1, = ax.plot(time, left_fe, 'k-', lw=2)
   pl2, = ax.plot(time, right_fe, 'r-', lw=2)
   ax.grid(lw=1)

   # Now plot the difference if requested
   if opt.diff:
      ax2 = ax.twinx()
      pl3, = ax2.plot(time[opt.lag:], (left_fe-right_fe)[opt.lag:], 'g--', lw=2)
      ax2.set_ylabel('Difference (kcal mol$^{-1}$)', color='g',
                     fontdict=fontdict2)
      ax.legend((pl1, pl2, pl3), ('Forward FEP', '-Reverse FEP', 'Difference'),
                loc=1)
      # Change the color of the tics and its labels
      for tic in ax2.get_yticklabels(): tic.set_color('g')
   else:
      ax.legend((pl1, pl2), ('Forward FEP', '-Reverse FEP'), loc=1)
   if opt.output:
      fig.savefig(opt.output)
      print ('Saved %s. Done!' % opt.output)
   else:
      plt.show()
