#!/usr/bin/env python
from __future__ import division

from mdoutremd import HRemMdout
import numpy as np
from optparse import OptionParser
from remd import HRemLog
import sys

usage='%prog [Options] | [mdout1] [mdout2] [mdout3] ... [mdoutN]'
parser = OptionParser(usage=usage)
parser.add_option('-l', '--log', dest='input_log', metavar='FILE', default=None,
                  help='Input log file to parse')
opt, arg = parser.parse_args()

if opt.input_log is None:
   if len(arg) == 0:
      parser.print_help()
      quit(1)
   remlog = HRemMdout(arg)
else:
   remlog = HRemLog(opt.input_log)

running_left_fes = np.zeros(len(remlog.reps))
running_right_fes = np.zeros(len(remlog.reps))
left_fes = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
right_fes = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
left_runs = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
right_runs = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
left_nexch = [[0 for i in range(remlog.numexchg)] for i in range(len(remlog.reps))]
right_nexch = [[0 for i in range(remlog.numexchg)] for i in range(len(remlog.reps))]
n_r_exch = [0 for i in range(len(remlog.reps))]
n_l_exch = [0 for i in range(len(remlog.reps))]
#fls = [open('python.00%d' % d, 'w') for d in range(len(remlog.reps))]

#print '   PotEne 1   PotEne 2'
for i, rep in enumerate(remlog.reps):
   for j in range(remlog.numexchg):
      nrep = remlog.reps[rep.neighbor_index[j]]
#     fls[i].write('='*80 + '\n')
      if (j % 2 != i % 2):
         n_l_exch[i] += 1
         running_left_fes[i] += np.exp((rep.potene1[j] - nrep.potene2[j])
                              * 503.01 / 300)
         left_fes[i][j] = 1 / 503.01 * 300 * np.log(running_left_fes[i] /
                                                     n_l_exch[i])
         left_runs[i][j] = running_left_fes[i]
         left_nexch[i][j] = n_l_exch[i]
         try:
            right_fes[i][j] = right_fes[i][j-1]
         except IndexError:
            pass
      else:
         n_r_exch[i] += 1
         running_right_fes[i] += np.exp((rep.potene1[j] - nrep.potene2[j])
                               * 503.01 / 300)
         right_fes[i][j] = 1 / 503.01 * 300 * np.log(running_right_fes[i] /
                                                      n_r_exch[i])
         right_runs[i][j] = running_right_fes[i]
         right_nexch[i][j] = n_r_exch[i]
         try:
            left_fes[i][j] = left_fes[i][j-1]
         except IndexError:
            pass

for j in range(remlog.numexchg):
   for i in range(len(remlog.reps)):
      print ('Replica %4d: Left FE = %10.4f (%10.4f) Right FE = %10.4f (%10.4f)'
             % (i, left_fes[i][j], remlog.reps[i].left_fe[j], right_fes[i][j],
                remlog.reps[i].right_fe[j])
            )
#  print 'Group 1: num_right_exchg = %10d %10d' % (right_nexch[i][0],
#                                                  right_nexch[i][1])
#  print 'Group 1: num_left_exchg  = %10d %10d' % (left_nexch[i][0],
#                                                  left_nexch[i][1])
#  print 'Group 1: total_right_fe  = %s %s' % (right_runs[i][0],
#                                              right_runs[i][1])
#  print 'Group 1: total_left_fe   = %s %s' % (left_runs[i][0],
#                                              left_runs[i][1])
   print '='*80
