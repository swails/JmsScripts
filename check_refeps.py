#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
from remd import HRemLog

parser = OptionParser()
parser.add_option('-l', '--log', dest='input_log', metavar='FILE', default=None,
                  help='Input log file to parse')
opt, arg = parser.parse_args()

if opt.input_log is None:
   parser.print_help()
   quit(1)

remlog = HRemLog(opt.input_log)

running_left_fes = np.zeros(len(remlog.reps))
running_right_fes = np.zeros(len(remlog.reps))
left_fes = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
right_fes = [np.zeros(remlog.numexchg) for i in range(len(remlog.reps))]
n_r_exch = [0 for i in range(len(remlog.reps))]
n_l_exch = [0 for i in range(len(remlog.reps))]

#print '   PotEne 1   PotEne 2'
for i, rep in enumerate(remlog.reps):
   for j in range(remlog.numexchg):
      nrep = remlog.reps[rep.neighbor_index[j]]
      if (j % 2 == i % 2):
         n_l_exch[i] += 1
         running_left_fes[i] += np.exp(-(rep.potene1[j] - nrep.potene2[j])
                              * 503.01 / 300)
         left_fes[i][j] = -1 / 503.01 * 300 * np.log(running_left_fes[i] /
                                                     n_l_exch[i])
      else:
         n_r_exch[i] += 1
         running_right_fes[i] += np.exp(-(rep.potene1[j] - nrep.potene2[j])
                               * 503.01 / 300)
         right_fes[i][j] = -1 / 503.01 * 300 * np.log(running_right_fes[i] /
                                                      n_r_exch[i])
#     if rep.index[j] == 0:
#        print '%10.3f %10.3f' % (rep.potene1[j], nrep.potene2[j])
#        print 'right: %d left: %d' % (n_r_exch[i], n_l_exch[i])
#        print 'Left FE: %.4f Right FE: %.4f' % (left_fes[i][j], right_fes[i][j])

for j in range(remlog.numexchg):
   for i in range(len(remlog.reps)):
      print 'Replica %4d: Left FE = %10.4f Right FE = %10.4f' % \
            (i, left_fes[i][j], right_fes[i][j])
   print '='*80
