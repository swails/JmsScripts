#!/usr/bin/env python

from optparse import OptionParser
from remd import pHRemLog
import sys, math, os

debug_printlevel = 1

def excepthook(exception_type, exception_value, tb):
   """ Replaces sys.excepthook so fatal exceptions kill all MPI threads and
       we can control the printing of tracebacks. Those are helpful for
       debugging purposes, but may be unsightly to users. debug_printlevel
       set above controls this behavior
   """
   import traceback
   if debug_printlevel > 0: traceback.print_tb(tb)
   sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
   sys.stderr.write('Exiting. All files have been retained.' + os.linesep)

sys.excepthook = excepthook

epilog = '''This program will calculate what would be the success ratio between
every pair of pH values in a given pH REM simulation'''

parser = OptionParser(epilog=epilog)
parser.add_option('-l', '--remlog', dest='input', metavar='FILE', default=None,
                  help='Input rem.log file')
parser.add_option('-d', '--debug', dest='debug', metavar='INT', default=0,
                  help='Set to 1 to turn on debugging (traceback) outputs')

opt, args = parser.parse_args()

if not opt.input or args:
   parser.print_help()
   sys.exit(1)

debug_printlevel = opt.debug

remlog = pHRemLog(opt.input)

numreps = len(remlog.reps)
numexchg = remlog.numexchg
LN_TO_LOG = math.log(10)

pH_vals = remlog.values[:]

# probs is the array of probability arrays. Each element is an array that stores
# the probability of exchange success with every replica that has a higher pH
# (i.e. we are assuming a symmetric matrix)
probs = [[] for i in range(numreps-1)]
for i in range(numreps-1): probs[i] = [0 for j in range(numreps-i-1)]

# Now let's go through each step and calculate the probability of exchange of
# every pH with every other pH
for i in range(numexchg):
   indexes = [rep.index[i] for rep in remlog.reps]
   prot_cnts = [rep.prot_cnt[i] for rep in remlog.reps]
   old_phs = [rep.old_pH[i] for rep in remlog.reps]
   new_phs = [rep.new_pH[i] for rep in remlog.reps]
#  print 'Indexes are: ', indexes
#  print prot_cnts
#  print old_phs
#  print new_phs
   for j in range(numreps):
      idx1 = indexes[j]
      for k in range(numreps):
         idx2 = indexes[k]
         if idx2 <= idx1: continue
#        print 'idx1,2 = ', idx1, idx2
#        print 'delta pH, delta prot = ', pH_vals[idx1]-pH_vals[idx2], prot_cnts[j]-prot_cnts[k]
         # Exchange probability is exp(DELTA prot_cnt * DELTA pH * LN_TO_LOG)
         delta = -(pH_vals[idx1]-pH_vals[idx2])*(prot_cnts[j]-prot_cnts[k]) * LN_TO_LOG
#        print idx1, idx2, len(probs), idx2-idx1-1, len(probs[idx1])
         probs[idx1][idx2-idx1-1] += min(1.0,math.exp(-delta))
#        print 'changing ',idx1, ', ',idx2-idx1-1, ' by (', delta,') ', min(1.0,math.exp(-delta))
#  print 'Probs are: ', probs


# Now we can report the averages:

print "   pH 1 |   pH 2 |   Avg. Acceptance Probability"
print "------------------------------------------------"
for i in range(numreps-1):
   for j in range(i+1, numreps):
      print "%7.2f |%7.2f |%30.6f" % (pH_vals[i], pH_vals[j], 
                                      probs[i][j-i-1]/numexchg)
