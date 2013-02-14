#!/usr/bin/env python

"""
This program will read in a pH stat file divided up into multiple chunks
and it outputs the standard deviation of the fraction deprotonated
"""

from __future__ import division

import sys, math
from optparse import OptionParser

from running_titration_curve import pHStatFile

epilog = '''
This program will read in a pH stat file divided up into multiple chunks
and it outputs the standard deviation of the fraction deprotonated
'''

parser = OptionParser(epilog=epilog)

parser.add_option('-i', '--input', dest='input_file', metavar='FILE',
                  help='Input file with the pH data. Default standard input',
                  default=None)
parser.add_option('-p', '--frac-prot', dest='frac_prot', default=False,
                  action='store_true', help='Analyze fraction protonated ' +
                  'instead of fraction deprotonated? NOT default behavior')

opt, arg = parser.parse_args()

if opt.input_file is None:
   infile = sys.stdin
else:
   try: 
      infile = open(opt.input_file, 'r')
   except IOError:
      sys.stderr.write("Error opening %s!\n" % opt.input_file)
      sys.exit(1)

statfile = pHStatFile(infile)
statfile.my_re = statfile.chunkre

# if we want fraction protonated, we need to do nothing
def frac_prot(val): return val
# if we want fraction deprotonated, we need to subtract from 1
def frac_deprot(val): return 1-val

# Calculating our stdev
def stdev(sum1, sum2, num_vals):
   # sum1 is the running sum of a sequence
   # sum2 is the running sum squared of that sequence
   # num_vals is the number of numbers in that sequence
   return math.sqrt(abs(sum2 / num_vals - (sum1 / num_vals)**2))

# Alias function -- avoids conditional eval in ever iteration of a tight loop
if opt.frac_prot:
   transform = frac_prot
else:
   transform = frac_deprot

sum_vals = {}    # Running sum of fraction (de)protonated
sum2_vals = {}   # Running sum of fraction (de)protonated squared (for stddev)
num_samples = 0

for key in statfile.list_of_residues:
   sum_vals[key] = 0.0
   sum2_vals[key] = 0.0
nexres = statfile.get_next_residue()
print nexres
printed = True
while nexres:
   key = '_'.join([nexres[0], str(nexres[1])])
   sum_vals[key] += transform(nexres[4])
   sum2_vals[key] += transform(nexres[4])**2
   if statfile.list_of_residues[0] == key: 
      num_samples += 1
   nexres = statfile.get_next_residue()
   if printed:
      print nexres
      printed = False

# We're all done... Now calculate all stdevs.
print 'I found %d total points' % num_samples
for resid in statfile.list_of_residues:
   print '%10s: %f' % (resid, 
                       stdev(sum_vals[resid], sum2_vals[resid], num_samples))
