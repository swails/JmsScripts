#!/usr/bin/env python

# Finds the average of a specific column of a given file

from optparse import OptionParser, OptionGroup
from os.path import exists
import sys

parser = OptionParser(usage='%prog [Options] file1 [file2 [file3 [... ] ] ]',
                      epilog='Calculates some statistics for a given column ' +
                             'of data. No input file means read from STDIN')
parser.add_option('-c', '--column', dest='col', type='int', default=1,
                  help='Which column of data to analyze (Default %default)')
parser.add_option('-d', '--delimiter', dest='delim', default=' ',
                  help='What to consider a column delimiter. Default is a space')
group = OptionGroup(parser, 'Verbose options', 'These options control how ' +
                    'much is printed as output.')
group.add_option('-l', '--loud', dest='verbose', default=False,
                 action='store_true', help='Print helpful words and the ' +
                 'standard deviation. NOT default behavior')
group.add_option('-q', '--quiet', dest='verbose', action='store_false',
                 help='Just print the average as a single value. ' +
                 'Default behavior')
parser.add_option_group(group)

opt, args = parser.parse_args()

run_sum = 0
num_vals = 0
run_sum2 = 0

if args:
   exist_array = [exists(fname) for fname in args]
   if False in exist_array:
      print >> sys.stderr, 'Error: %s does not exist!' % (
               args[exist_array.index(False)])
      parser.print_help()
      sys.exit(1)
   
   for fname in args:
      datafile = open(fname, 'r')
      
      for line in datafile:
         words = line.split(opt.delim)
         try:
            val = float(words[opt.col - 1])
         except ValueError: continue
         except IndexError: continue
         else:
            run_sum += val
            run_sum2 += val * val
            num_vals += 1
      
      datafile.close()

else: # read from stdin

   for line in sys.stdin:
      words = line.split(opt.delim)
      try:
         val = float(words[opt.col - 1])
      except ValueError: continue
      except IndexError: continue
      else:
         run_sum += val
         run_sum2 += val
         num_vals += 1

# Compute avg, stdev
ave = run_sum / number_entries
ave2 = run_sum2 / number_entries
std = math.sqrt(ave2 - ave * ave)
      
if opt.verbose:
   print 'The average of column %d is: %15.5g' % (opt.col, ave)
   print 'The standard deviation of column %d is: %15.8g' % (opt.col, std)
else:
   print ave
