#! /usr/bin/env python

"""
This program effectively works like 'awk' and prints out the given columns from
a text file with a given delimiter
"""

from optparse import OptionParser
import sys

parser = OptionParser(usage='%prog [Options] column1 [column2 [column3 [...] ] ]',
                      epilog="This program effectively works like 'awk' and " +
                      "prints out the given columns from a text file with a " +
                      "given delimiter")
parser.add_option('-d', '--delim', dest='delimiter', default=' ',
                  help='Delimiter to split strings on. Default is a space')
parser.add_option('-i', '--input', dest='input_file', default=None,
                  help='Input file to get data from. Defaults from stdin')
parser.add_option('-o', '--output-delimiter', dest='out_delim', default=' ',
                  help='Delimiter to separate variables in the output. ' +
                  'Default is a space')
parser.add_option('-n', '--numbers-only', dest='numbers_only', default=False,
                  action='store_true', help='Only print out numbers')
parser.add_option('-a', '--all', dest='numbers_only', action='store_false',
                  help='Print out numbers and words (Default behavior)')
opt, args = parser.parse_args()

cols = [int(i) for i in args]

if opt.input_file:
   infile = open(opt.input_file, 'r')
else:
   infile = sys.stdin

if not args:
   # In this case, open up the file, print out the first 5 lines of text,
   # and separate the columns
   lines = infile.readlines()
   for i in range(max(len(lines)-1,5)):
      print ' || '.join(lines[i].strip().split(opt.delimiter))

else:
   for line in infile:
      words = line.split(opt.delimiter)
      try:
         out_words = [words[i-1] for i in cols]
      except IndexError:
         continue
      if opt.numbers_only:
         try:
            test_words = [float(i) for i in out_words]
         except ValueError:
            continue
      print opt.out_delim.join(out_words)
