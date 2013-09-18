#!/usr/bin/env python

# Finds the sum of a specific column of a given file

def do_sum(f, c):
   mysum = 0
   for line in f:
      try:
         mysum += float(line.split()[c])
      except (ValueError, IndexError):
         pass
   return mysum

import sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--input', dest='input', metavar='FILE', default=[],
                  help='Input file. Default is read from stdin', nargs='*')
parser.add_argument('-n' ,'--column', dest='col', metavar='INT', default=1,
                  type=int, help='Column with data you want to sum.')
opt = parser.parse_args()

sum = 0
if not opt.input:
   sum += do_sum(sys.stdin, opt.col-1)
else:
   for fname in opt.input:
      with open(fname, 'r') as f:
         sum += do_sum(sys.stdin, opt.col-1)

print sum
