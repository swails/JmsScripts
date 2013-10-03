#!/usr/bin/env python
from __future__ import division, print_function

import os, sys
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-f', '--input-file', dest='inputs', metavar='FILE',
                    nargs='*', default=None, help='''Input files. Default from
                    stdin''')
parser.add_argument('-c', '--column', dest='column', metavar='INT', type=int,
                    help='Column number (1st column is `1\', the default)',
                    default=1)

opt = parser.parse_args()

if opt.inputs is None:
   files = [sys.stdin]
else:
   files = [open(f, 'r') for f in opt.inputs]

col = opt.column - 1

for file in files:
   for line in file:
      try:
         num = float(line.split()[col])
      except ValueError:
         continue
      try:
         _max = max(_max, num)
      except NameError:
         _max = num
      try:
         _min = min(_min, num)
      except NameError:
         _min = num

prog = os.path.basename(sys.argv[0]).lower()
if prog.startswith('max'):
   print(_max)
elif prog.startswith('min'):
   print(_min)
else:
   print('%g\t%g' % (_min, _max))
