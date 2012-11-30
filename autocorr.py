#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
import sys

parser = OptionParser()
parser.add_option('-f', '--file', dest="data_file", 
                  help="Input data file", default=None)
parser.add_option('-o', '--output', dest="output_file",
                  help="Autocorrelation output file", default=None)
parser.add_option('-n', '--normalize', dest="normalize", action="store_true",
                  default=False, help="Normalize the histograms")
parser.add_option('-d', '--delimiter', dest="delimiter", default=None,
               help="The column delimiter (defaults to any kind of whitespace)")
parser.add_option('-c', '--column', dest='column', default=1, type="int",
                  help="Which column to pull the data from")

opt, arg = parser.parse_args()

# read in data, check for existing file
if opt.data_file is None:
   input_data = sys.stdin
else:
   try:
      input_data = open(opt.data_file,'r')
   except IOError:
      print 'Error: data file ' + opt.data_file + ' not found!'
      sys.exit(1)

if opt.output_file is None:
   outfile = sys.stdout
else:
   outfile = open(opt.output_file, 'w')

if opt.delimiter is None:
   splitmethod = lambda x: str.split(x)
else:
   splitmethod = lambda x: str.split(x, opt.delimiter)

data = []
for line in input_data:
   # Skip comments
   if line.startswith('#'): continue
   # Skip over non-qualifying lines
   try:
      data.append(float(splitmethod(line)[opt.column-1].strip()))
   except (ValueError, IndexError):
      pass

# Now we have our data, convert to ndarray
data = np.asarray(data)

# If we normalize it then subtract off the mean, divide by the standard
# deviation, and divide by sqrt of the length
if opt.normalize:
   data -= data.mean()
   data /= data.std()
   data2 = data.copy() / data.shape[0]
   acor = np.correlate(data, data2, 'full')
else:
   acor = np.correlate(data, data2, 'full')

# Only keep the positive lags
acor = acor[len(acor)//2:]

# Dump the data
for i, val in enumerate(acor):
   outfile.write('%12s %15s\n' % (i, val))
