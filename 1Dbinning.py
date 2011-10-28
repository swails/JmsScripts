#!/usr/bin/env python

from optparse import OptionParser
import math, sys, utilities, time

ttotstart = time.time()

clopts = OptionParser()
clopts.add_option('-f', '--file', dest="data_file", 
                  help="Input data file", default='none')
clopts.add_option('-o', '--output', dest="output_file",
                  help="Binned output file", default=None)
clopts.add_option('-b', '--bins', dest="bins", help="Number of bins to use",
                  default=0, type="int")
clopts.add_option('-r', '--binrange', dest="binrange",
                  help="Range of bins to use")
clopts.add_option('-n', '--normalize', dest="normalize", action="store_true",
                  default=False, help="Normalize the histograms")
clopts.add_option('-g', '--gnuplot', dest="script_name", default='',
                  help="Gnuplot script name")
clopts.add_option('-d', '--delimiter', dest="delimiter", default=None,
               help="The column delimiter (defaults to any kind of whitespace)")
clopts.add_option('-c', '--column', dest='column', default=1, type="int",
                  help="Which column to pull the data from")
(opts, args) = clopts.parse_args()

if not opts.output_file:
   clopts.print_help()
   sys.exit(1)

# set variables with default values
bins = opts.bins 
if opts.binrange:
   binrange = [int(opts.binrange.split('-')[0].strip()),
               int(opts.binrange.split('-')[1].strip())]
else:
   binrange = [0,0]
data_file = opts.data_file
output_file = opts.output_file
normalize = opts.normalize
script_name = opts.script_name
delimiter = opts.delimiter
column = opts.column
phipsibins = []                    # the list of all bins in both directions (dimension bins)
data = []                          # array of 2-element arrays that contains every pair of points
discarded = 0                      # count number of points that have been discarded
pointweight = 1.0                  # how much a single point is worth. 1 if not normalized.

# read in data, check for existing file
if data_file == 'none':
   input_data = sys.stdin
else:
   try:
      input_data = open(data_file,'r')
   except IOError:
      print 'Error: data file ' + data_file + ' not found!'
      sys.exit(1)

# load data into file
for line in input_data:
   if delimiter: words = line.split(delimiter)
   else: words = line.split()
   if len(words) < column: # skip any lines that don't have enough values
      continue
   try: # skip any lines where value is not a float
      data.append(float(words[column-1]))
   except ValueError:
      continue

input_data.close()

# Set up default ranges
if (binrange[0] == 0 and binrange[1] == 0):
   xmaxminholder = utilities.minmax(data)
   binrange[0] = math.floor(xmaxminholder[0])
   binrange[1] = math.ceil(xmaxminholder[1])

# Set up default number of bins according to "Scott's Choice"
if bins == 0:
   inttmp = 3.5 * utilities.stdev(data,'no') / float(len(data)) ** (1/3)
   bins = int(math.ceil(binrange[1] - binrange[0]/inttmp))

if normalize: 
   pointweight /= float(len(data)) * (binrange[1] - binrange[0])/bins

interval = (binrange[1] - binrange[0])/bins

# create a large 1-D array with every bin (x1y1 x1y2 ... x1yN x2y1 x2y2 ... ... xNyN)
for x in range(bins):
   phipsibins.append(0)

for x in range(len(data)):

   xval = data[x]

   if xval > binrange[1] or xval < binrange[0]:
      discarded += 1
      continue

   xval -= binrange[0]

   binnum = int(math.floor(xval/interval))

   try:
      phipsibins[binnum] += pointweight
   except:
      phipsibins[binnum-1] += pointweight

xprintval = binrange[0]
outputfile = open(output_file,'w')
for x in range(len(phipsibins)):
   outputfile.write(str(xprintval) + ' ' + str(phipsibins[x]) + '\n')
   xprintval += interval
   
outputfile.close()

if script_name != '':
   script = open(script_name,'w')
   script.write("plot '{0}' w l".format(output_file))
   script.close()


ttotend = time.time()

print 'Binning Results:\n'
print 'Time Taken:      {0:.3f} sec.'.format(ttotend - ttotstart)
print 'Data file:       ' + data_file
print 'Output file:     ' + output_file
if script_name != '':
   print 'GNUPLOT script:  ' + script_name
print 'Numbers of Bins: ' + str(bins)
print 'Bin X-Range:     ' + str(binrange[0]) + ' - ' + str(binrange[1])
print 'Points Omitted:  ' + str(discarded)
print 'Done!'
