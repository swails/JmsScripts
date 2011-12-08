#!/usr/bin/env python

import math, sys, utilities, time, re
from optparse import OptionParser
from os import path

ttotstart = time.time()

parser = OptionParser()
parser.add_option('-f', '--file', dest='data_file', help='Input data file',
                  default=None)
parser.add_option('-o', '--output', dest='output_file', default=None,
                  help='Binned output file')
parser.add_option('-b', '--bins', dest='bins', default=None,
                  help='Number of bins [NxM]')
parser.add_option('-x', '--xbinrange', dest='xbinrange', default=None,
                  help='Range for bins in the x-direction')
parser.add_option('-y', '--ybinrange', dest='ybinrange', default=None,
                  help='Range for bins in the y-direction')
parser.add_option('-n', '--normalize', dest='normalize', default=False,
                  action='store_true', help='Normalize the histograms')
parser.add_option('-g', '--gnuplot', dest='gnuplot', default=None,
                  help='Name of gnuplot script to write.')
parser.add_option('-c', '--columns', dest='columns', default=None,
                  help='Which two columns to analyze [##,##] starting from 1')
parser.add_option('-d', '--delimiter', dest='delimiter', default=' ',
                  help='The string/character that delimits fields in the ' +
                  'input data file')
opt, args = parser.parse_args()

# Set up regular expressions to parse the input commands.
numbins_re = re.compile(r'(\d+)x(\d+)')
binrange_re = re.compile(r'([+-]{0,1}\d+.*\d*)-([+-]{0,1}\d+.*\d*)')
columns_re = re.compile(r'(\d+),(\d+)')

# Process the command-line flags
if not opt.output_file:
   parser.print_help()
   sys.exit(0)

if opt.bins:
   optmatch = numbins_re.match(opt.bins)
   if not optmatch:
      print 'Bad format for -b/--bins: INTxINT'
      sys.exit(1)
   bins = [int(i) for i in optmatch.groups()]
else:
   bins = [0,0]

if opt.xbinrange:
   optmatch = binrange_re.match(opt.xbinrange)
   if not optmatch:
      print 'Bad format for -x/--xbinrange: FLOAT-FLOAT'
      sys.exit(1)
   binrange = [float(x) for x in optmatch.groups()]
else:
   binrange = [0,0]

if opt.ybinrange:
   optmatch = binrange_re.match(opt.ybinrange)
   if not optmatch:
      print 'Bad format for -y/--ybinrange: FLOAT-FLOAT'
      sys.exit(1)
   binrange.extend([float(x) for x in optmatch.groups()])
else:
   binrange.extend([0,0])

if opt.columns:
   optmatch = columns_re.match(opt.columns)
   col1, col2 = tuple([int(i) for i in optmatch.groups()])
   # Adjust for indexing
   col1, col2 = col1 - 1, col2 - 1
else:
   col1 = 0
   col2 = 1

maxcol = max(col1, col2) + 1

if not opt.data_file:
   input_data = sys.stdin
else:
   if not path.exists(opt.data_file):
      print 'Error: Could not find %s!' % opt.data_file
      sys.exit(1)
   input_data = open(opt.data_file)

phipsibins = []   # the list of all bins in both directions
data1 = []        # array of 2-element arrays that contains every pair of points
data2 = []        # array of 2-element arrays that contains every pair of points
discarded = 0     # count number of points that have been discarded
pointweight = 1.0 # how much a single point is worth. 1 if not normalized.
script_name = opt.gnuplot

# load data into file
for line in input_data:
   words = line.split(opt.delimiter)
   if len(words) < maxcol: # skip any lines that don't have enough values
      continue
   try: # skip any lines where at least 1 of first 2 columns are not floats
      data1.append(float(words[col1]))
      data2.append(float(words[col2]))
   except ValueError:
      continue

input_data.close()

# Set up default ranges
if (binrange[0] == 0 and binrange[1] == 0) or \
            (binrange[2] == 0 and binrange[3] == 0):
   xmaxminholder = utilities.minmax(data1)
   ymaxminholder = utilities.minmax(data2)
   binrange[0] = math.floor(xmaxminholder[0])
   binrange[1] = math.ceil(xmaxminholder[1])
   binrange[2] = math.floor(ymaxminholder[0])
   binrange[3] = math.ceil(ymaxminholder[1])

# Set up default number of bins according to "Scott's Choice"
if bins[0] == 0 or bins[1] == 0:
   xinttmp = 3.5 * utilities.stdev(data1,'no') / float(len(data1)) ** (1/3)
   yinttmp = 3.5 * utilities.stdev(data2,'no') / float(len(data2)) ** (1/3)
   bins[0] = int(math.ceil(binrange[1] - binrange[0]/xinttmp))
   bins[1] = int(math.ceil(binrange[3] - binrange[2]/yinttmp))
if opt.normalize: 
   pointweight /= float(len(data1)) * (
            binrange[1]-binrange[0])*(binrange[3]-binrange[2])/(bins[0]*bins[1])

xinterval = (binrange[1] - binrange[0])/bins[0]
yinterval = (binrange[3] - binrange[2])/bins[1]

# create a large 1-D array with every bin 
# (x1y1 x1y2 ... x1yN x2y1 x2y2 ... ... xNyN)
for x in range(bins[0]):
   for y in range(bins[1]):
      phipsibins.append(0)

for x in range(len(data1)):

   xval = data1[x]
   yval = data2[x]

   if xval > binrange[1] or xval < binrange[0] or yval > binrange[3] or \
               yval < binrange[2]:
      discarded += 1
      continue

   xval -= binrange[0]
   yval -= binrange[2]

   binnum = int(math.floor(xval/xinterval) * bins[1] + math.floor(yval/yinterval))

   phipsibins[binnum] += pointweight

xprintval = binrange[0]
yprintval = binrange[2]
outputfile = open(opt.output_file,'w')
for x in range(len(phipsibins)):
   if x != 0 and x % bins[1] == 0:
      xprintval += xinterval
      yprintval = binrange[2]
      outputfile.write('\n')
   outputfile.write('%s %s %s\n' % (xprintval,yprintval,phipsibins[x]))
   yprintval += yinterval
   
outputfile.close()

if script_name:
   script = open(script_name,'w')
   script.write("""unset surface
set contour base
set cntrparam levels 20
set cntrparam bspline
set cntrparam order 7
set view 0,0
unset ztics\n""")
   script.write("splot '%s' w l\n" % opt.output_file)
   script.close()


ttotend = time.time()

print 'Binning Results:\n'
print 'Time Taken:      %.4f sec.' % (ttotend - ttotstart)
print 'Data file:       ' + opt.data_file
print 'Output file:     ' + opt.output_file
if script_name:
   print 'GNUPLOT script:  ' + script_name
print 'Numbers of Bins: %d x %d' % (bins[0], bins[1])
print 'Bin X-Range:     %.4f - %.4f' % (binrange[0], binrange[1])
print 'Bin Y-Range:     %.4f - %.4f' % (binrange[2], binrange[3])
print 'Points Omitted:  %d' % discarded
print 'Done!'
