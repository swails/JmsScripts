#! /apps/python/26/bin/python

import math, sys, utilities, time

ttotstart = time.time()

# function to print the usage statement
def printusage():
   import sys
   print 'Usage: 1Dbinning.py -f data_file        \\'
   print '                    -o output_file      \\'
   print '                   {-bins NUM}          \\'
   print '                   {-binrange NUM-NUM}  \\'
   print '                   {-normalize}         \\'
   print '                   {-gnuplot script_name}'
   sys.exit()

# answer any call for help
if len(sys.argv) == 1 or sys.argv[1].startswith('-h') or sys.argv[1].startswith('--h'):
   printusage()

# set variables with default values
bins = 0                       # number of bins in both dimensions
binrange = [0,0]               # range for each set of bins
data_file = 'none'             # file with the initial data, must be 1 column
output_file = 'binned_output.dat'  # output file to contain the binned data plottable by gnuplot
phipsibins = []                    # the list of all bins in both directions (dimension bins)
data = []                         # array of 2-element arrays that contains every pair of points
discarded = 0                      # count number of points that have been discarded
normalize = False                  # whether or not to normalize
pointweight = 1.0                  # how much a single point is worth. 1 if not normalized.
script_name = ''

# My traditional parser
try:
   for x in range(len(sys.argv)):
      if sys.argv[x] == '-bins':
         bins = int(sys.argv[x+1].split('x')[0].strip())
      elif sys.argv[x] == '-binrange':
         binrange[0] = float(sys.argv[x+1].split('-')[0].strip())
         binrange[1] = float(sys.argv[x+1].split('-')[1].strip())
      elif sys.argv[x] == '-f':
         data_file = sys.argv[x+1]
      elif sys.argv[x] == '-o':
         output_file = sys.argv[x+1]
      elif sys.argv[x] == '-normalize':
         normalize = True
      elif sys.argv[x] == '-gnuplot':
         script_name = sys.argv[x+1]
      elif sys.argv[x].startswith('-'):
         print 'Error: Unknown flag ' + sys.argv[x]
         printusage()
except IndexError:
   print 'Error: Command line error!'
   printusage()
except ValueError:
   print 'Error: Bins must be integers and binrange must be floats!'
   printusage()

# read in data, check for existing file
if data_file == 'none':
   input_data = sys.stdin
else:
   try:
      input_data = open(data_file,'r')
   except IOError:
      print 'Error: data file ' + data_file + ' not found!'
      printusage()

# load data into file
for line in input_data:
   words = line.split()
   if len(words) < 1: # skip any lines that don't have enough values
      continue
   try: # skip any lines where value is not a float
      data.append(float(words[0]))
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
print 'Time Taken:      {0:.3f} min.'.format((ttotend - ttotstart) / 60)
print 'Data file:       ' + data_file
print 'Output file:     ' + output_file
if script_name != '':
   print 'GNUPLOT script:  ' + script_name
print 'Numbers of Bins: ' + str(bins)
print 'Bin X-Range:     ' + str(binrange[0]) + ' - ' + str(binrange[1])
print 'Points Omitted:  ' + str(discarded)
print 'Done!'
