#!/usr/bin/env python

import math, sys, utilities, time

ttotstart = time.time()

# function to print the usage statement
def printusage():
   import sys
   print 'Usage: 2Dbinning.py -f data_file               \\'
   print '                    -o output_file             \\'
   print '                   {-bins NUMxNUM}             \\'
   print '                   {-binrange NUM-NUM NUM-NUM} \\'
   print '                   {-normalize}                \\'
   print '                   {-gnuplot script_name}'
   sys.exit()

# answer any call for help
if len(sys.argv) == 1 or sys.argv[1].startswith('-h') or sys.argv[1].startswith('--h'):
   printusage()

# set variables with default values
bins = [0,0]                       # number of bins in both dimensions
binrange = [0,0,0,0]               # range for each set of bins
data_file = 'dihedrals.dat'        # file with the initial data, must be 2 columns
output_file = 'binned_output.dat'  # output file to contain the binned data plottable by gnuplot
phipsibins = []                    # the list of all bins in both directions (dimension bins[0]xbins[1])
data1 = []                         # array of 2-element arrays that contains every pair of points
data2 = []                         # array of 2-element arrays that contains every pair of points
discarded = 0                      # count number of points that have been discarded
normalize = False                  # whether or not to normalize
pointweight = 1.0                  # how much a single point is worth. 1 if not normalized.
script_name = ''

# My traditional parser
try:
   for x in range(len(sys.argv)):
      if sys.argv[x] == '-bins':
         bins[0] = int(sys.argv[x+1].split('x')[0].strip())
         bins[1] = int(sys.argv[x+1].split('x')[1].strip())
      elif sys.argv[x] == '-binrange':
         binrange[0] = float(sys.argv[x+1].split('-')[0].strip())
         binrange[1] = float(sys.argv[x+1].split('-')[1].strip())
         binrange[2] = float(sys.argv[x+2].split('-')[0].strip())
         binrange[3] = float(sys.argv[x+2].split('-')[1].strip())
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
try:
   input_data = open(data_file,'r')
except IOError:
   print 'Error: data file ' + data_file + ' not found!'
   printusage()

# load data into file
for line in input_data:
   words = line.split()
   if len(words) < 2: # skip any lines that don't have enough values
      continue
   try: # skip any lines where at least 1 of first 2 columns are not floats
      data1.append(float(words[0]))
      data2.append(float(words[1]))
   except ValueError:
      continue

input_data.close()

# Set up default ranges
if (binrange[0] == 0 and binrange[1] == 0) or (binrange[2] == 0 and binrange[3] == 0):
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
   if normalize: pointweight /= float(len(data1))

xinterval = (binrange[1] - binrange[0])/bins[0]
yinterval = (binrange[3] - binrange[2])/bins[1]

# create a large 1-D array with every bin (x1y1 x1y2 ... x1yN x2y1 x2y2 ... ... xNyN)
for x in range(bins[0]):
   for y in range(bins[1]):
      phipsibins.append(0)

for x in range(len(data1)):

   xval = data1[x]
   yval = data2[x]

   if xval > binrange[1] or xval < binrange[0] or yval > binrange[3] or yval < binrange[2]:
      discarded += 1
      continue

   xval -= binrange[0]
   yval -= binrange[2]

   binnum = int(math.floor(xval/xinterval) * bins[1] + math.floor(yval/yinterval))

   phipsibins[binnum] += pointweight

xprintval = binrange[0]
yprintval = binrange[2]
outputfile = open(output_file,'w')
for x in range(len(phipsibins)):
   if x != 0 and x % bins[1] == 0:
      xprintval += xinterval
      yprintval = binrange[2]
      outputfile.write('\n')
   outputfile.write(str(xprintval) + ' ' + str(yprintval) + ' ' + str(phipsibins[x]) + '\n')
   yprintval += yinterval
   
outputfile.close()

if script_name != '':
   script = open(script_name,'w')
   script.write("""unset surface
set contour base
set cntrparam levels 20
set cntrparam bspline
set cntrparam order 7
set view 0,0
unset ztics\n""")
   script.write('splot \'' + output_file + '\' w l')
   script.close()


ttotend = time.time()

print 'Binning Results:\n'
print 'Time Taken:      {0:.3f} min.'.format((ttotend - ttotstart) / 60)
print 'Data file:       ' + data_file
print 'Output file:     ' + output_file
if script_name != '':
   print 'GNUPLOT script:  ' + script_name
print 'Numbers of Bins: ' + str(bins[0]) + ' x ' + str(bins[1])
print 'Bin X-Range:     ' + str(binrange[0]) + ' - ' + str(binrange[1])
print 'Bin Y-Range:     ' + str(binrange[2]) + ' - ' + str(binrange[3])
print 'Points Omitted:  ' + str(discarded)
print 'Done!'
