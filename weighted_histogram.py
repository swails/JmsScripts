#!/usr/bin/env python

"""
This will (hopefully) analyze a weighted histogram using a density of states
as input, which it will then convert to a Temperature-dependent distribution
based on the input temperature.
"""

from optparse import OptionParser, OptionGroup
import os, sys, math
import numpy as np

KB = 0.0019872041
debug_printlevel = 0

# Set up the new excepthook to control how fatal exceptions are handled

def excepthook(exception_type, exception_value, tb):
   """ Replaces sys.excepthook so fatal exceptions kill all MPI threads and
       we can control the printing of tracebacks. Those are helpful for 
       debugging purposes, but may be unsightly to users. debug_printlevel
       set above controls this behavior
   """
   import traceback
   if debug_printlevel > 0: traceback.print_tb(tb)
   sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
   sys.stderr.write('Goodbye.\n')
   sys.exit(1)

sys.excepthook = excepthook # replace sys.excepthook with my definition above

class WeightedHistogramError(Exception):
   """ For errors in this script """

parser = OptionParser()

group = OptionGroup(parser, 'Input Options', 'Specify input variables')
group.add_option('-d', '--density-of-states', default=None, dest='density',
                 metavar='FILE', help='File containing the density of states.')
group.add_option('-p', '--property-file', default=None, dest='property',
                 metavar='FILE', help='File containing synchronized time-' +
                 'series of the property of interest with the energy.')
group.add_option('-t', '--temperature', dest='temp', metavar='FLOAT',
                 default=0.0, type='float', help='Temperature to get ' +
                 'distribution for. No default. Must be > 0.')
group.add_option('--debug', default=False, dest='debug', action='store_true',
                 help='Print debugging information and full tracebacks')
parser.add_option_group(group)

group = OptionGroup(parser, 'Histogram Options', 'Specify histogram settings')
group.add_option('-n', '--no-normalize', default=False, action='store_true',
                 help='Do not normalize the data', dest='normalize')
group.add_option('--min', dest='min', default=None, metavar='FLOAT',
                 help='Minimum value for histogram -- defaults to minimum of ' +
                 'data set')
group.add_option('--max', dest='max', default=None, metavar='FLOAT',
                 help='Maximum value for histogram -- defaults to maximum of ' +
                 'data set')
group.add_option('-s', '--spacing', dest='spacing', metavar='FLOAT',
                 default=None, help='Spacing between histograms.')
group.add_option('-b', '--number-of-bins', dest='bins', metavar='FLOAT',
                 default=None, help='Number of bins. Cannot be used with ' +
                 '--spacing')
parser.add_option_group(group)

group = OptionGroup(parser, 'Output Options', 'Specify output variables')
group.add_option('-o', '--output', dest='output', metavar='FILE', default=None,
                 help='Output file. Defaults to stdout.')
parser.add_option_group(group)

opt,arg = parser.parse_args()

# Load the debug print level
debug_printlevel = int(opt.debug)

# Make sure we have a density file
if not opt.density:
   raise WeightedHistogramError('Missing density of states file!')
elif not os.path.exists(opt.density):
   raise WeightedHistogramError('Could not find density of states file %s!' %
                                opt.density)

# Make sure we have our properties
if not opt.property:
   raise WeightedHistogramError('Missing property file!')
elif not os.path.exists(opt.density):
   raise WeightedHistogramError('Could not find property file %s!' %
                                opt.property)

# Make sure we have a valid temperature
if opt.temp <= 0.0:
   raise WeightedHistogramError('No (or bad) temperature specified! ' +
         'Temperature must be > 0 K')
KbT = KB * opt.temp

# Make sure our histogramming defaults are OK
if opt.spacing and opt.bins:
   raise WeightedHistogramError('Bin spacing and bin numbers are mutually ' +
         'exclusive!')

# Open our output file
if opt.output:
   output = open(opt.output, 'w')
else:
   output = sys.stdout

# Now read in our density of states { ln (P[Q] from Dan's program) }
density = open(opt.density, 'r')

# Count how many lines we have then rewind the file
num_density_bins = 0
while density.readline(): num_density_bins += 1
density.seek(0)

# Initialize our numpy arrays
density_of_states = np.zeros((num_density_bins,2))

# Fill our density of states. Then determine its range and interval, assuming
# the list is ordered and has regular intervals
for i, line in enumerate(density):
   words = line.split()
   density_of_states[i][0] = float(words[0])
   density_of_states[i][1] = float(words[1])
density_of_states.minE = density_of_states[0][0]
density_of_states.maxE = density_of_states[num_density_bins-1][0]
density_of_states.interval = ((density_of_states.maxE - density_of_states.maxE)/
                               num_density_bins)

# Now read in all of the properties and load them into numpy arrays
property = open(opt.property, 'r')
num_properties = 0
while property.readline(): num_properties += 1
property.seek(0)

property_data = np.zeros((num_properties,2))

for i, line in enumerate(property):
   words = line.split()
   property_data[i][0] = float(words[0])
   property_data[i][1] = float(words[1])

# Now build the log of the denominator -- this is given by the formula
#         N                                  N  
#  log( sum (a_i) )  = log(a_0) + log (1 + sum ( exp[log(a_i)-log(a_0)] ) )
#         i=0                                i=1

log_a0 = density_of_states[0][1]
sum2 = 0
for i in range(1,len(density_of_states)):
   sum2 += math.exp( density_of_states[i][1]-density_of_states[0][1] )
denominator = log_a0 + math.log(1+sum2)

if opt.debug: print >> sys.stderr, 'Log of denominator is: %f' % denominator

# Now we have to set up the property's histograms. We will define the ranges
# using defaults if the user didn't provide us information

# Now we have a bunch of properties that correspond with the energies for our
# expanded ensemble. What we want to do is reweight each point properly within
# the ensemble defined by our new temperature. We do this using the formula:

# O(E) n(E) exp(-E/KbT)
# ---------------------
#     denominator

# We have to find which energy "bin" the Property (O(E)) belongs to, and add
# the weighted factor to that bin. So, what we have to do now is to set up the
# property histogram

# Determine the max/min
if opt.min:
   mymin = opt.min
else:
   mymin = min([prop[1] for prop in property_data])
if opt.max:
   mymax = opt.max
else:
   mymax = max([prop[1] for prop in property_data])

# Make the max and min nice if we set the defaults
power = int(math.log10(abs(mymax)))
if power > 1: power = 1
if not opt.min:
   mymin = math.floor(mymin*10**power) / 10**power
   if opt.debug: print >> sys.stderr, 'Default minimum is %f' % mymin
if not opt.max:
   mymax = math.ceil(mymin*10**power) / 10**power
   if opt.debug: print >> sys.stderr, 'Default maximum is %f' % mymax

binrange = mymax - mymin

# Determine the number of bins
if opt.bins:
   nbins = opt.bins
   spacing = binrange / nbins
else:
   if opt.spacing:
      # Make sure we don't have a stupid spacing
      if opt.spacing >= mymax - mymin:
         raise WeightedHistogramError('Spacing is larger than your data range!')
      spacing = opt.spacing
   else:
      # Use "Scott's Choice" for spacing
      spacing = property_data.std(0)[0] * 3.5 / len(property_data)**(1/3.0)
   nbins = (mymax - mymin) // spacing

# Load the property histogram with energies as the first element in each row
prop_hist = np.zeros(nbins)

normal_fac = 0 # normalization factor
omitted_nums = 0

# Now go through every property and apply the proper weight based on the energy
# (the value in the first dimension of the property_data).

for i, item in enumerate(property_data):
   # Determine which energy bin we are (add lower than minimum to 0 bin and
   # higher than maximum to last bin)
   if item[0] < density_of_states.minE:
      energy_bin = density_of_states[0]
   elif item[0] > density_of_states.maxE:
      energy_bin = density_of_states[len(density_of_states)-1]
   else:
      energy_bin = density_of_states[int((item[0]-density_of_states.minE)/
                                          density_of_states.interval)]
   if opt.debug:
      print >> sys.stderr, 'Energy diff is %f' % (energy_bin[0]-item[0])
   wt = math.exp(energy_bin[1] + KbT*energy_bin[0] - denominator)
   normal_fac += wt

   # Now add this to the property histogram
   if item[1] > mymax or item[1] < mymin:
      omitted_nums += 1
      continue
   prop_bin = (item[1]-mymin)//spacing
   prop_hist[prop_bin] += wt

print >> sys.stderr, '# Done histogramming! Normalization factor is %f' % (
         normal_fac)
if opt.normalize:
   print >> sys.stderr, '# Normalizing...'
   for i, val in enumerate(prop_bin): prop_bin[i] = val / normal_fac
else:
   print >> sys.stderr, '# Not normalizing'

for i, val in range(prop_bin):
   output.write('%f %f\n' % (mymin+i*spacing, val))
