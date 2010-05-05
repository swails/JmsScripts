#!/usr/bin/env python

###################################################################
#                                                                 #
# This program will integrate a series of points using the trap-  #
# ezoid rule (average of left and right rectangle method). It is  #
# less accurate than the Simpsons method or midpoint method, but  #
# it is good for TI when both ends (lambda = 0 and lambda = 1)    #
# are in the data set.                                            #
#                                                                 #
# Last update: 11/21/2009   (by jms)                              # 
#  by Jason Swails                                                #
#                                                                 #
################################################################### 

# import needed modules

import math, sys

########################### USER-DEFINED FUNCTIONS ################

def PrintUsage():
   print 'Usage: Integrate_Trap.py <data_file>\n'
   print '<data_file> must be 2 columns space-delimited <x> <y>'

###################################################################


# get command-line arguments

if len(sys.argv) != 2 or '-help' in sys.argv[1]:
   PrintUsage
   sys.exit()

# open the file, load the data into memory, and then close the file
try:
   datafile = open(sys.argv[1],'r')
except IOError:
   print >> sys.stderr, "Error: Data file '" + sys.argv[1] + "' does not exist!"
   PrintUsage
   sys.exit()

lines = datafile.readlines()
datafile.close()
# end reading data

# strip the lines so that the first two and last two lines are both data.
end = len(lines) - 1
secondlast = len(lines) - 2

while len(lines[0].split()) != 2:
   lines.pop(0)
while len(lines[1].split()) != 2:
   lines.pop(1)
while len(lines[end].split()) != 2:
   lines.pop(end)
   end = len(lines) - 1
   secondlast = len(lines) - 2
while len(lines[secondlast].split()) != 2:
   lines.pop(secondlast)
   secondlast = len(lines) - 2

# initialize the sums
line1 = lines[len(lines)-1].split()
line2 = lines[len(lines)-2].split()

right = float(line2[1]) * (float(line2[0]) - float(line1[0]))

line1 = lines[0].split()
line2 = lines[1].split()

left = float(line1[1]) * (float(line2[0]) - float(line1[0]))


# now cycle through all of the data
for x in range(1, len(lines) - 1):

   # Get rid of any blank or comment lines
   while len(lines[x].split()) != 2:
      lines.pop(x)

   # load in the lines
   line1 = lines[x-1].split()
   line2 = lines[x].split()
   line3 = lines[x+1].split()

   left = left + float(line2[1]) * (float(line3[0]) - float(line2[0]))
   right = right + float(line2[1]) * (float(line2[0]) - float(line1[0]))

print 'Left-rectangle approximation  : ' + str(left)
print 'Right-rectangle approximation : ' + str(right)
print 'Trapezoidal appoximation      : ' + str((left + right) / 2)
