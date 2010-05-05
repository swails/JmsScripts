#!/usr/bin/env python

# Finds the average of a specific column of a given file

import sys, math

if len(sys.argv) < 3:
   print 'This utility computes the average and standard deviation of a given'
   print 'column of data in given data files.\n'
   print 'Usage: AvgCol.py file1 file2 ... filen column'
   sys.exit()

   
try:
   column = int(sys.argv[len(sys.argv)-1])
except ValueError:
   print 'Usage: AvgCol.py file1 file2 ... filen column'
   print str(sys.argv[len(sys.argv)-1]) + ' is not an integer!'
   sys.exit()

sum = 0
number_entries = 0
stdev = 0
sumsquared = 0

for x in range(len(sys.argv)-2):
   y = x + 1
   try:
      file = open(sys.argv[y],'r')
   except IOError:
      print 'Usage: AvgCol.py file1 file2 ... filen column'
      print sys.argv[y] + ' could not be opened!'
      continue
   
   for line in file:
      words = line.split()
      if len(words) < column:
         continue
      try:
         sum = sum + (float(words[column - 1]))
         sumsquared = sumsquared + (float(words[column - 1]) * float(words[column - 1]))
         number_entries = number_entries + 1
      except ValueError:
         continue
   
   file.close()
   
ave = sum / number_entries
ave2 = sumsquared / number_entries
std = math.sqrt(ave2 - ave * ave)
   
print 'The average of column ' + str(column) + ' is: ' + str(ave)
print 'The standard deviation of column ' + str(column) + ' is: ' + str(std)
