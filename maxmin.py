#!/usr/bin/env python

import sys

def printusage():
   print 'Usage: maxmin.py -f file1 file2 file3 ... -n column'
   sys.exit()

if len(sys.argv) < 2 or sys.argv[1].startswith('-h') or sys.argv[1].startswith('--h'):
   printusage()

file = []
try:
   for x in range(1,len(sys.argv)):
      if sys.argv[x] == '-f':
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith('-'):
            file.append(sys.argv[x])
            x += 1
         x -= 1
      elif sys.argv[x] == '-n':
         col = int(sys.argv[x+1])
      elif sys.argv[x].startswith('-'):
         print 'Unknown flag ' + sys.argv[x]
         printusage()

except IndexError:
   print >> sys.stderr, 'Command line error!'
   printusage()
except ValueError:
   print >> sys.stderr, 'Column number must be an integer!'
   printusage()

unset = True
for i in range(len(file)):
   try:
      data = open(file[i],'r')
   except IOError:
      print >> sys.stderr, 'Input file not found!'
      printusage()

   for line in data:
      words = line.split()
      if len(words) < col:
         continue
      
      try:
         val = float(words[col-1])
      except ValueError:
         continue
   
      if unset:
         max = val
         min = val
         unset = False
      if val > max:
         max = val
      if val < min:
         min = val
   
print 'Max is: ' + str(max)
print 'Min is: ' + str(min)
