#!/usr/bin/env python

import sys

def printusage():
   print 'Usage: maxmin.py -f file -n column'
   sys.exit()

if len(sys.argv) < 2 or sys.argv[1].startswith('-h') or sys.argv[1].startswith('--h'):
   printusage()

try:
   for x in range(len(sys.argv)):
      if sys.argv[x] == '-f':
         file = sys.argv[x+1]
      elif sys.argv[x] == '-n':
         col = int(sys.argv[x+1])
      elif sys.argv[x].startswith('-h'):
         print 'Unknown flag ' + sys.argv[x]
         printusage()

except IndexError:
   print >> sys.stderr, 'Command line error!'
   printusage()
except ValueError:
   print >> sys.stderr, 'Column number must be an integer!'
   printusage()

try:
   data = open(file,'r')
except IOError:
   print >> sys.stderr, 'Input file not found!'
   printusage()

unset = True
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
