#!/usr/bin/env python

##########################################################################
#
# This program emulates the -nxy option of xmgrace and prepares a gnuplot
# script that can be immediately read into gnuplot (and added upon) to 
# plot the data (1st column is always x, each subsequent column is the y
# value for another data set)
#
#   Jason Swails 8/16/2009
#
##########################################################################

import sys

gnuscriptname = 'gnu.in'

for x in range(len(sys.argv)):
   if '-help' in sys.argv[x].lower():
      print 'Usage: gnuplotnxy.py <file.dat> {<gnuplot script>}'
      sys.exit()

if len(sys.argv) == 1 or len(sys.argv) > 3:
   print 'Make sure that you have data on one of your first three lines'
   print 'Usage: gnuplotnxy.py <file.dat> {<gnuplot script>}'
   sys.exit()

try:
   datafile = open(sys.argv[1],'r')
except IOError:
   print 'Outputfile ' + sys.argv[1] + ' does not exist!'
   print 'Usage: gnuplotnxy.py <file.dat> {<gnuplot script>}'
   sys.exit()

if len(sys.argv) == 3:
   gnuscriptname = sys.argv[2]

lines = [len(datafile.readline().split()), len(datafile.readline().split()), len(datafile.readline().split())]
datafile.close()

maxlen = max(lines)
gnuin = open(gnuscriptname, 'w')

gnuin.write('unset key\nplot ')
for x in range(2,maxlen+1):

   if x != maxlen:
      gnuin.write('\'' + sys.argv[1] + '\' u 1:' + str(x) + ' with lines, \\\n')
   else:
      gnuin.write('\'' + sys.argv[1] + '\' u 1:' + str(x) + ' with lines\n')

gnuin.close()
