#!/usr/bin/env python

# cut up the cpout files into various chunks

import sys, os

try:
   help = sys.argv[1]
except IndexError:
   help = '-help'

if '-help' in help:
   print 'Usage: pkachunks.py interval prefix cpin cpout1 cpout2 ... cpoutN'
   sys.exit()

interval = int(sys.argv[1])
counter = 0
prefix = sys.argv[2]
cpouts = []
filecounter = 0
cpin = sys.argv[3]

for x in range(4,len(sys.argv)):
   cpouts.append(open(sys.argv[x],'r'))

filecounter = filecounter + 1
file = open(prefix + '.' + str(filecounter), 'w')
for x in range(len(cpouts)):
   for line in cpouts[x]:
      if 'Solvent' in line:
         counter = counter + 1
         if counter % interval == 0:
            file.close()
            filecounter = filecounter + 1
            file = open(prefix + '.' + str(filecounter),'w')
      file.write(line)

file.close()

for x in range(filecounter):
   y = x + 1
   os.system('calcpka.pl ' + cpin + ' ' + prefix + '.' + str(y) + ' > ' + prefix + '.' + str(y) + '.pka')

file = open(prefix + '.total.dat' , 'w')


for x in range(filecounter):
   y = x + 1
   pkafile = open(prefix + '.' + str(y) + '.pka','r')
   pkafilelines = pkafile.readlines()
   pkafile.close()
   
   if x == 0:
      numres = len(pkafilelines) - 3

   file.write(str(y) + ' ')
   for z in range(numres):
      pkawords = pkafilelines[z].split()
      file.write(pkawords[5] + ' ')
   file.write('\n')

