#!/usr/bin/env python

################################################################################
#                                                                              #
#  This script will integrate a ramachandran plot within 2 ranges.             #
#                                                                              #
################################################################################

import sys, os, math

def printusage():
   import sys
   print 'integrateRama.py -i Ramachandran_datafile -phi NxN -psi NxN'
   sys.exit()

if len(sys.argv) != 7:
   printusage()

phi = [0,0]
psi = [0,0]

for x in range(len(sys.argv)):
   if sys.argv[x] == '-i':
      datafile = sys.argv[x+1]
   elif sys.argv[x] == '-phi':
      phi[0] = float(sys.argv[x+1].split('x')[0].strip())
      phi[1] = float(sys.argv[x+1].split('x')[1].strip())
   elif sys.argv[x] == '-psi':
      psi[0] = float(sys.argv[x+1].split('x')[0].strip())
      psi[1] = float(sys.argv[x+1].split('x')[1].strip())

sum = 0

data = open(datafile,'r')

for line in data:
   if len(line.strip()) == 0:
      continue
   words = line.split()
   for x in range(len(words)):
      words[x] = float(words[x])

   if words[0] > phi[0] and words[0] < phi[1] and words[1] > psi[0] and words[1] < psi[1]:
      sum += words[2]

print sum
