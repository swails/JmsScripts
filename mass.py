#!/usr/bin/env python

####################################################
#                                                  #
# This utility prints the total mass of a system   #
# described by the prmtop.                         #
#                                                  #
#    Written by Jason Swails, 11/21/2009           #
#                                                  #
####################################################


import sys, utilities

if len(sys.argv) != 2 or '-help' in sys.argv[1]:
   print 'mass.py <prmtop>'
   sys.exit()

mass = 0

try:
   prmtop = open(sys.argv[1],'r')
except IOError:
   print >> sys.stderr, "Error: Topology file " + sys.argv[1] + " does not exist!"
   sys.exit()

lines = prmtop.readlines()
prmtop.close()
residues = []

for x in range(len(lines)):
   if 'MASS' in lines[x]:
#     skip any comment lines
      while not 'FORMAT' in lines[x]:
         x = x + 1
#     skip past the FORMAT statement
      x = x + 1
#     read in all of the residues and append them to the 'residue' list
      while not '%FLAG' in lines[x]:
         words = lines[x].split()
         for y in range(len(words)):
            mass += float(words[y])
         x = x + 1
#     after this while loop, all of the residues are loaded into the list
      break


print >> sys.stdout, 'The mass of ' + sys.argv[1] + ' is ' + str(mass) + ' g/mol.'
