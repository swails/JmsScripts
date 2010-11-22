#!/usr/bin/env python

#++++++++++++++++++++++++++++++++++++++++++
# import modules
#++++++++++++++++++++++++++++++++++++++++++

import sys

try:
   from cpin_data import *
except:
   print >> sys.stderr, 'Error importing cpin_data.py.'
   sys.exit(1)

#++++++++++++++++++++++++++++++++++++++++++
# user-defined functions
#++++++++++++++++++++++++++++++++++++++++++

def round(number):
   """ Rounds a number to the nearest integer """
   number = float(number)
   if number >= 0:
      return int(number+0.5)
   else:
      return int(number-0.49999)

#++++++++++++++++++++++++++++++++++++++++++
# define necessary variables
#++++++++++++++++++++++++++++++++++++++++++

residues = TITRATABLE.strip().split()
problems = []

#++++++++++++++++++++++++++++++++++++++++++
# loop through residues to check
#++++++++++++++++++++++++++++++++++++++++++

print >> sys.stdout, 'Checking residues...'

for i in range(len(residues)): # i loops through all of the titr. res.

   data = getData(residues[i], 2)
   charge_count = [[],[]] # 2-D array that stores proton_count, charge
   for j in range(len(data)): # j loops through all of the states
      
      charge_count[0].append(data[j][1]) # proton count is the 2nd element after relative energy
      totcharge = sum(data[j][2:]) # charges are from the 3rd element to the end.
      charge_count[1].append(totcharge)

      if abs(totcharge - round(totcharge)) >= 0.0001:
         problems.append('The net charge of %s in state %d is %.5f' % (residues[i], j, totcharge))

   difference = charge_count[0][0] - charge_count[1][0]
   for j in range(1,len(charge_count[0])):
      if abs(charge_count[0][j] - charge_count[1][j] - difference) > 0.0001:
         problems.append('Residue %s, state %s has a proton count of %s corresponding to a charge of %.5f.\n' % (residues[i], j, charge_count[0][j], charge_count[1][j]) +
                  '\tState 0 has proton count %s with charge %s' % (charge_count[0][0], charge_count[1][0]))


#++++++++++++++++++++++++++++++++++++++++++
# print summary
#++++++++++++++++++++++++++++++++++++++++++

if len(problems) == 0:
   print >> sys.stdout, "cpin_data.py has no glaring problems."
   sys.exit(0)
else:
   print >> sys.stderr, "Issues exist with cpin_data.py...\n"
   for i in range(len(problems)):
      print >> sys.stderr, problems[i]
   sys.exit(1)
