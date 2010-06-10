#!/usr/bin/env python

###############################################################################
#                                                                             #
#  This script will locate salt bridges in a trajectory.                      #
#                                                                             #
###############################################################################

import sys, utilities, os

acceptorlist = {'ASP' : 'CG', 'AS4' : 'CG', 'GLU' : 'CD', 'GL4' : 'CD'}
donorlist = {'AS4' : 'CG', 'ASH' : 'CG', 'GL4' : 'CD', 'GLH' : 'CD', 'LYS' : 'NZ', 'ARG' : 'CZ'}
tolerance = 5.000 # how close in angstroms atoms have to be to be called a salt bridge
fraction = 0.25  # what fraction of the time the salt bridge must be there

mdcrds = []
prmtop = 'prmtop'
output = 'salt_bridge_output'

acceptor_inds = []
donor_inds = []
pairs = []

salt_bridges = []

def clean():
   os.system('rm -f _FSB_*')

def printusage(): # usage statement
   print 'FindSaltbridge.py -p prmtop -y mdcrd1 mdcrd2 ... -o output_list -percent percent_population'
   clean()
   sys.exit()

if len(sys.argv) > 1 and (sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--help'): # ask for help
   printusage()

try:
   for x in range(len(sys.argv)):
      if sys.argv[x] == '-p':
         prmtop = sys.argv[x+1]
      elif sys.argv[x] == '-o':
         output = sys.argv[x+1]
      elif sys.argv[x] == '-percent':
         fraction = float(sys.argv[x+1])
      elif sys.argv[x] == '-y':
         x += 1
         while x < len(sys.argv) and not sys.argv[x].startswith('-'):
            mdcrds.append(sys.argv[x])
            x += 1
except IndexError:
   print 'Error: Improper command!'
   printusage()
except ValueError:
   print 'Error: "percent" must be a floating point decimal!'
   printusage()

if len(mdcrds) == 0: # if no mdcrd supplied, give default
   mdcrds.append('mdcrd')

ptraj = utilities.which('ptraj') # look for ptraj

if ptraj == 'none': 
   print 'Error: ptraj needed for FindSaltbridge.py!'
   sys.exit() # quit if not found

if utilities.fileexists(prmtop) == -1: # check for prmtop existence
   printusage()

for x in range(len(mdcrds)): # check for mdcrd existences
   if utilities.fileexists(mdcrds[x]) == -1:
      printusage()

residues = utilities.getallresinfo(prmtop,'RESIDUE_LABEL') # get all residue names

for x in range(len(residues)): # build acceptor list and donor list
   try:
      tmp = acceptorlist[residues[x]]
      acceptor_inds.append(x+1)
   except KeyError:
      pass

   try:
      tmp = donorlist[residues[x]]
      donor_inds.append(x+1)
   except KeyError:
      pass

for x in range(len(acceptor_inds)): # build list of pairs
   for y in range(len(donor_inds)):
      if acceptor_inds[x] != donor_inds[y]:
         pairs.append([acceptor_inds[x],donor_inds[y],0])

file = open('_FSB_dist.ptraj','w')

for x in range(len(mdcrds)):
   file.write('trajin {0}\n'.format(mdcrds[x]))

file.write('\n')

for x in range(len(pairs)):
   acc_index = pairs[x][0]
   acc_atom = acceptorlist[residues[acc_index-1]]
   don_index = pairs[x][1]
   don_atom = donorlist[residues[don_index-1]]

   file.write('distance {0}_{1} :{0}@{2} :{1}@{3} out _FSB_{0}_{1}.dat\n'.format(acc_index, don_index, acc_atom, don_atom))

file.close()

os.system('ptraj {0} _FSB_dist.ptraj'.format(prmtop))

print 'Done with ptraj...'

outputfile = open(output,'w',0)

for x in range(len(pairs)):
   file = open('_FSB_{0}_{1}.dat'.format(pairs[x][0],pairs[x][1]),'r')
   within = 0.0
   total = 0.0
   for line in file:
      words = line.split()
      if float(words[1]) < tolerance:
         within += 1.0
      total += 1.0
   file.close()
   if within / total >= fraction:
      outputfile.write('{0} {1} - {2} {3} : Fraction {4:.3f}\n'.format(residues[pairs[x][0]-1].ljust(3), str(pairs[x][0]).rjust(3),
                 residues[pairs[x][1]-1].ljust(3), str(pairs[x][1]).rjust(3), within / total))

outputfile.close()
