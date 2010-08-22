#!/usr/bin/env python

###############################################################################
#                                                                             #
# This utility will (hopefully) return a state energy value for a specific    #
# residue for use in constant pH simulations.                                 #
#                                                                             #
###############################################################################

import math, os, sys, time

tstart = time.time()

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)

###############################################################################

#  User-defined functions
def fileexists(filename): # returns logical True if file exists, False if not
   try:
      file = open(filename,'r')
   except IOError:
      return False
   file.close()
   return True

def printusage():
   print 'Usage: GetEnergy.py -igb <igb> -pka <pka> -resname <residue_name> -maxcyc <number_iterations> \\'
   print '                    -ns <nano-s of simulation per rep> -cpin <cpin name> -reps <repetitions>'
   sys.exit()

###############################################################################

if len(sys.argv) > 1:
   if sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--help':
      printusage()

os.system('rm -f _TITR_*')

# Edit these files for your system

igb           = 5       # GB model to use
pka           = 4.0     # The pKa of the residue
residuename   = "AS4"   # The name of the residue to titrate
repetitions   = 5       # how many times to titrate and average
tolerance     = 0.008   # how tolerant we should be before we call it a success
maxcycles     = 10      # how many times to iterate through to find the statene
ns            = 2       # how many ns each titration is
cpinname      = "cpin"  # the name of the cpin file

for x in range(1,len(sys.argv)):
   if sys.argv[x] == '-igb':
      igb = int(sys.argv[x+1])
   elif sys.argv[x] == '-pka':
      pka = float(sys.argv[x+1])
   elif sys.argv[x] == '-resname':
      residuename = sys.argv[x+1]
   elif sys.argv[x] == '-maxcyc':
      maxcycles = int(sys.argv[x+1])
   elif sys.argv[x] == '-ns':
      ns = float(sys.argv[x+1])
   elif sys.argv[x] == '-cpin':
      cpinname = sys.argv[x+1]
   elif sys.argv[x] == '-reps':
      repetitions = int(sys.argv[x+1])
   elif sys.argv[x].startswith('-'):
      print 'Unknown command-line argument ' + sys.argv[x]
      printusage()

# Not frequently edited

system_prefix = residuename.lower()

# Very unlikely to change

statenes = []   # list of the state energies
protcnts = []   # list of the corresponding protonation counts
ratiohist= []   # the history of the ratios
enehist  = []   # the history of the energies
ratio    = 1    # the ratio of protonated:deprotonated -- we want it to be 0.5
ene_line = 0    # which line of the cpin has the state energies
nstlim   = int(ns * 0.5 * 1e6) # number of iterations
kboltz   = 0.00199 # boltzmann's constant in kcal/mol K
temp     = 300  # temperature
beta     = 1 / (kboltz * temp) # 1 / k T
zero     = 1e-8

# The text for the leap script to make the files
if amino_acid:
   leap_script = """source leaprc.constph

l = sequence {ACE """ + residuename + " NME}" + """

saveamberparm l {0}.prmtop {0}.inpcrd
quit
""".format(system_prefix)
else:
   leap_script = """source leaprc.constph

l = sequence {MOC """ + residuename + " CH3}" + """

saveamberparm l {0}.prmtop {0}.inpcrd
quit
""".format(system_prefix)

titrate_mdin = """Implicit solvent constant pH molecular dynamics
&cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1000, ntwx=1000, nstlim={0},
   dt=0.002, ntt=3, tempi=300,
   temp0=300, tautp=2.0, ig=-1,
   ntp=0, ntc=2, ntf=2, cut=30,
   ntb=0, igb={1}, saltcon=0.1,
   nrespa=1, tol=0.000001, icnstph=1,
   solvph=0, ntcnstph=5, gamma_ln=2.0,
   ntwr=500, ioutfm=1,
/
""".format(nstlim,igb)

# Run leap really quick and create the mdin file

leapin = open('_TITR_leap.in','w')
leapin.write(leap_script)
leapin.close()

if not fileexists('{0}.prmtop'.format(system_prefix)) and \
     not fileexists('{0}.inpcrd'.format(system_prefix)):
   os.system('tleap -f _TITR_leap.in')

mdin = open('_TITR_.mdin','w')
mdin.write(titrate_mdin)
mdin.close()

# create the cpin file

if not fileexists(cpinname):
   os.system('cpinutil.py -p {0}.prmtop -igb {1} > {2}'.format(system_prefix, igb,
             cpinname) )

# read the cpin file

cpin = open(cpinname,'r')
cpinlines = cpin.readlines()
cpin.close()

# Get the state energies and protonation counts
for x in range(len(cpinlines)):
   if cpinlines[x][0:9] == ' STATENE=': # we've found it
      ene_line = x # store which line it occurs on
      statenes = cpinlines[x][9:].split(',') # store the state energies
      statenes.pop()

   elif cpinlines[x][0:9] == ' PROTCNT=': 
      protcnts = cpinlines[x][9:].split(',') # store the protonation counts
      protcnts.pop()

for x in range(len(statenes)): # make sure they're all floats and ints
   statenes[x] = float(statenes[x].strip())
   protcnts[x] = int(protcnts[x].strip())


# find the number of deprotonated states and the number of protonated states:
if len(statenes) != len(protcnts):
   print 'Error: statene list is not equal to protcnt list in the CPIN file!'
   sys.exit()

min,max = 99999,-99999     # min and max number of protons in all states
protnum,deprotnum = [],[]  # the array of indices that correspond to prot, deprot states

for x in range(len(protcnts)): # find the deprotonated and protonated forms states by max/min num of protons
   if protcnts[x] < min:
      min = protcnts[x]
   if protcnts[x] > max:
      max = protcnts[x]

for x in range(len(protcnts)): # build list of prot, deprot state indices
   if protcnts[x] == min:
      deprotnum.append(x)
   if protcnts[x] == max:
      protnum.append(x)

# make sure that either protonated or deprotonated has energy of zero

zeroene = 'deprot'
for x in range(len(deprotnum)):
   if statenes[deprotnum[x]] != 0.0:
      zeroene = 'prot'
for x in range(len(protnum)):
   if zeroene == 'prot' and statenes[protnum[x]] != 0:
      print 'Error: Either deprot or prot must have zero energies!'
      sys.exit()

enehist.append(statenes)
# debug printing
print statenes
print protcnts

# set up some numbers we need for calculating new energies
nd = float(len(deprotnum))
np = float(len(protnum))

# We now have enough to start the iteration
step = 0
while (step < maxcycles and abs(ratio - 0.5) > tolerance):

   print '\n\n Starting step {0}'.format(step+1)

   ratio = 0

   for x in range(repetitions): # loop through number of samples

      # run sander
      print 'Running sander...'
      os.system('sander -O -i _TITR_.mdin -o _TITR_.mdout -inf _TITR_.mdinfo ' +
         '-x _TITR_.mdcrd -cpin {0} -cpout _TITR_.cpout '.format(cpinname) + 
         '-cprestrt _TITR_.cprestrt -p {0}.prmtop -c {0}.inpcrd'.format(system_prefix) +
         ' -r _TITR_.restrt')

      # calculate the pKa to get the ratio of protonated to unprotonated
      print 'Running calcpka...'
      os.system("calcpka {0} _TITR_.cpout | head -n 2 | tail -n 1 |".format(cpinname) + \
                " awk '{print $9}' > _TITR_pka.dat")

      # get the new ratio from that file
      get_ratio = open('_TITR_pka.dat','r')
      tmp = float(get_ratio.readline().strip())
      get_ratio.close()

      # print the ratio and add it to ratio for averaging
      print tmp 
      ratio += tmp

   ratio /= float(repetitions)
   print 'New ratio is {0}'.format(ratio)
   ratiohist.append(ratio)

   if abs(ratio - 0.5) < tolerance: # don't calc new energy if converged
      break

   print 'Calculating new state energies...'

   subfactor = beta * math.log(ratio/0.5) * nd / np

   if zeroene == 'deprot':
      for x in range(len(protnum)):
         statenes[protnum[x]] -= subfactor
   else:
      for x in range(len(deprotnum)):
         statenes[deprotnum[x]] += subfactor

   print 'New state energy array is '
   print statenes
   enehist.append(statenes)

   print 'Writing new cpin file...'

   os.system('rm -f {0}'.format(cpinname))
   cpin = open(cpinname,'w')

   stateneline = ' STATENE='
   for x in range(len(statenes)):
      stateneline += '{0:.5f},'.format(statenes[x])
   stateneline += '\n'

   for x in range(len(cpinlines)):
      if x != ene_line:
         cpin.write(cpinlines[x])
      else:
         cpin.write(stateneline)

   cpin.close()

   step += 1


statenestring = 'State Energies: '
for x in range(len(statenes)):
   statenestring += '{0:.5f}, '.format(statenes[x])

output = open("FINAL_RESULTS.txt","w")
output.write(statenestring + 'Final Ratio: {0}\n'.format(ratio))
output.close()

tend = time.time()

print '\n'
print enehist
print ratiohist
print statenestring + 'Final Ratio: {0}'.format(ratio)

print '\nTiming: {0:.2f} hr.'.format((tend-tstart)/3600)
