#!/usr/bin/env python

###############################################################################
#                                                                             #
# This utility will (hopefully) return a state energy value for a specific    #
# residue for use in constant pH simulations.                                 #
#                                                                             #
###############################################################################

import math, os, sys

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)

def fileexists(filename):
   try:
      file = open(filename,'r')
   except IOError:
      return False
   file.close()
   return True

os.system('rm _TITR_*')

# Edit these files for your system

igb           = 5       # GB model to use
pka           = 4.0     # The pKa of the residue
residuename   = "AS4"   # The name of the residue to titrate
repetitions   = 5       # how many times to titrate and average
tolerance     = 0.008   # how tolerant we should be before we call it a success
maxcycles     = 1       # how many times to iterate through to find the statene
ns            = 0.5     # how many ns each titration is

# Not frequently edited

cpinname      = "cpin"  # the name of the cpin file
system_prefix = residuename.lower()

# Very unlikely to change

statenes = []   # list of the state energies
protcnts = []   # list of the corresponding protonation counts
ratio    = 1    # the ratio of protonated:deprotonated -- we want it to be 0.5
ene_line = 0    # which line of the cpin has the state energies
nstlim   = int(ns * 0.5 * 1e6) # number of iterations

# The text for the leap script to make the files
leap_script = """source leaprc.constph

l = sequence {ACE """ + residuename + " NME}" + """

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


print statenes
print protcnts

# We now have enough to start the iteration
step = 0
while (step < maxcycles and ratio - 0.5 > tolerance):

   print '\n\n Starting step {0}'.format(step+1)

   ratio = 0

   for y in range(repetitions): # loop through number of samples

      # run sander
      print 'Running sander...'
      os.system('sander -O -i _TITR_.mdin -o _TITR_.mdout -inf _TITR_.mdinfo ' +
         '-x _TITR_.mdcrd -cpin {0} -cpout _TITR_.cpout '.format(cpinname) + 
         '-cprestrt _TITR_.cprestrt -p {0}.prmtop -c {0}.inpcrd'.format(system_prefix) +
         ' -r _TITR_.restrt')

      # calculate the pKa to get the ratio of protonated to unprotonated
      print 'Running calcpka.pl...'
      os.system("calcpka.pl {0} _TITR_.cpout | head -n 1 |".format(cpinname) + \
                " awk '{print $9}' > _TITR_pka.dat")

      # get the new ratio from that file
      get_ratio = open('_TITR_pka.dat','r')
      ratio += float(get_ratio.readline().strip())
      get_ratio.close()

   ratio /= float(repetitions)
   print 'New ratio is {0}'.format(ratio)

   print 'Calculating new state energies...'

   print 'New state energy array is '
   print statenes

   step += 1


statenestring = 'State Energies: '
for x in range(len(statenes)):
   statenestring += str(statenes[x]) + ', '

output = open("FINAL_RESULTS.txt","w")
output.write(statenestring + 'Final Ratio: {0}'.format(ratio))
output.close()

print statenestring + 'Final Ratio: {0}'.format(ratio)
