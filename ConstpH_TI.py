#!/usr/bin/env python

###########################################################
#                                                         #
#  This script will generate a TI result for protonation  #
#  transitions from one state to another.                 #
#                                                         #
###########################################################

# modules
from cpin_data import getData
from chemistry.amber.readparm import amberParm
from utilities import which

from math import fsum
from time import time
import sys, os

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

tstart = time()

# usage instructions
def printusage():
   print >> sys.stderr, 'Usage: ConstpH_TI.py -igb <igb_value> -resname <resname> \\'
   print >> sys.stderr, '                     -states <state1> <state2> {-na || -aa}'
   sys.exit()

# Get the MPI command from DO_PARALLEL, as is custom with AMBER
try:
   mpi_cmd = os.environ["DO_PARALLEL"]
except KeyError:
   print >> sys.stderr, 'Error: You must set DO_PARALLEL to run sander.MPI with 2 threads!'
   sys.exit()

# obvious errors or call for help
if len(sys.argv) < 3:
   printusage()

if sys.argv[1] == '-h' or sys.argv[1] == '--help':
   printusage()

amino_acid = True
state1 = 0
state2 = 1

# Input parser
try:
   for i in range(len(sys.argv)):
      if sys.argv[i] == "-igb":
         igb = int(sys.argv[i+1])
      elif sys.argv[i] == "-resname":
         resname = sys.argv[i+1].strip('"').strip("'")
      elif sys.argv[i] == '-na':
         amino_acid = False
      elif sys.argv[i] == '-aa':
         amino_acid = True
      elif sys.argv[i] == '-states':
         state1 = int(sys.argv[i+1])
         state2 = int(sys.argv[i+2])
      elif sys.argv[i].startswith('-'):
         print >> sys.stderr, 'Unrecognized option, %s!' % sys.argv[i]
         printusage()
except ValueError:
   print >> sys.stderr, 'Error: IGB, STATE1, and STATE2 must be integers!'
   sys.exit()
except IndexError:
   print >> sys.stderr, 'Error: Command line error!'
   sys.exit()

# MDIN file for TI calculations

TI_mdin = """TI calculation
&cntrl
   nstlim =5000000, nscm=2000,
   ntx={0}, irest={1}, ntpr=1000,
   tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,
   ntb=0, igb={2}, cut=999.0,
   dt=0.001, 
   ntc=2, ntf=2, saltcon=0.1,
   ntwr = 10000, ntwx=1000, 
   icfe=1, clambda={3},
/
"""

Min_mdin = """Minimization to relax initial bad contacts, implicit solvent constant pH
&cntrl
   imin=1,
   ncyc=1000,
   maxcyc=5000,
   ntpr=50,
   ntb=0,
   cut=1000,
   igb={0},
/
""".format(igb)


# Groupfiles -- 1 is for the first run, 2 is for each one after
TI_groupfile1 = """-O -i mdin -p {0}0.prmtop -c {0}0.inpcrd -o 0_0.mdout -r 0_0.restrt -x 0_0.mdcrd
-O -i mdin -p {0}1.prmtop -c {0}1.inpcrd -o 0_1.mdout -r 0_1.restrt -x 0_1.mdcrd
""".format(resname.lower())

TI_groupfile2 = """-O -i mdin -p {0}0.prmtop -c {1}_0.restrt -o {2}_0.mdout -r {2}_0.restrt -x {2}_0.mdcrd
-O -i mdin -p {0}1.prmtop -c {1}_1.restrt -o {2}_1.mdout -r {2}_1.restrt -x {2}_1.mdcrd
"""

# Find sander and tleap or quit
sander = which("sander.MPI")
sandermin = which("sander")
tleap = which("tleap")
if igb == 8:
   changerad = which("ChangeParmRadii.py")
if tleap == "none":
   tleap = which("sleap")

if sander == "none" or tleap == "none" or (igb == 8 and changerad == 'none'):
   print >> sys.stderr, "Error: You need sander.MPI and tleap or sleap in your PATH to run this program!"
   sys.exit()

# Make tleap script and build the residue
if amino_acid:
   tleap_script = """source leaprc.constph
l = sequence {ACE %s NME}
saveamberparm l %s0.prmtop %s0.inpcrd
quit
""" % (resname, resname.lower(), resname.lower())
else:
   tleap_script = """source leaprc.constph
l = sequence {MOC %s CH3}
saveamberparm l %s0.prmtop %s0.inpcrd
quit
""" % (resname, resname.lower(), resname.lower())

file = open("tleap.in","w")
file.write(tleap_script)
file.close()
os.system(tleap + " -f tleap.in > tleap.log")
if igb == 8:
   os.system(changerad + ' -p {0}0.prmtop -r mbondi3'.format(resname.lower()))

# Get the data for the residue and load the amberParm object
charges = getData(resname, igb)

if charges == -1:
   print >> sys.stderr, 'Error: Could not find %s in cpin_data.py! Add the residue and re-run.'
   sys.exit()

prm = amberParm("%s0.prmtop" % resname.lower())
prm.overwrite = True # allow the prmtop to be overwritten

# check validity of prmtop

if not prm.valid:
   print >> sys.stderr, "Error: Invalid prmtop (%s)!" % prm
   sys.exit()

try:
   test = prm.parm_data["POINTERS"][0]
except KeyError:
   print >> sys.stderr, "Error: Prmtop is not valid! Check tleap.log for errors"
   sys.exit()

# Create the new prmtops and print out net charge of each state, then delete the prmtop object

start_res = prm.parm_data["RESIDUE_POINTER"][1] - 1 # starts at 0

if (prm.parm_data["RESIDUE_POINTER"][2] - prm.parm_data["RESIDUE_POINTER"][1] != len(charges[state1]) - 2):
   print >> sys.stderr, 'Warning: Residue not expected size (comparing cpin_data and prmtop).'
   sys.exit()

print >> sys.stdout, "Charge of state {0}: {1:8.4f}".format(state1, fsum(charges[state1][2:]))
print >> sys.stdout, "Charge of state {0}: {1:8.4f}".format(state2, fsum(charges[state2][2:]))

for i in range(len(charges[state1])-2):
   prm.parm_data["CHARGE"][start_res+i] = charges[state1][2+i]

prm.writeParm("%s0.prmtop" % resname.lower())

for i in range(len(charges[state2])-2):
   prm.parm_data["CHARGE"][start_res+i] = charges[state2][2+i]

prm.writeParm("%s1.prmtop" % resname.lower())
# Now it's time to minimize
file = open('min.mdin','w')
file.write(Min_mdin)
file.close()

print >> sys.stdout, "Minimizing structure..."
os.system('%s -i min.mdin -o _rm.mdout -inf _rm.mdinfo -r min.restrt -p %s0.prmtop -c %s0.inpcrd' % 
               (sandermin, resname.lower(), resname.lower()))
print >> sys.stdout, "Done minimizing. Starting TI..."
os.system('rm _rm.*')
os.system('mv min.restrt {0}0.inpcrd'.format(resname.lower()))
os.system('cp {0}0.inpcrd {0}1.inpcrd'.format(resname.lower()))

del prm

# Now it's time to start the TI

for i in range(11):
   if i == 0:
      file = open('mdin','w')
      file.write(TI_mdin.format(1, 0, igb, 0.1 * i))
      file.close()
      file = open('groupfile','w')
      file.write(TI_groupfile1)
      file.close()
   else:
      file = open('mdin','w')
      file.write(TI_mdin.format(5, 1, igb, 0.1 * i))
      file.close()
      file = open('groupfile','w')
      file.write(TI_groupfile2.format(resname.lower(),i-1,i))
      file.close()

   os.system('{0} {1} -ng 2 -groupfile groupfile'.format(mpi_cmd, sander))

# Now it's time to construct the final data file. This is hacked together from what
# I used to do with simple bash commands

os.system('rm -f profile.dat')

for i in range(11):
   os.system('grep -A 8 "A V E R A G E" %s_0.mdout | tail -1 | awk \'{print $3}\' >> profile.dat' % i)

file = open('profile.dat','r')
profile_lines = file.readlines()
file.close()
file = open('profile.dat','w')

for i in range(11):
   file.write("{0} {1}".format(i * 0.1, profile_lines[i]))

file.close()

tend = time()

print >> sys.stdout, "Done with your state energy! Integrate profile.dat to find the reference energy."
print >> sys.stdout, ""
print >> sys.stdout, "It took {0} minutes".format((tend - tstart) / 60)
