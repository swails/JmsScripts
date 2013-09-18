#!/usr/bin/env python
from __future__ import division

###########################################################
#                                                         #
#  This script will generate a TI result for protonation  #
#  transitions from one state to another.                 #
#                                                         #
###########################################################

import sys, os
sys.path.append(os.path.join(os.getenv('AMBERHOME'), 'bin'))

# modules
from ParmedTools.ParmedActions import changeprotstate, netcharge
from chemistry.amber.readparm import AmberParm
from utilities import which
from optparse import OptionParser, OptionGroup

from time import time

try:
   import numpy as np
   from scipy import integrate
except ImportError:
   integrate = None

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

tstart = time()

epilog = '''This program will run TI calculations to determine a good estimate 
of the reference energy for a would-be titratable residue (between two specific
states). You must set the environment variable DO_PARALLEL such that sander.MPI
can be run with at least 2 threads (e.g., 'mpirun -np 2')'''

parser = OptionParser(epilog=epilog)
group = OptionGroup(parser, 'Residue Info', 'Options that pertain to the ' +
                    'nature of the residue you want to make titratable')
group.add_option('-r', '--resname', dest='resname', metavar='RESIDUE_NAME',
                 help='Name of the residue you wish to titrate. No default',
                 default=None)
group.add_option('-0', '--state0', dest='state0', metavar='INT', default=0,
                 help='The original protonation state you want as lambda=0. ' +
                 '(Default %default)', type='int')
group.add_option('-1', '--state1', dest='state1', metavar='INT', default=1,
                 help='The original protonation state you want as lambda=1. ' +
                 '(Default %default)', type='int')
group.add_option('-n', '--nucleic-acid', dest='amino_acid', default=True,
                 action='store_false', help='This residue is a nucleic acid')
group.add_option('-a', '--amino-acid', dest='amino_acid', action='store_true',
                 help='This residue is an amino acid. (Default)')
group.add_option('--left-res', dest='left_res', metavar='RES', default=None,
                 help='Which residue comes to the left in the model compound ' +
                 'Defaults to the basic caps for AA or NA')
group.add_option('--right-res', dest='right_res', metavar='RES', default=None,
                 help='Which residue comes to the right in the model ' +
                 'compound. Defaults to the basic caps for AA or NA')
group.add_option('-i', '--isolated-residue', dest='isolated', default=False,
                 action='store_true', help='This is an isolated ligand-like ' +
                 'residue, so do not cap it. NOT default.')
group.add_option('-l', '--off-lib', dest='off', default=None, metavar='FILE',
                 help='For custom residues, the OFF file with the residue ' +
                 'definition')
group.add_option('-f', '--frcmod', dest='frcmod', default=None, metavar='FILE',
                 help='For custom residues, the frcmod file with the extra ' +
                 'parameters necessary for this residue')
parser.add_option_group(group)
group = OptionGroup(parser, 'Calculation Info', 'Options that pertain to the ' +
                    'conditions that you wish to titrate your residue under.')
group.add_option('-g', '--igb', dest='igb', metavar='INT', default=0,
                 help='Amber igb (GB) model to titrate for/with. No default',
                 type='int')
group.add_option('-t', '--nstlim', dest='nstlim', default=5000000, type='int',
                 metavar='INT', help='Number of time steps for simulations.')
group.add_option('-d', '--intdiel', dest='intdiel', default=1.0, type='float',
                 metavar='FLOAT', help='Internal dielectric constant to use. '
                 'Default is %default')
parser.add_option_group(group)

opt, arg = parser.parse_args()

# Initial setup
amino_acid = opt.amino_acid
isolated = opt.isolated
state1, state2 = opt.state0, opt.state1
igb = opt.igb
intdiel = opt.intdiel
resname = opt.resname.strip('"').strip("'")

if not opt.igb or not opt.resname:
   print 'Missing arguments.'
   parser.print_help()
   sys.exit(1)
if arg:
   print 'Extra arguments [%s]' % ', '.join(arg)
   sys.exit(1)

# Get the MPI command from DO_PARALLEL, as is custom with AMBER
mpi_cmd = os.getenv('DO_PARALLEL')
if mpi_cmd is None:
   sys.exit('Error: You must set DO_PARALLEL to run sander.MPI with 2 threads!')

# MDIN file for TI calculations
TI_mdin = """TI calculation
&cntrl
   nstlim = {4}, nscm=2000,
   ntx={0}, irest={1}, ntpr=1000,
   tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,
   ntb=0, igb={2}, cut=999.0,
   dt=0.001, intdiel={5},
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
   intdiel={1},
/
""".format(igb, intdiel)


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

if sander is None or tleap is None or sandermin is None:
   sys.exit("Error: You need sander.MPI and tleap to run this program!")

# Make tleap script and build the residue
extras = ''
if opt.off:
   extras += 'loadOFF %s\n' % opt.off
if opt.frcmod:
   extras += 'loadAmberParams %s\n' % opt.frcmod

if isolated:
   prev_str, next_str = '', ''
   tleap_script = """source leaprc.constph
%s
set default pbradii %%s
saveamberparm %s %s0.prmtop %s0.inpcrd
quit
""" % (extras, resname, resname.lower(), resname.lower())

   if igb == 1:
      tleap_script = tleap_script % 'mbondi'
   elif igb == 7:
      tleap_script = tleap_script % 'bondi'
   elif igb == 8:
      tleap_script = tleap_script % 'mbondi3'
   else:
      tleap_script = tleap_script % 'mbondi2'

elif amino_acid:
   prev_str, next_str = opt.left_res, opt.right_res
   if opt.left_res is None: prev_str = 'ACE'
   if opt.right_res is None: next_str = 'NME'
   tleap_script = """source leaprc.constph
%s
l = sequence {%s %s %s}
set default pbradii %%s
saveamberparm l %s0.prmtop %s0.inpcrd
quit
""" % (extras, prev_str, resname, next_str, resname.lower(), resname.lower())
   if igb == 1:
      tleap_script = tleap_script % 'mbondi'
   elif igb == 7:
      tleap_script = tleap_script % 'bondi'
   elif igb == 8:
      tleap_script = tleap_script % 'mbondi3'
   else:
      tleap_script = tleap_script % 'mbondi2'

else:
   next_str, prev_str = opt.left_res, opt.right_res
   if opt.left_res is None: next_str = 'CH3'
   if opt.right_res is None: prev_str = 'MOC'

   tleap_script = """source leaprc.constph
%s
l = sequence {%s %s %s}
set default pbradii %%s
saveamberparm l %s0.prmtop %s0.inpcrd
quit
""" % (extras, prev_str, resname, next_str, resname.lower(), resname.lower())

   if igb == 1:
      tleap_script = tleap_script % 'mbondi'
   elif igb == 7:
      tleap_script = tleap_script % 'bondi'
   elif igb == 8:
      tleap_script = tleap_script % 'mbondi3'
   else:
      tleap_script = tleap_script % 'mbondi2'

open("tleap.in","w").write(tleap_script)
os.system(tleap + " -f tleap.in > tleap.log")

prm = AmberParm("%s0.prmtop" % resname.lower())
prm.overwrite = True # allow the prmtop to be overwritten

# check validity of prmtop

if not prm.valid:
   sys.exit('Error: Invalid prmtop (%s)!' % prm)

if isolated or not prev_str:
   act = changeprotstate(prm, ':1 %s' % opt.state0)
else:
   act = changeprotstate(prm, ':2 %s' % opt.state0)
act.execute()
prm.writeParm(str(prm))
chg1 = netcharge(prm, '*').execute()

if isolated or not prev_str:
   act = changeprotstate(prm, ':1 %s' % opt.state1)
else:
   act = changeprotstate(prm, ':2 %s' % opt.state1)
act.execute()
prm.writeParm('%s1.prmtop' % resname.lower())
chg2 = netcharge(prm, '*').execute()

# Create the new prmtops and print out net charge of each state, then delete the prmtop object

idx = prm.parm_data['RESIDUE_LABEL'].index(opt.resname)
start_res = prm.parm_data['RESIDUE_POINTER'][idx] - 1

print 'State 0: Net charge is %.4f' % chg1
print 'State 1: Net charge is %.4f' % chg2

# Now it's time to minimize
open('min.mdin','w').write(Min_mdin)

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
      open('mdin','w').write(TI_mdin.format(1,0,igb,0.1*i,opt.nstlim,intdiel))
      open('groupfile','w').write(TI_groupfile1)
   else:
      open('mdin','w').write(TI_mdin.format(5,1,igb,0.1*i,opt.nstlim,intdiel))
      open('groupfile','w').write(TI_groupfile2.format(resname.lower(),i-1,i))

   os.system('{0} {1} -ng 2 -groupfile groupfile'.format(mpi_cmd, sander))

# Now it's time to construct the final data file. This is hacked together from what
# I used to do with simple bash commands

os.system('rm -f profile.dat')

for i in range(11):
   os.system('grep -A 8 "A V E R A G E" %s_0.mdout | tail -1 | awk \'{print $3}\' >> profile.dat' % i)

profile_lines = open('profile.dat','r').readlines()
ofile = open('profile.dat','w')

if integrate is not None:
   xdata = np.zeros(11)
   ydata = np.zeros(11)

for i in range(11):
   ofile.write("{0} {1}".format(i * 0.1, profile_lines[i]))
   if integrate is not None:
      xdata[i] = i * 0.1
      ydata[i] = profile_lines[i]

# Now integrate the data if possible
if integrate is not None:
   integral = integrate.simps(ydata, x=xdata)
   open('FINAL_INTEGRAL', 'w').write("Final integral (via Simpsons' rule) is %f"
                                     % integral + '\n')

ofile.close()

tend = time()

print "Done! Integrate profile.dat to find the reference energy.\n"
print "It took {0} minutes".format((tend - tstart) / 60)
