#!/usr/bin/env python

# This script will run constant pH simulations in explicit solvent using typical
# Metropolis criteria to evaluate protonation changes

import math, sys, os, random, copy

try:
   from mdin import mdin
   from readparm import amberParm
   from pbsjob import PBSjob
   from utilities import which
except ImportError:
   print >> sys.stderr, 'Error: Could not import mdin.py, readparm.py, and/or pbsjob.py! Place these'
   print >> sys.stderr, '       files in a directory referenced by PYTHONPATH.'
   sys.exit()

KB     = 0.00199
TEMP   = 300
FACTOR = math.log(10) * KB * TEMP

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# usage statement
def usage():
   print >> sys.stderr, 'cph_pmemd.py -t <total_time_steps> -dt <monte_carlo_timestep> -pH <pH value> -cpin <cpin> \\'
   print >> sys.stderr, '             -p <prmtop> -c <inpcrd> -o <mdout> -cpout <cpout> -r <restrt> -x <mdcrd>'
   sys.exit()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# evaluate metropolis monte carlo critera -- returns true to accept or false to reject
def monteCarlo(dvdl, refene, state1, state2):
   from random import random

   return (math.exp((dvdl - refene + FACTOR * (state2-state1))/(KB * TEMP)) < random() or \
         dvdl - refene + FACTOR * (state2-state1) < 0)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Parse the cpin file and load array data
def parseCpin(cpin_file, chrgdat, resstates, num_states, first_state, first_charge, \
              first_atom, num_atoms, statene, protcnt):
   cpinlines = cpin_file.readlines()
   cpin_fields = []
   inBlock = False
   for i in range(len(cpinlines)):
      if not inBlock and cpinlines[i].strip() != '&CNSTPH':
         continue
      elif not inBlock and cpinlines[i].strip() == '&CNSTPH':
         inBlock = True
      elif inBlock and (cpinlines[i].strip() == '/' or cpinlines[i].strip() == '&end'):
         inBlock = False
      elif inBlock:
         cpin_fields.extend(cpinlines[i].strip().split(','))
      else:
         return -1
   block = ''
   for i in range(len(cpin_fields)):
      if '=' in cpin_fields[i]:
         block = cpin_fields[i].split('=')[0].strip().lower()
         if block == 'chrgdat':
            chrgdat.append(float(cpin_fields[i].split('=')[1]))
         elif block == 'protcnt':
            protcnt.append(float(cpin_fields[i].split('=')[1]))
         elif block == 'resstate':
            resstates.append(int(cpin_fields[i].split('=')[1]))
         elif block == 'statene':
            statene.append(float(cpin_fields[i].split('=')[1]))
         elif block == 'trescnt':
            trescnt = int(cpin_fields[i].split('=')[1])
         elif block.endswith('first_atom'):
            first_atom.append(int(cpin_fields[i].split('=')[1]))
         elif block.endswith('first_charge'):
            first_charge.append(int(cpin_fields[i].split('=')[1]))
         elif block.endswith('first_state'):
            first_state.append(int(cpin_fields[i].split('=')[1]))
         elif block.endswith('num_atoms'):
            num_atoms.append(int(cpin_fields[i].split('=')[1]))
         elif block.endswith('num_states'):
            num_states.append(int(cpin_fields[i].split('=')[1]))
      elif len(cpin_fields[i].strip()) == 0:
         continue
      elif block == 'chrgdat':
         chrgdat.append(float(cpin_fields[i].strip()))
      elif block == 'protcnt':
         protcnt.append(int(cpin_fields[i].strip()))
      elif block == 'statene':
         statene.append(float(cpin_fields[i].strip()))
   return trescnt
         
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def updateParm(prmtop, first_charge, chrgdat, state, num_atoms, first_atom):
   for i in range(num_atoms):
#     print 'Charge of atom %s, %s is going from %8.4f to %8.4f' % (first_atom + i, prmtop.parm_data["AMBER_ATOM_TYPE"][first_atom-1+i],
#           prmtop.parm_data["CHARGE"][first_atom-1+i], chrgdat[first_charge+state*num_atoms+i])
      prmtop.parm_data["CHARGE"][first_atom-1+i] = chrgdat[first_charge+state*num_atoms+i]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def mdEnergy(mden):
   file = open(mden,'r')
   lines = file.readlines()
   file.close()

   return float(lines[16].split()[2])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def minEnergy(mdout):
   file = open(mdout,'r')
   lines = file.readlines()
   file.close()

   energy = 0
   for i in range(len(lines)):
      if lines[i].strip().startswith("BOND"):
         words = lines[i].split()
         energy += float(words[2]) + float(words[5]) + float(words[8])
         words = lines[i+1].split()
         energy += float(words[2]) + float(words[5]) + float(words[8])
         words = lines[i+2].split()
         energy += float(words[3]) + float(words[7]) + float(words[10])
         return energy

   return -1
         
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# load default values
total_time = 10
dt         = 10
pH         = 7
cpin       = 'cpin'
prmtop     = 'prmtop'
inpcrd     = 'inpcrd'
mdout      = 'mdout'
cpout      = 'cpout'
restrt     = 'restrt'
mdcrd      = 'mdcrd'

# answer calls for help
if len(sys.argv) > 1:
   if sys.argv[1] == '-h' or sys.argv[1] == '--help':
      usage()

# parser for command-line arguments
i = 1
try:
   while i < len(sys.argv):
      if sys.argv[i] == '-t':
         i += 1
         total_time = int(sys.argv[i])
      elif sys.argv[i] == '-dt':
         i += 1
         dt = int(sys.argv[i])
      elif sys.argv[i] == '-pH' or sys.argv[i] == '-ph':
         i += 1
         pH = float(sys.argv[i])
      elif sys.argv[i] == '-cpin':
         i += 1
         cpin = sys.argv[i]
      elif sys.argv[i] == '-p':
         i += 1
         prmtop = sys.argv[i]
      elif sys.argv[i] == '-c':
         i += 1
         inpcrd = sys.argv[i]
      elif sys.argv[i] == '-o':
         i += 1
         mdout = sys.argv[i]
      elif sys.argv[i] == '-cpout':
         i += 1
         cpout = sys.argv[i]
      elif sys.argv[i] == '-r':
         i += 1
         restrt = sys.argv[i]
      elif sys.argv[i] == '-x':
         i += 1
         mdcrd = sys.argv[i]
      else:
         print >> sys.stderr, 'Error: Unrecognized input option %s!' % sys.argv[i]
         usage()
      i += 1

except IndexError:
   print >> sys.stderr, 'Error: Expecting more command line arguments!'
   usage()
except ValueError:
   print >> sys.stderr, 'Error: Unexpected data types obtained (i.e. float, int, etc.)'
   usage()

FACTOR *= pH # put the pH in the factor

# Variable declaration
pmemd_input = mdin('pmemd') # we're creating a pmemd input file
pmemd_energy = mdin('pmemd') # we're creating another pmemd input file
parm = amberParm(prmtop)
resstates = []
chrgdat = []
protcnt = []
num_states = []
first_state = []
first_charge = []
first_atom = []
num_atoms = []
statene = []
trescnt = 0

# Check for errors and parse CPIN file.
if not parm.exists:
   print >> sys.stderr, 'Error: Prmtop %s does not exist!' % prmtop
   usage()

try:
   cpin_file = open(cpin, 'r')
except IOError:
   print >> sys.stderr, 'Error: Cpin %s does not exist!' % cpin
   usage()

try:
   test = open(inpcrd, 'r')
except IOError:
   print >> sys.stderr, 'Error: Inpcrd %s does not exist!' % inpcrd
   usage()

trescnt = parseCpin(cpin_file, chrgdat, resstates, num_states, first_state, first_charge, first_atom, num_atoms, statene, protcnt)

if trescnt == -1:
   print >> sys.stderr, 'Error: Cannot parse cpin file %s!' % cpin
   sys.exit()

cpin_file.close()

# Set initial states
for i in range(trescnt):
   updateParm(parm, first_charge[i], chrgdat, resstates[i], num_atoms[i], first_atom[i])

parm.writeParm('current.%s' % prmtop)

# Get first random state
randres = random.randint(0,trescnt-1)
randstate = random.randint(0,num_states[randres]-2)

if randstate >= resstates[randres]:
   randstate += 1

updateParm(parm, first_charge[randres], chrgdat, randstate, num_atoms[randres], first_atom[randres])

parm.writeParm('proposed.%s' % prmtop)

# write the mdin file
pmemd_input.constTemp()
pmemd_input.constPressure()
pmemd_input.restart()
pmemd_input.SHAKE()
pmemd_input.change('cntrl','dt',0.002)
pmemd_input.change('cntrl','nstlim',dt)
pmemd_input.change('cntrl','ioutfm',1)
pmemd_input.change('cntrl','ntwx',dt)
pmemd_input.change('cntrl','ntpr',dt)
pmemd_input.change('cntrl','ntwe',dt)

pmemd_input.write('pmemd_pH.mdin')

pmemd_energy.minimization(maxcyc=1)
pmemd_energy.change('cntrl','ntpr',1)
pmemd_energy.write('pmemd_min.mdin')

# Now it's time to enter the actual MD loop
numsteps = int(total_time / dt)

try:
   mpi_cmd = os.environ['DO_PARALLEL']
except KeyError:
   mpi_cmd = 'none'

if mpi_cmd == 'none':
   pmemd = which('pmemd')
   exe = '%s' % pmemd
else:
   pmemd = which('pmemd.MPI')
   exe = '%s %s' % (mpi_cmd, pmemd)

if pmemd == 'none':
   print >> sys.stderr, 'Error: PMEMD cannot be found! It must be in your PATH for cph_pmemd.py to work.'
   sys.exit()

cpout_file = open(cpout, 'w')
cpout_file.write("""Solvent pH: %8.5f
Monte Carlo step size: %8i
Time step: %8i
Time: %10.3f
""" % (pH, dt, 0, 0))

for i in range(trescnt):
   cpout_file.write('Residue %4i State: %2i\n\n' % (i, resstates[i]))

for i in range(numsteps):
   if os.system('%s -i pmemd_pH.mdin -p %s -c %s -o %s -r %s -e mden' %  (exe, 'current.'+prmtop, inpcrd, mdout, restrt)) != 0:
      print >> sys.stderr, 'Error: %s did not exit cleanly!' % exe
      sys.exit()
   os.system('mv %s %s' % (restrt, inpcrd))

   if os.system('%s -O -i pmemd_min.mdin -p %s -c %s -r %s -o _min.mdout' % (exe, 'proposed.'+prmtop, inpcrd, restrt)) != 0:
      print >> sys.stderr, 'Error: %s did not exit cleanly!' % exe
      sys.exit()

   e_cur = mdEnergy('mden')
   e_pro = minEnergy('_min.mdout')

   print 'Energies : %12.5f, %12.5f' % (e_cur, e_pro)

   if monteCarlo(e_pro - e_cur, statene[first_state[randres]+resstates[randres]]-statene[first_state[randres]+randstate],
                 protcnt[first_state[randres]+resstates[randres]], protcnt[first_state[randres]+randstate]): # accept change
      cpout_file.write('Residue %4i State: %2i\n\n' % (randres, randstate))
      resstates[randres] = randstate
      os.system('mv proposed.%s current.%s' % (prmtop, prmtop))
      print 'Accepted!'
   else: # reject the change
      cpout_file.write('Residue %4i State: %2i\n\n' % (randres, resstates[randres]))
      updateParm(parm, first_charge[randres], chrgdat, resstates[randres], num_atoms[randres], first_atom[randres])
      print 'Rejected!'


