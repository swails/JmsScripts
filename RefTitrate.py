#! /usr/bin/env python

""" 
This program will titrate model compounds and determine
how their energies should be adjusted in order to get
the correct populations
"""

# Should all be here
import sys, os, math
from optparse import OptionParser
from subprocess import Popen

from mpi4py import MPI
from utilities import which

# Set up MPI
commworld = MPI.COMM_WORLD
commrank = commworld.Get_rank()
commsize = commworld.Get_size()
master = commrank == 0

# Set up the parser and add options
parser = OptionParser()

parser.add_option('-r', '--residue', dest='res', help='Residue name to titrate')
parser.add_option('-p', '--pH', dest='pH', type='float', help='Target pH')
parser.add_option('-f', '--fine-resolution', dest='fineres', default=0.2, 
                  type='float', help='pH increments near pKa')
parser.add_option('-c', '--coarse-resolution', dest='coarseres', default=1.0,
                  type='float', help='pH increments further from pKa')
parser.add_option('-g', '--igb', dest='igb', default=5, type='int', 
                  help='GB value to parameterize for')
parser.add_option('-t', '--nstlim', dest='nstlim', default=2000000, type='int',
                  help='How long to run each window')
parser.add_option('-a', '--amino-acid', dest='aa', action='store_true', default=True)
parser.add_option('-n', '--nucleic-acid', dest='aa', action='store_false', default=True)

(options, args) = parser.parse_args()

# Make sure we have enough arguments
if options.res == None or options.pH == None:
   parser.print_help()
   sys.exit()

# Now determine where required programs are
sander = which('sander')
tleap = which('tleap')
cpinutil = which('cpinutil')
if options.igb == 8:
   converter = which('ChangeParmRadii.py')

if 'none' in [sander, tleap, cpinutil]:
   print >> sys.stderr, 'sander, tleap, and cpinutil are all necessary!'
   sys.exit()

if options.igb == 8 and converter == 'none':
   print >> sys.stderr, 'ChangeParmRadii.py is needed for igb = 8!'
   sys.exit()

# Keep a log of all stdout
log = open('%s.log' % os.path.split(sys.argv[0])[1].strip('.py'), 'w')

mdin = """Mdin file for titrating stuff
 &cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1000, nstlim=%s,
   dt=0.002, ntt=3, tempi=300,
   temp0=300, tautp=2.0, ig=-1,
   ntp=0, ntc=2, ntf=2, cut=999,
   ntb=0, igb=%s, saltcon=0.1,
   nrespa=1, tol=0.000001, icnstph=1,
   solvph=%%s, ntcnstph=5, gamma_ln=2.0,
   ntwr=500, ioutfm=1,
 /
""" % (options.nstlim, options.igb)

min_mdin = """Minimization to relax initial bad contacts, implicit solvent
 &cntrl
   imin=1,
   ncyc=100,
   maxcyc=200,
   ntpr=50,
   ntb=0,
   cut=1000,
 /
"""

if options.aa:
   caps = ['ACE', 'NME']
else:
   caps = ['MOC', 'CH3']

tleapin = """source leaprc.constph
l = sequence {%s %s %s}
saveamberparm l %s.parm7 %s.rst7
quit
""" % (caps[0], options.res, caps[1], options.res, options.res)

if master:
   # First it's time to create the prmtop
   file = open('tleap.in', 'w')
   file.write(tleapin)
   file.close()

   proc_return = Popen(['tleap', '-f', 'tleap.in'], executable=tleap, stdout=log)

   if proc_return != 0:
      print >> sys.stderr, 'tleap error!'
      sys.exit()
   
   # If we're doing igb = 8, do the prmtop conversion
   proc_return = Popen(['ChangeParmRadii.py','-p','%s.parm7' % options.res, '-r', 'mbondi3'],
                       executable=converter, stdout=log)

   # Create the cpin
   cpin = open(options.res + '.cpin', 'w')
   proc_return = Popen(['cpinutil', '-p', '%s.parm7' % options.res, '-igb', options.igb],
                       executable=cpinutil, stdout=cpin, stderr=log)
   cpin.close()

   if proc_return != 0:
      print >> sys.stderr, 'cpinutil error!'
      sys.exit()

   # Now it's time to minimize the structure
   mdin = open('mdin.min', 'w')
   mdin.write(min_mdin)
   mdin.close()
   
   proc_return = Popen(['sander', '-i', 'mdin.min', '-c', '%s.rst7' % options.res, '-p', 
                        '%s.parm7' % options.res,'-o', 'min.mdout', '-r', '%s.min.rst7' % options.res], 
                        executable=sander)
   
   if proc_return != 0:
      print >> sys.stderr, 'sander minimization error!'
      sys.exit()

   os.system('rm -f min.mdout leap.log')

# Now we're done with the system, so it's time to start each process

# First determine which pH values have to be simulated

starting_pH = math.floor(options.pH * 10) / 10
lower_pHs = starting_pH - options.fineres
higher_pHs = starting_pH + options.fineres
pH_sims = [lower_pHs, starting_pH, higher_pHs]

# Now go 2 above and 2 below the pH at both resolutions

if starting_pH < options.pH:
   lower_pHs -= options.coarseres
   pH_sims.append(lower_pHs)
   lower_pHs -= options.coarseres
   pH_sims.append(lower_pHs)
   higher_pHs += options.fineres
   pH_sims.append(higher_pHs)
   higher_pHs += options.coarseres
   pH_sims.append(higher_pHs)
   higher_pHs += options.coarseres
   pH_sims.append(higher_pHs)
else:
   lower_pHs -= options.coarseres
   pH_sims.append(lower_pHs)
   lower_pHs -= options.coarseres
   pH_sims.append(lower_pHs)
   lower_pHs -= options.coarseres
   pH_sims.append(lower_pHs)
   higher_pHs += options.fineres
   pH_sims.append(higher_pHs)
   higher_pHs += options.coarseres
   pH_sims.append(higher_pHs)

# Now split up the load
xtras = len(pH_sims) % commsize

num_frames = int(pH_sims / commsize)
start_frame = num_frames * rank

if rank < xtras:
   num_frames += 1
   start_frame += rank
else:
   start_frame += xtras

end_frame = start_frame + num_frames

# Now it's time to write the MDIN file
for i in range(startframe, endframe):
   mdin = open('mdin.%s' % commrank, 'w')
   mdin.write(mdin % pH_sims[i])
   mdin.close()

   proc_return = Popen(['sander', '-i', 'mdin.%s' % commrank, '-c', '%s.min.rst7' % options.res, '-p', 
                        '%s.parm7' % options.res,'-o', '%s_pH%s.mdout' % (options.res, pH_sims[i]), '-r', 
                        '%s.md.rst7' % options.res, '-inf', commrank + '.mdinfo', '-cpin', options.res + '.cpin',
                        '-cpout', '%s_pH%s.cpout' % (options.res, pH_sims[i]), '-cprestrt', commrank + '.cprestrt'], 
                        executable=sander)
