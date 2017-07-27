#! /usr/bin/env python

""" 
This program will titrate model compounds and determine
how their energies should be adjusted in order to get
the correct populations
"""

from __future__ import division, print_function

# Should all be here
from optparse import OptionParser, OptionGroup
import os
from subprocess import Popen
import sys
from utilities import which

import parmed as pmd

# Set up the parser and add options
parser = OptionParser()

group = OptionGroup(parser, 'Residue Info', 'Information about the residue you '
                    'wish to titrate')
group.add_option('-r', '--residue', dest='res', metavar='RESIDUE_NAME',
                 help='Residue name to titrate')
group.add_option('-p', '--pKa', dest='pKa', type='float', metavar='FLOAT',
                 help='pKa of the model compound. No default', default=None)
group.add_option('-a', '--amino-acid', dest='aa', action='store_true',
                default=True, help='Use amino acid caps for reference compound')
group.add_option('-u', '--nucleic-acid', dest='aa', action='store_false',
              default=True, help='Use nucleic acid caps for reference compound')
parser.add_option_group(group)

group = OptionGroup(parser, 'Advanced Residue Options', 'If you wish to '
                    'titrate a non-standard residue, you may need to specify '
                    'the residue sequence of the reference compound as well '
                    'as the residue name itself. These options allow you to '
                    'specify custom parameters and capping residues')
group.add_option('-l', '--left-residue', dest='leftres', default=None,
                 help='Which residue to cap with on the left terminus')
group.add_option('-i', '--right-residue', dest='rightres', default=None,
                 help='Which residue to cap with on the right terminus')
group.add_option('--frcmod', dest='frcmod', metavar='FILE', default=None,
                 help='File with additional parameters for the compound')
group.add_option('--lib', dest='lib', metavar='FILE', default=None,
                 help='File with the residue definition (OFF file)')
parser.add_option_group(group)

group = OptionGroup(parser, 'Titration Options', 'Options controlling the '
                    'titration of the reference compounds')
group.add_option('-f', '--fine-resolution', dest='fineres', default=0.2,
                 type='float', metavar='FLOAT', help='pH increments near pKa')
group.add_option('-c', '--coarse-resolution', dest='coarseres', default=1.0,
                 type='float', metavar='FLOAT',
                 help='pH increments further from pKa')

parser.add_option_group(group)

group = OptionGroup(parser, 'Simulation Options', 'Options controlling the '
                    'simulation parameters.')
group.add_option('-g', '--igb', dest='igb', default=5, type='int', 
                 help='GB value to parameterize for')
group.add_option('-t', '--nstlim', dest='nstlim', default=2000000, type='int',
                 help='How long to run each window')
group.add_option('-e', '--eaf', dest='eaf', default=100, type='int',
                 metavar='INT', help='Number of steps between replica exchange '
                 'attempts.')
group.add_option('-n', '--num-replicas', default=6, metavar='INT', dest='nreps',
                 help='Number of replicas to run in REMD calcs (must be even!) '
                 'Default %default', type='int')
group.add_option('-m', '--mpi-cmd', dest='mpi_cmd', metavar='STRING',
                 help='MPI Command to run MPI programs on your machine. ('
                 'Default "%default")', default='mpiexec -n 6')
group.add_option('-d', '--intdiel', dest='intdiel', type='float', default=1.0,
                 metavar='FLOAT', help='Internal dielectric to use. Default is '
                 '1.0')
parser.add_option_group(group)

(options, args) = parser.parse_args()

# Make sure we have enough arguments
if options.res is None or options.pKa is None:
   parser.print_help()
   sys.exit(1)

# Now determine where required programs are
sander = which('sander')
sanderMPI = which('sander.MPI')
tleap = which('tleap')
cpinutil = which('cpinutil.py')
converter = which('parmed')

if options.nreps % 2 != 0:
   print('Error: Even number of replicas required!', file=sys.stderr)
   sys.exit(1)

if None in [sander, sanderMPI, tleap, cpinutil]:
   print('sander, tleap, and cpinutil.py are all necessary!', file=sys.stderr)
   sys.exit(1)

if options.igb == 8 and converter is None:
   print('parmed is needed for igb = 8!', file=sys.stderr)
   sys.exit(1)

print(" Found necessary programs!")

# Keep a log of all stdout
log = open('%s.log' % os.path.split(sys.argv[0])[1].strip('.py'), 'w')

nstlim = options.eaf
numexchg = options.nstlim // options.eaf

md_mdin = """Mdin file for titrating stuff
 &cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1000, nstlim=%d,
   dt=0.002, ntt=3, tempi=300,
   temp0=300, tautp=2.0, ig=-1,
   ntp=0, ntc=2, ntf=2, cut=999,
   ntb=0, igb=%s, saltcon=0.1,
   nrespa=1, tol=0.000001, icnstph=1,
   solvph=%%s, ntcnstph=5, gamma_ln=2.0,
   ntwr=500, ioutfm=1, numexchg=%d,
   intdiel=%s,
 /
""" % (nstlim, options.igb, numexchg, options.intdiel)

min_mdin = """Minimization to relax initial bad contacts, implicit solvent
 &cntrl
   imin=1,
   ncyc=100,
   maxcyc=200,
   ntpr=50,
   ntb=0,
   cut=1000,
   igb=%d,
 /
""" % options.igb

if options.leftres is not None:
   left_term = options.leftres
elif options.aa:
   left_term = 'ACE'
else:
   left_term = 'MOC'

if options.rightres is not None:
   right_term = options.rightres
elif options.aa:
   right_term = 'NME'
else:
   right_term = 'CH3'

tleapin = "source leaprc.constph\n"

if options.frcmod is not None:
   tleapin += "loadamberparams %s\n" % options.frcmod
if options.lib is not None:
   tleapin += "loadoff %s\n" % options.lib

tleapin += """l = sequence {%s %s %s}
set default pbradii %%s
saveamberparm l %s.parm7 %s.rst7
quit
""" % (left_term, options.res, right_term, options.res, options.res)

if options.igb == 8:
   tleapin = tleapin % 'mbondi3'
elif options.igb == 7:
   tleapin = tleapin % 'bondi'
elif options.igb == 1:
   tleapin = tleapin % 'mbondi'
else:
   tleapin = tleapin % 'mbondi2'

f = open('tleap.in', 'w')
f.write(tleapin)
f.close()

# First it's time to create the prmtop
print("\n Making topology file")
file = open('tleap.in', 'w')
file.write(tleapin)
file.close()

proc_return = Popen(['tleap', '-f', 'tleap.in'], executable=tleap, stdout=log).wait()

if proc_return != 0:
   print('tleap error!', file=sys.stderr)
   sys.exit(1)

print(" Successfully created topology file %s.parm7" % options.res)

#print(' Changing the carboxylate radii to 1.3 A')
#parm = pmd.load_file('%s.parm7' % options.res)
#pmd.tools.change(parm, 'RADII', '@OD=,OE=', 1.3).execute()
#parm.write_parm('%s.parm7' % options.res)

# Create the cpin
print("\n Creating cpin file")
cpin = open(options.res + '.cpin', 'w')
proc_return = Popen([cpinutil, '-p', '%s.parm7' % options.res, '-igb', 
                     str(options.igb), '-intdiel', str(options.intdiel)],
                     stdout=cpin, stderr=log).wait()
cpin.close()
print(" Finished making cpin file")

if proc_return != 0:
   print('cpinutil error!', file=sys.stderr)
   sys.exit(1)

# Now it's time to minimize the structure
mdin = open('mdin.min', 'w')
mdin.write(min_mdin)
mdin.close()

print("\n Minimizing initial structure")
proc_return = Popen([sander, '-O', '-i', 'mdin.min', '-c',
                   '%s.rst7' % options.res, '-p', '%s.parm7' % options.res,
                   '-o', 'min.mdout', '-r', '%s.min.rst7' % options.res]).wait()

if proc_return != 0:
   print('sander minimization error!', file=sys.stderr)
   sys.exit(1)

print(" Structure minimized")

os.unlink('min.mdout')
os.unlink('leap.log')

# Now we're done with the system, so it's time to start each process

# First determine which pH values have to be simulated
half_nreps = options.nreps // 2
ncoarse = half_nreps // 2
nfine = half_nreps - ncoarse

offsets = [options.fineres]
for i in range(1, nfine):
   offsets.append(options.fineres * (i + 1))

for i in range(ncoarse):
   offsets.append(options.fineres * nfine + options.coarseres * (i + 1))

offlist = [-i for i in offsets]
offlist.reverse()
offlist.extend(offsets)

phlist = [options.pKa + i for i in offlist]
del offlist, offsets

# open up the groupfile for writing
grpfile = open('groupfile', 'w')
opts = {'sys' : options.res}
# Write out all of the MDIN files and each line of the groupfile
for i, ph in enumerate(phlist):
   mdin = open('mdin.%d' % i, 'w')
   mdin.write(md_mdin % ph)
   mdin.close()
   opts['num'] = i; opts['ph'] = ph
   grpfile.write(('-O -i mdin.%(num)s -c %(sys)s.min.rst7 -p %(sys)s.parm7 '
                  '-o %(sys)s_pH%(ph)s.mdout -r %(sys)s.md.rst7.%(num)s '
                  '-inf %(num)s.mdinfo -cpin %(sys)s.cpin -cpout '
                  '%(sys)s_pH%(ph)s.cpout -cprestrt %(num)d.cprestrt -rem 4 '
                  '-remlog %(sys)s_rem.log\n') % opts)

grpfile.close()
   
# Now run the simulation!
print('Beginning titration of %d replicas...' % options.nreps)
print('\tpKa is %f' % options.pKa)
print('\tSimulating pH values ' + ', '.join([str(i) for i in phlist]))
if os.system('%s %s -ng %d -groupfile groupfile' % (options.mpi_cmd, sanderMPI,
                                                    options.nreps)):
   print('Error during calculation!')
   sys.exit(1)

print('Done! Don\'t forget to process output as pH-REMD')
