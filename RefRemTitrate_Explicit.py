#! /usr/bin/env python

""" 
This program will titrate model compounds and determine
how their energies should be adjusted in order to get
the correct populations
"""

from __future__ import division

# Should all be here
from optparse import OptionParser, OptionGroup
import os
from subprocess import Popen
import sys
from utilities import which

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
group.add_option('-b', '--box-size', type='float', dest='box', default=10.0,
                 help='Size of the solvent box. Default %default Angstroms')
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
group.add_option('--ntcnstph', metavar='INT', dest='ntcnstph', default=50,
                 type='int', help='Number of steps between protonation state '
                 'change attempts. (Default %default)')
parser.add_option_group(group)

group = OptionGroup(parser, 'Simulation Options', 'Options controlling the '
                    'simulation parameters.')
group.add_option('-g', '--igb', dest='igb', default=5, type='int', 
                 help='GB value to parameterize for')
group.add_option('--ntrelax', dest='ntrelax', type='int', metavar='INT',
                 default=50,
                 help='Water relaxation steps to run (Default %default)')
group.add_option('-d', '--intdiel', dest='dielc', type='float', metavar='FLOAT',
                 help='Internal dielectric constant to use for protonation '
                 'state change evaluations. Default 1.', default=1.0)
group.add_option('-t', '--nstlim', dest='nstlim', default=2000000, type='int',
                 help='How long to run each window')
group.add_option('--heating-steps', dest='heatsteps', metavar='INT', type='int',
                 default=100000, help='How many heating steps to run (Default '
                 '%default)')
group.add_option('--equil-steps', dest='equisteps', metavar='INT', type='int',
                 default=1000000, help='How many equilibration steps to run '
                 '(Default %default)')
group.add_option('-e', '--eaf', dest='eaf', default=100, type='int',
                 metavar='INT', help='Number of steps between replica exchange '
                 'attempts.')
group.add_option('-n', '--num-replicas', default=6, metavar='INT', dest='nreps',
                 help='Number of replicas to run in REMD calcs (must be even!) '
                 'Default %default', type='int')
group.add_option('-m', '--mpi-cmd', dest='mpi_cmd', metavar='STRING',
                 help='MPI Command to run MPI programs on your machine. ('
                 'Default "%default")', default='mpiexec -n 6')
parser.add_option_group(group)

(options, args) = parser.parse_args()

# Make sure we have enough arguments
if options.res is None or options.pKa is None:
   parser.print_help()
   sys.exit(1)

if options.dielc != 1 and options.dielc != 2:
   sys.exit('--intdiel must be 1 or 2!')
# Now determine where required programs are
sander = which('sander.MPI')
pmemd = which('pmemd.MPI')
pmemd_cuda = which('pmemd.cuda')
tleap = which('tleap')
cpinutil = which('cpinutil.py')
converter = which('parmed.py')

if options.nreps % 2 != 0:
   print >> sys.stderr, 'Error: Even number of replicas required!'
   sys.exit(1)

if 'none' in [sander, tleap, cpinutil]:
   print >> sys.stderr, 'sander, tleap, and cpinutil.py are all necessary!'
   sys.exit(1)

if options.igb == 8 and converter == 'none':
   print >> sys.stderr, 'parmed.py is needed for igb = 8!'
   sys.exit(1)

print " Found necessary programs!"

# Keep a log of all stdout
log = open('%s.log' % os.path.split(sys.argv[0])[1].strip('.py'), 'w')

nstlim = options.eaf
numexchg = options.nstlim // options.eaf

md_mdin = """Mdin file for titrating stuff
 &cntrl
   imin=0, irest=1, ntx=5,
   ntpr=1000, nstlim=%d,
   dt=0.002, ntt=3, tempi=300,
   temp0=300, tautp=2.0, ig=-1,
   ntp=0, ntc=2, ntf=2, cut=8.0,
   saltcon=0.1,
   nrespa=1, tol=0.000001, icnstph=2,
   solvph=%%f, ntcnstph=%d, gamma_ln=2.0,
   ntwr=500, ioutfm=1, numexchg=%d,
   ntrelax=%d, ntwx=1000,
 /
""" % (nstlim, options.ntcnstph, numexchg, options.ntrelax)

min_mdin = """Minimization to relax initial bad contacts, explicit solvent
 &cntrl
   imin=1,
   ncyc=100,
   maxcyc=1000,
   ntpr=50,
   cut=8,
 /
"""

heat_mdin = """Slow heating in explicit solvent
 &cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1000, ntwx=1000, nstlim=%s,
   dt=0.002, ntt=3, gamma_ln=5.0, ig=-1,
   ntc=2, ntf=2, cut=8, ntb=2, ntp=1,
   iwrap=1, ioutfm=1, nmropt=1,
 /
 &wt
   TYPE='TEMP0', ISTEP1=0, ISTEP2=100000,
   VALUE1=50.0, VALUE2=300.0,
 /
 &wt TYPE='END' /
""" % (options.heatsteps)

equil_mdin = """Constant pressure equilibration dynamics
 &cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1000, ntwx=1000, nstlim=%s,
   dt=0.002, ntt=3, tempi=300,
   temp0=300, gamma_ln=1.0, ig=-1,
   ntp=1, ntc=2, ntf=2, cut=8,
   ntb=2, iwrap=1, ioutfm=1,
 /
""" % (options.equisteps)

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
solvateoct l TIP3PBOX %s
saveamberparm l %s.parm7 %s.rst7
quit
""" % (left_term, options.res, right_term, options.box, options.res,
       options.res)

f = open('tleap.in', 'w')
f.write(tleapin)
f.close()

# First it's time to create the prmtop
print "\n Making topology file"
file = open('tleap.in', 'w')
file.write(tleapin)
file.close()

proc_return = Popen([tleap, '-f', 'tleap.in'], stdout=log).wait()

if proc_return != 0:
   print >> sys.stderr, 'tleap error!'
   sys.exit(1)

print " Successfully created solvated topology file %s.parm7" % options.res

print "\n Setting prmtop radii"
# If we're doing igb = 8, do the prmtop conversion
if options.igb == 8:
   file = open('__TMP__','w')
   file.write('changeradii mbondi3\nsetoverwrite True\nparmout %s.parm7\n'
              % options.res)
   file.close()
   file = open('__TMP__','r')
   proc_return = Popen([converter,'%s.parm7' % options.res], 
                       stdout=log, stdin=file).wait()
   file.close()
   os.remove('__TMP__')
   print " Set prmtop radii to mbondi3"
else:
   file = open('__TMP__','w')
   file.write('changeradii mbondi2\nsetoverwrite True\nparmout %s.parm7\n'
              % options.res)
   file.close()
   file = open('__TMP__','r')
   proc_return = Popen([converter,'%s.parm7' % options.res],
                       stdout=log, stdin=file).wait()
   file.close()
   os.remove('__TMP__')
   print " Set prmtop radii to mbondi2"

# Create the cpin
print "\n Creating cpin file"
cpin = open(options.res + '.cpin', 'w')
proc_return = Popen([cpinutil, '-p', '%s.parm7' % options.res, '-igb', 
                     '%d' % options.igb, '-intdiel', '%s' % options.dielc],
                     stdout=cpin, stderr=log).wait()
cpin.close()
print " Finished making cpin file"

if proc_return != 0:
   print >> sys.stderr, 'cpinutil error!'
   sys.exit(1)

# Now it's time to minimize the structure
mdin = open('mdin.min', 'w')
mdin.write(min_mdin)
mdin.close()
mdin = open('heat.mdin', 'w')
mdin.write(heat_mdin)
mdin.close()
mdin = open('equil.mdin', 'w')
mdin.write(equil_mdin)
mdin.close()

print "\n Minimizing initial structure"
proc_return = os.system(options.mpi_cmd + ' ' + pmemd + ' ' +
                        ('-O -i mdin.min -c %s.rst7 -p %s.parm7 -o min.mdout '
                        '-r %s.min.rst7') % (options.res, options.res,
                        options.res))

if proc_return != 0:
   print >> sys.stderr, 'pmemd minimization error!'
   sys.exit(1)

print " Structure minimized"

os.unlink('min.mdout')
os.unlink('leap.log')

print '\n Heating structure'
if pmemd_cuda is not None:
   proc_return = os.system(pmemd_cuda + ' ' +
                        ('-O -i heat.mdin -c %s.min.rst7 -p %s.parm7 '
                        '-o heat.mdout -r %s.heat.rst7') % (options.res,
                        options.res, options.res))
else:
   proc_return = os.system(options.mpi_cmd + ' ' + pmemd + ' ' + 
                        ('-O -i heat.mdin -c %s.min.rst7 -p %s.parm7 '
                        '-o heat.mdout -r %s.heat.rst7') % (options.res,
                        options.res, options.res))

if proc_return != 0:
   print >> sys.stderr, 'pmemd heating error!'
   sys.exit(1)

print '\n Equilibrating structure'
if pmemd_cuda is not None:
   proc_return = os.system(pmemd_cuda + ' ' + 
                        ('-O -i equil.mdin -c %s.heat.rst7 -p %s.parm7 '
                        '-o equil.mdout -r %s.equil.rst7') % (options.res,
                        options.res, options.res))
else:
   proc_return = os.system(options.mpi_cmd + ' ' + pmemd + ' ' + 
                        ('-O -i equil.mdin -c %s.heat.rst7 -p %s.parm7 '
                        '-o equil.mdout -r %s.equil.rst7') % (options.res,
                        options.res, options.res))

if proc_return != 0:
   print >> sys.stderr, 'pmemd equilibration error!'
   sys.exit(1)

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
   grpfile.write(('-O -i mdin.%(num)s -c %(sys)s.equil.rst7 -p %(sys)s.parm7 '
                  '-o %(sys)s_pH%(ph)s.mdout -r %(sys)s.md.rst7.%(num)s '
                  '-inf %(num)s.mdinfo -cpin %(sys)s.cpin -cpout '
                  '%(sys)s_pH%(ph)s.cpout -cprestrt %(num)d.cprestrt -rem 4 '
                  '-remlog %(sys)s_rem.log -x %(sys)s_pH%(ph)s.nc\n') % opts)

grpfile.close()
   
# Now run the simulation!
print 'Beginning titration of %d replicas...' % options.nreps
print '\tpKa is %f' % options.pKa
print '\tSimulating pH values ' + ', '.join([str(i) for i in phlist])
if os.system('%s %s -ng %d -groupfile groupfile' % (options.mpi_cmd, sander,
                                                    options.nreps)):
   print 'Error during calculation!'
   sys.exit(1)

print 'Done! Don\'t forget to process output as pH-REMD'
