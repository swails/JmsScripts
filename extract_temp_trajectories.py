#!/usr/bin/env python

""" 
This program will take a list of REMD trajectory files and temperatures and it
will extract temperature-specific trajectories from those trajectory files for
each temperature provided.
"""
from optparse import OptionParser, OptionGroup
from subprocess import Popen, PIPE
import sys, time

start_time = time.time()

def which(program):
   import os
   def is_exe(fpath):
      return os.path.exists(fpath) and os.access(fpath, os.X_OK)

   fpath, fname = os.path.split(program)
   if fpath:
      if is_exe(program):
         return program
   else:
      for path in os.environ["PATH"].split(os.pathsep):
         exe_file = os.path.join(path, program)
         if is_exe(exe_file):
            return exe_file
   return None

parser = OptionParser('usage: %prog [Options] mdcrd1 [mdcrd2 [mdcrd3 ...] ]')
parser.add_option('--temperatures', dest='temp_list', help='Comma-separated ' +
                  'list of temperatures to extract trajectories for')
parser.add_option('--prefix', dest='prefix', help='Prefix for trajectory ' +
                  'file names. They will be named PREFIX.temp.nc')
parser.add_option('--prmtop', dest='prmtop', help='Prmtop for your system')
group = OptionGroup(parser, 'RMSd Options', 'If a REFSTRUCT is specified, ' +
                    'these options will be used to RMS fit a structure')
group.add_option('--rmsd', dest='refstruct', help='Reference structure for ' +
                 'RMS fitting (if desired)')
group.add_option('--rmsd-mask', dest='rmsmask', default='@CA,C,O,N,H',
                 help='Mask to RMS fit around. Default=@CA,C,O,N,H')
group.add_option('--rmsd-output', dest='rmsout', default='rms',
                 help='Output file name prefix for RMSd data dump. ' +
                 'Default rms.out')
parser.add_option_group(group)
opt,trajs = parser.parse_args()

if not opt.temp_list or not opt.prefix or not opt.prmtop or not trajs:
   parser.print_help()
   sys.exit(1)

cpptraj = which('cpptraj')

assert(cpptraj)

temps = opt.temp_list.split(',')

for temp in temps:
   try: temp=float(temp)
   except ValueError: continue

   cpptraj_str = ''
   for traj in trajs:
      cpptraj_str += 'trajin %s remdtraj remdtrajtemp %s\n' % (traj, temp)

   if opt.refstruct: 
      cpptraj_str += 'reference %s\n' % opt.refstruct
      cpptraj_str += 'rmsd data %s reference out %s_%s.dat\n' % (opt.rmsmask,
                      opt.rmsout, temp)

   cpptraj_str += 'trajout %s.%f.nc netcdf' % (opt.prefix, temp)

   process = Popen([cpptraj, opt.prmtop], stdin=PIPE)

   process.communicate(cpptraj_str)

   process.wait()


print '\n\nThis took %f min.' % ((time.time() - start_time) / 60)
