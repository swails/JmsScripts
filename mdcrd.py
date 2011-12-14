"""
This module contains code useful for analyzing Amber trajectory files and
generating useful data sets for analyzing them.
"""

import numpy as np
from subprocess import Popen, PIPE
from utilities import which
import re, sys, os

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class TrajError(Exception):
   def __init__(self, msg='Trajectory error!'):
      self.msg = msg
   def __str__(self):
      return self.msg

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class RMSError(Exception): pass

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class AmberTraj(object):
   " This is a class to analyze trajectory files (a series of them, or just 1) "
   
   def __init__(self, parm, traj_list, start=1, stride=1, end=99999999,
                logfile=None, overwrite=False):
      # Get cpptraj
      self.cpptraj = which('cpptraj')
      if not self.cpptraj:
         raise TrajError('Could not find cpptraj!')

      # Load instance data
      self.parm = str(parm)

      # Order matters
      self.traj_name_list = []

      # Start, end, stride
      self.start, self.stride, self.end = [], [], []

      # traj_list is a dictionary that matches the trajectory name to the
      # number of frames in that trajectory
      self.traj_list = {}
      if type(traj_list).__name__ in ['list', 'tuple']:
         for item in traj_list:
            self.traj_list[item] = -1
            self.traj_name_list.append(item)
            self.start.append(start)
            self.stride.append(stride)
            self.end.append(end)
      elif type(traj_list).__name__ == 'str':
         self.traj_list[traj_list] = -1
         self.traj_name_list.append(traj_list)
         self.start.append(start)
         self.stride.append(stride)
         self.end.append(end)
      else:
         raise TypeError('AmberTraj: trajectory(s) must be a list/tuple or str')

      # Now set up num_frames
      self._query()

      # Can we overwrite files?
      self.overwrite = overwrite

      # what is our logfile?
      if type(logfile).__name__ == 'str':
         if os.path.exists(logfile) and not self.overwrite:
            raise TrajError('Cannot overwrite %s' % logfile)
         self.logfile = open(logfile, 'w')
      elif not logfile:
         self.logfile = sys.stdout
      else:
         self.logfile = logfile

      # Start keeping track of the commands we want to run
      self._cpptraj_commands = ''

   #===================================================

   def add_trajectory(self, name, start=1, stride=1, end=99999999):
      """ Adds a trajectory to the list """
      if type(name).__name__ == 'str':
         self.traj_list[name] = -1
         self.traj_name_list.append(name)
         self.start.append(start)
         self.stride.append(stride)
         self.end.append(end)
      elif type(name).__name__ in ['list', 'tuple']:
         for i in name:
            self.traj_list[i] = -1
            self.traj_name_list.append(i)
            self.start.append(start)
            self.stride.append(stride)
            self.end.append(end)
      self._query()

   #===================================================

   def _query(self):
      " Determine how many frames are in each trajectory "
      get_num_frames = re.compile(r' *\[.+\] contains (\d+) frames.')
      for traj in self.traj_list:
         # Skip over any frames we've already determined
         if self.traj_list[traj] != -1: continue
         # Now launch a subprocess where we just trajin the file and parse
         # the output to find out how many frames are present
         process = Popen([self.cpptraj, self.parm], stdin=PIPE, stdout=PIPE,
                         stderr=PIPE)
         out, err = process.communicate('trajin %s' % traj)
         # This shouldn't be reached, since cpptraj doesn't bail out with a
         # non-zero exit code, it just prints errors.
         if process.wait():
            raise TrajError('Bad trajectory file %s:\nOutput: %s\nError: %s' % 
                            (traj, out, err))
         nframes = get_num_frames.search(out)
         # This really means we put in a bad trajectory file
         if len(nframes.groups()) != 1:
            raise TrajError('Bad trajectory file (%s):\nOutput: %s\nError: %s' %
                            (traj, out, err))
         # Otherwise, we got our number of frames
         self.traj_list[traj] = int(nframes.groups()[0])

   #===================================================

   def rmsd(self, setname='rmsd', mask=':*', ref=None, outfile=None, fit=True):
      """ Calculates the RMSD value for all of the trajectories """
      # Make sure we aren't overwriting something:
      if outfile and os.path.exists(outfile) and not self.overwrite:
         raise TrajError('Cannot overwrite %s' % outfile)

      # See if we RMSD to a reference
      if ref:
         self._cpptraj_commands += 'reference %s\n' % ref
      # Now start the RMS commands
      self._cpptraj_commands += 'rmsd %s %s ' % (setname, mask)

      # Do we have a reference structure? If not, use first
      if ref:
         self._cpptraj_commands += 'reference ref %s ' % ref
      else:
         self._cpptraj_commands += 'first '

      # Dump to a file?
      if outfile:
         self._cpptraj_commands += 'out %s ' % outfile

      # Terminate the command
      self._cpptraj_commands += '\n'

   #===================================================

   def run(self):
      """ This runs cpptraj with the given commands """
      # adjust the ends to be either the highest frame # or the value of end
      self.end = [min(self.end[i], self.traj_list[j]) 
                  for i,j in enumerate(self.traj_name_list)]
      cmd_str = ''
      for i, traj in enumerate(self.traj_name_list):
         cmd_str += 'trajin %s %d %d %d \n' % (traj, self.start[i], self.end[i],
                                               self.stride[i])
      
      process = Popen([self.cpptraj, self.parm], stdin=PIPE)
#                     stdout=self.logfile, stderr=self.logfile)
      
      print >> self.logfile, 'Running cpptraj:'
      process.communicate(cmd_str + self._cpptraj_commands)

      if process.wait():
         print >> self.logfile, 'Running cpptraj failed.'
      else:
         print >> self.logfile, 'cpptraj ran successfully!'

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

class RmsdData(object):
   """ 
   Collects the RMSD data from a cpptraj-generated data file and loads it into
   numpy arrays
   """
   def __init__(self, fname, col=1):
      from utilities import linecount
      if not os.path.exists(fname):
         raise RMSError('RmsdData: Non-existent file %s' % fname)
      data_re = re.compile(r'\s*(\d+)')
      nums_re = re.compile(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
      # Find the number of lines in the 
      numframes = linecount(fname)
      # Load the numpy array
      self.data = np.empty(numframes)
      idx = 0
      for line in open(fname, 'r'):
         if not data_re.match(line): continue
         # The header line doesn't match
         cols = nums_re.findall(line)
         # cols is every number in the line, separated by column. The first
         # number is the frame number
         self.data[idx] = float(cols[col])
         idx += 1

#~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~

if __name__ == '__main__':
   # Test suite
   from optparse import OptionParser
   import sys

   parser = OptionParser('python %prog [Options] mdcrd1 [mdcrd2 [...] ]')
   parser.add_option('-p', '--parm', dest='prmtop', metavar='FILE',
                   default=None, help='Topology file matching the trajectories')
   opt, args = parser.parse_args()

   if not args or not opt.prmtop:
      print 'Missing information!'
      parser.print_help()
      sys.exit(1)

   mytraj = AmberTraj(opt.prmtop, args, logfile='mdcrd_py.log', overwrite=True)

   # Print out the sizes of all of the trajectories
   for traj in mytraj.traj_name_list:
      print 'Trajectory %20s has %10d frames.' % (traj, mytraj.traj_list[traj])

   # separator
   print ''

   # Test the RMSd
   mytraj.rmsd(outfile='AmberTraj_RMSD.dat')

   mytraj.run()

   # Now test the RMSd class

   my_rmsd = RmsdData('AmberTraj_RMSD.dat.check')

   # Print out some statistics about the RMSD data
   print 'The average RMSD is ', np.average(my_rmsd.data)
   print 'The minimum RMSD is ', my_rmsd.data.min()
   print 'The maximum RMSD is ', my_rmsd.data.max()
   print 'The stdev of the RMSD is ', my_rmsd.data.std()

