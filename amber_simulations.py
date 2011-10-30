#!/usr/bin/env python
""" 
This contains classes for a bunch of different calculation types
"""

from subprocess import Popen, PIPE
from chemistry.amber.readparm import AmberParm
from utilities import which
from mdin import mdin
from optparse import OptionParser
import os, sys

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def get_mpi_cmd():
   """ Returns the MPI command stored in ~/.mpiexec_cmd """
   if not os.path.exists(os.path.join(os.getenv('HOME'), '.mpiexec_cmd')):
      raise ProgramError('~/.mpiexec_cmd does not exist!')
   mpi_file = open(os.path.join(os.getenv('HOME'), '.mpiexec_cmd'), 'r')
   ret_str = mpi_file.read()
   mpi_file.close()
   return ret_str

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class InputError(Exception):
   """ Raised in the instance of bad input """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ProgramError(Exception):
   """ Raised in the instance of bad input """
   pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberSystem(object):
   """ An Amber system with a topology file and input coordinate file """

   # ==============================

   def __init__(self, prmtop, inpcrd):
      self.prmtop = AmberParm(prmtop, inpcrd)
      self.inpcrd = inpcrd

   # ==============================

   def periodic(self):
      """ Returns whether or not a topology file is set up for PBC """
      # PBCs are indicated by a non-zero IFBOX pointer
      return bool(self.prmtop.ptr('ifbox'))

   # ==============================

   def query_radii(self):
      """ Returns the radii set in the topology file """
      return self.prmtop.parm_data['RADIUS_SET'][0]

   # ==============================

   def set_radii(self, igb=0, radius_set='mbondi'):
      """ Sets the radius set to a new one using parmed """

      # If someone sets an igb, change to the appropriate radius set

      if igb:
         if igb == 1: radius_set = 'mbondi'
         elif igb == 2: radius_set = 'mbondi2'
         elif igb == 5: radius_set = 'mbondi2'
         elif igb == 7: radius_set = 'bondi'
         elif igb == 8: radius_set = 'mbondi3'

      if not radius_set in ['mbondi', 'mbondi2', 'mbondi3', 'bondi', 'amber6']:
         raise InputError('Bad radius set! Choose from ' +
                          'mbondi, mbondi2, mbondi3, bondi, and amber6')

      parmed = which('parmed.py')
      change_str = ("setOverwrite True\n" +
                    "changeRadii %s\n" % radius_set +
                    "parmout %s\n" % self.prmtop +
                    "go\n")

      process = Popen([parmed, '-q', '-n', str(self.prmtop)], stdin=PIPE,
                      stderr=PIPE, stdout=PIPE)

      (output, error) = process.communicate(change_str)

      if process.wait():
         raise ProgramError('parmed.py failed to change radii!')

      # Reload our topology file now that we've changed radius sets
      self.prmtop = AmberParm(str(self.prmtop))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BaseType(object):
   """ Base input file type """

   # ==============================

   def write_mdin(self, name):
      """ Writes the mdin file """
      self.mdin.write(name)

   # ==============================

   # Make changes to the minimization the same as the changes to the mdin file
   change = mdin.change
   
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Minimization(BaseType):
   """ Sets up a minimization """

   program = 'pmemd'

   # ==============================

   def __init__(self, amber_sys, num_steps=1000, igb=5, restrained=True,
                rst_wt=10.0, rst_mask='@CA,C,O,N'):
      """ Sets up a minimization input file """
      # Type checking
      if type(num_steps).__name__ != 'int':
         raise TypeError('num_steps must be an integer!')
      if type(igb).__name__ != 'int':
         raise TypeError('igb must be an integer!')
      # Create the mdin instance
      self.mdin = mdin(self.program) 
      # Set up periodic/non-periodic settings
      if amber_sys.periodic():
         self.mdin.change('cntrl', 'ntb', 1)
      else:
         self.mdin.change('cntrl', 'ntb', 0)
         if not igb in [1, 2, 5, 7, 8]:
            raise ValueError('For non-periodic systems, igb must be ' +
                             '1, 2, 5, 7, or 8')
         self.mdin.change('cntrl', 'igb', igb)
         self.mdin.change('cntrl', 'cut', 1000.0)
      self.mdin.change('cntrl', 'imin', 1)
      self.mdin.change('cntrl', 'maxcyc', num_steps)
      self.mdin.change('cntrl', 'ntr', int(restrained))
      if restrained:
         self.mdin.change('cntrl', 'restraint_wt', rst_wt)
         self.mdin.change('cntrl', 'restraintmask', rst_mask)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Heating(BaseType):
   """ Handles heating """
   
   program = 'pmemd'

   # ==============================

   def __init__(self, amber_sys, nstlim=1000, igb=5, restrained=True,
                rst_wt=2.0, rst_mask='@CA,C,O,N', temp0=300.0, tempi=10.0,
                slow=False, thermostat='langevin', thermostat_param=5.0,
                ntpr=1000, ntwx=1000, ntwr=10000):
      """ Sets up a heating input file """
      # Type checking
      if type(nstlim).__name__ != 'int':
         raise TypeError('nstlim must be an integer!')
      # Create the mdin instance
      self.mdin = mdin(self.program) 
      # Set up periodic/non-periodic settings
      if amber_sys.periodic():
         self.mdin.change('cntrl', 'ntb', 1)
      else:
         # Type-check igb
         if type(igb).__name__ != 'int':
            raise TypeError('igb must be an integer!')
         self.mdin.change('cntrl', 'ntb', 0)
         if not igb in [1, 2, 5, 7, 8]:
            raise ValueError('For non-periodic systems, igb must be ' +
                             '1, 2, 5, 7, or 8')
         self.mdin.change('cntrl', 'igb', igb)
         self.mdin.change('cntrl', 'cut', 1000.0)
      self.mdin.change('cntrl', 'ioutfm', 1) # always NetCDF!
      if not slow: self.mdin.change('cntrl', 'tempi', tempi)
      self.mdin.change('cntrl', 'temp0', temp0)
      self.mdin.change('cntrl', 'ntr', int(restrained))
      self.mdin.change('cntrl', 'ig', -1)
      self.mdin.change('cntrl', 'ntc', 2)
      self.mdin.change('cntrl', 'ntf', 2)
      self.mdin.change('cntrl', 'ntpr', ntpr)
      self.mdin.change('cntrl', 'ntwx', ntwx)
      self.mdin.change('cntrl', 'ntwr', ntwr)
      self.mdin.change('cntrl', 'dt', 0.002)
      if thermostat.lower() == 'langevin': 
         self.mdin.change('cntrl', 'ntt', 3) 
         self.mdin.change('cntrl', 'gamma_ln', thermostat_param)
      elif thermostat.lower() == 'berendsen':
         self.mdin.change('cntrl', 'ntt', 1)
         self.mdin.change('cntrl', 'tautp', thermostat_param)
      else:
         raise ValueError('thermostat must be langevin or berendsen!')
      if restrained:
         self.mdin.change('cntrl', 'restraint_wt', rst_wt)
         self.mdin.change('cntrl', 'restraintmask', rst_mask)

      if slow:
         self.add_lines("&wt\n   TYPE='TEMP0', ISTEP1=0, ISTEP2=%d" %
                        ((nstlim * 2) // 3))
         self.add_lines("   VALUE1=%.2f, VALUE2=%.2f,\n/\n&wt TYPE='END' /" % (
                        tempi, temp0))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Production(BaseType):
   """ This is for running straight-out MD """

   program = 'pmemd'

   # ==============================

   def __init__(self, amber_sys, nstlim=1000, igb=5, restrained=False,
                rst_wt=2.0, rst_mask='@CA,C,O,N', temp0=300.0,
                thermostat='langevin', thermostat_param=5.0, barostat=None,
                barostat_param=5.0, ntpr=1000, ntwr=10000, ntwx=1000,
                restart=True):
      # Create the mdin instance
      self.mdin = mdin(self.program) 
      # Do we have velocities?
      if restart:
         self.mdin.change('cntrl', 'ntx', 5)
         self.mdin.change('cntrl', 'irest', 1)
      # Set up periodic/non-periodic settings
      if amber_sys.periodic():
         if barostat.lower() == 'berendsen':
            self.mdin.change('cntrl', 'ntb', 2)
            self.mdin.change('cntrl', 'ntp', 1)
            self.mdin.change('cntrl', 'taup', barostat_param)
         elif barostat:
            raise ValueError('barostat must be "berendsen" or None/False')
         else:
            self.mdin.change('cntrl', 'ntb', 1)
      else:
         # Type-check igb
         if type(igb).__name__ != 'int':
            raise TypeError('igb must be an integer!')
         self.mdin.change('cntrl', 'ntb', 0)
         if not igb in [1, 2, 5, 7, 8]:
            raise ValueError('For non-periodic systems, igb must be ' +
                             '1, 2, 5, 7, or 8')
         self.mdin.change('cntrl', 'igb', igb)
         self.mdin.change('cntrl', 'cut', 1000.0)
      self.mdin.change('cntrl', 'ioutfm', 1) # always NetCDF!
      self.mdin.change('cntrl', 'ntr', int(restrained))
      self.mdin.change('cntrl', 'ig', -1)
      self.mdin.change('cntrl', 'ntc', 2)
      self.mdin.change('cntrl', 'ntf', 2)
      self.mdin.change('cntrl', 'ntpr', ntpr)
      self.mdin.change('cntrl', 'ntwx', ntwx)
      self.mdin.change('cntrl', 'ntwr', ntwr)
      self.mdin.change('cntrl', 'dt', 0.002)
      if thermostat.lower() == 'langevin':
         self.mdin.change('cntrl', 'ntt', 3)
         self.mdin.change('cntrl', 'gamma_ln', thermostat_param)
      elif thermostat.lower() == 'berendsen':
         self.mdin.change('cntrl', 'ntt', 1)
         self.mdin.change('cntrl', 'tautp', thermostat_param)
      else:
         raise ValueError('thermostat must be langevin or berendsen!')
      if restrained:
         self.mdin.change('cntrl', 'restraint_wt', rst_wt)
         self.mdin.change('cntrl', 'restraintmask', rst_mask)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ConstantpH(Production):
   """ A constant pH simulation. It's a type of production """

   program = 'sander'

   # ==============================

   def __init__(self, amber_sys, nstlim=1000, igb=5, restrained=False,
                rst_wt=2.0, rst_mask='@CA,C,O,N', temp0=300.0,
                thermostat='langevin', thermostat_param=5.0, barostat=None,
                barostat_param=5.0, ntpr=1000, ntwr=10000, ntwx=1000,
                restart=True, ntcnstph=5, ntrelax=500, solvph=7.0):


      Production.__init__(self, amber_sys, nstlim, igb, restrained, rst_wt, 
                          rst_mask, temp0, thermostat, thermostat_param,
                          barostat, barostat_param, ntpr, ntwr, ntwx, restart)
      
      if amber_sys.periodic():
         self.mdin.change('cntrl', 'icnstph', 2)
         self.mdin.change('cntrl', 'ntrelax', ntrelax)
      else:
         self.mdin.change('cntrl', 'icnstph', 2)

      self.mdin.change('cntrl', 'ntcnstph', ntcnstph)
      self.mdin.change('cntrl', 'solvph', solvph)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':

   parser = OptionParser(epilog='This program sets up a ~/.mpiexec_cmd file' +
          ' that is used to define how an MPI program is run on this computer')
   
   parser.add_option('--mpiexec-string', dest='mpi_string', default=None,
                     help='mpiexec string to run MPI programs. Written to ' +
                     '~/.mpiexec_cmd')
   
   (opt, args) = parser.parse_args()

   if len(args) != 0 or not opt.mpi_string:
      print 'Bad command-line arguments!'
      parser.print_help()
      sys.exit(1)

   mpi_file = open(os.path.join(os.getenv('HOME'), '.mpiexec_cmd'), 'w')
   mpi_file.write(opt.mpi_string)
   mpi_file.close()
