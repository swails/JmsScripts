#!/usr/bin/env python
""" 
This contains classes for a bunch of different calculation types
"""

from subprocess import Popen, PIPE
from chemistry.amber.readparm import AmberParm
from utilities import which
from mdin import mdin
from optparse import OptionParser
from pbsjob import PBS_Script
import sys, os

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

# Main part of the program here.
if __name__ == '__main__':
   
   parser = OptionParser(usage="%prog [Options] <prmtop> <inpcrd>")
   parser.add_option('-m', '--min', dest='min_cyc', default=0, type='int',
                     help="Number of minimization steps to run. 0 means no " +
                     'minimization. Default 0')
   parser.add_option('-g', '--igb', dest='igb', default=5, type='int',
                     help='GB model to use for implicit solvent calculations.' +
                     ' Default is 5. Allowed 1, 2, 5, 7, or 8')
   parser.add_option('--restrain-min', dest='min_rst_wt', default=0,
                     type='float', help="Restraint weight for restrained " +
                     'minimization. 0 means no restraints. Default 0')
   parser.add_option('--min-mask', dest='min_mask', default='@CA,C,O,N',
                     help='Restraint mask for minimization. ' +
                     'Defaults to "@CA,C,O,N"')
   parser.add_option('--heat', dest='heat_nstlim', default=0, type='int',
                     help='Number of 2 fs time steps to heat your system. ' +
                     '0 means no heating. Default 0')
   parser.add_option('--heating-thermostat', dest='h_thermostat', 
                     default='langevin', help='Which thermostat to use for ' +
                     'heating step. Defaults to langevin. ' +
                     'Allowed berendsen or langevin (case-insensitive)')
   parser.add_option('--slow-heat', dest='slow_heat', default=False, 
                     action='store_true', help='Use NMR coupling to TEMP0 to ' +
                     'heat the system slowly. Default False.')
   parser.add_option('--restrain-heat', dest='heat_rst_wt', default=0,
                     type='float', help='Restraint weight for restrained ' +
                     'heating. 0 means no restraining. Default 0')
   parser.add_option('--heat-mask', dest='heat_mask', default='@CA,C,O,N',
                     help='Restraint mask for heating. Defaults to "@CA,C,O,N"')
   parser.add_option('--temp', dest='temp0', type='float', default=300.0,
                     help='Target temperature to run simulations at. This is ' +
                     'the temperature heated to during heating step ' +
                     '(if present), and that production MD is run at')
   parser.add_option('--equilibration', dest='equil_nstlim', default=0, 
                     type='int', help='Number of 2 fs time steps to run ' +
                     'constant pressure equilibration (for PBC)')
   parser.add_option('--restrain-equil', dest='equil_rstwt', default=0, 
                     type='float', help='Restraint weight for restrained ' +
                     'equilibration. 0 means no restraints. Default 0')
   parser.add_option('--equil-mask', dest='equil_mask', default='@CA,C,O,N',
                   help='Restraint mask for equilibration. Default "@CA,C,O,N"')
   parser.add_option('--prod-type', dest='prod_type', default='production',
                     help='Type of production simulation to run. Options are ' +
                     '"production" and "constant pH". Default "production"')
   parser.add_option('--nstlim', dest='nstlim', type='int', default=1000,
                     help='Number of 2 fs time steps to run production for. ' +
                     '0 means no production. Default 1000')
   parser.add_option('--ntpr', dest='ntpr', type='int', default=1000,
                     help='How frequently to print to mdout. Default 1000.')
   parser.add_option('--ntwr', dest='ntwr', type='int', default=10000,
                     help='How frequently to dump a restrt. Default 10000.')
   parser.add_option('--ntwx', dest='ntwx', type='int', default=1000,
                     help='How frequently to print to mdcrd. Default 1000.')
   parser.add_option('--ntcnstph', dest='ntcnstph', type='int', default=5,
                     help='How frequently to attempt protonation MC. Default 5')
   parser.add_option('--solvph', dest='solvph', type='float', default=7,
                     help='How frequently to attempt protonation MC. Default 5')
   parser.add_option('--thermostat',dest='prod_thermostat', default='berendsen',
                     help='Which thermostat to use for production MD. Allowed' +
                     ' berendsen, langevin, or none')
   parser.add_option('--t-couple', dest='temp_param', type='float', default=10,
                     help='Value of temperature coupling parameter. Default 10')
   parser.add_option('--barostat', dest='prod_barostat', default='none',
                     help='Barostat to use for production MD. Allowed ' +
                     'berendsen or none. Default none.')
   parser.add_option('--p-couple', dest='pres_param', type='float', default=10,
                     help='Value of pressure coupling parameter. Default 10')
   parser.add_option('--restrain-prod', dest='prod_rstwt', default=0, 
                     type='float', help='Restraint weight for restrained ' +
                     'production. 0 means no restraints. Default 0')
   parser.add_option('--prod-mask', dest='prod_mask', default='@CA,C,O,N',
                     help='Restraint mask for production. Default "@CA,C,O,N"')
   parser.add_option('--ask', dest='ask', action='store_true', default=True,
                     help='Print submission command and ask prior to ' +
                     'submitting. Default behavior. Overrides --force')
   parser.add_option('--force', dest='ask', action='store_false', default=True,
                     help='Do not ask for permission before submitting job')

   (opt, args) = parser.parse_args()

   if len(args) != 2:
      print >> sys.stderr, 'Bad command-line arguments.'
      parser.print_help()
      sys.exit(1)

   # First set up the minimization
