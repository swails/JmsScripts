# This module will create a sander/pmemd input

class mdin:

   import cntrl, ewald
   from sys import stderr, stdout

   # Variable list: These are dictionaries corresponding to the different
   # namelists that can be specified. They will be defined in __init__

   # cntrl
   # cntrl_defaults
   # ewald
   # ewald_defaults

   program = ''

   def __init__(self, program = 'sander'):
      if program == "sander":
         self.cntrl = cntrl.cntrl.sander
         self.cntrl_defaults = cntrl.cntrl.sander
         self.ewald = ewald.ewald.sander
         self.ewald_defaults = ewald.ewald.sander
      elif program == "pmemd":
         self.cntrl = cntrl.cntrl.pmemd
         self.cntrl_defaults = cntrl.cntrl.pmemd
         self.ewald = ewald.ewald.pmemd
         self.ewald_defaults = ewald.ewald.pmemd
      else:
         print >> stderr, 'Error: program (%s) unrecognized!' % program
