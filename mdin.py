# This module will create a sander/pmemd input

class mdin:

   from cntrl import cntrl
   from ewald import ewald
   from sys import stderr, stdout

   # Variable list: These are dictionaries corresponding to the different
   # namelists that can be specified. They will be defined in __init__

   cntrl_obj = cntrl
   cntrl_obj_defaults = cntrl
   ewald_obj = ewald
   ewald_obj_defaults = ewald
   
   cntrl_nml = {}
   ewald_nml = {}
   cntrl_nml_defaults = {}
   ewald_nml_defaults = {}

   def __init__(self, program = 'sander'):
      if program == "sander":
         self.cntrl_nml = cntrl_obj.sander
         self.cntrl_nml_defaults = cntrl_obj_defaults.sander
         self.ewald_nml = ewald.
         self.ewald_nml_defaults = ewald_obj_defaults.sander
      elif program == "pmemd":
         self.cntrl = cntrl_obj.pmemd
         self.cntrl_defaults = cntrl_obj_defaults.pmemd
         self.ewald = ewald_obj.pmemd
         self.ewald_defaults = ewald_obj_defaults.pmemd
      else:
         print >> stderr, 'Error: program (%s) unrecognized!' % program
   
   def write(self, filename = 'mdin'):
   
   
   def read(self, filename = 'mdin'):
