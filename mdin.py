# This module will create a sander/pmemd input
from cntrl import cntrl
from ewald import ewald
from sys import stderr, stdout

def addOn(line, string, file):
   if len(line.strip()) == 0:
      return line + string
   elif len(line) + len(string) > 40:
      file.write(line + '\n')
      return ' ' + string
   else:
      return line + string

class mdin:

   cntrl_obj = cntrl() # object with cntrl namelist vars in a dictionary
   ewald_obj = ewald() # object with ewald namelist vars in a dictionary

   verbosity = 0 # verbosity level: 0 -- print nothing
                 #                  1 -- print errors
                 #                  2 -- print errors and warnings
                 #                  3 -- print errors, warnings, and notes
   
   cntrl_nml = {}          # dictionary with cntrl namelist vars
   cntrl_nml_defaults = {} # dictionary with default cntrl namelist vars
   ewald_nml = {}          # dictionary with ewald namelist vars
   ewald_nml_defaults = {} # dictionary with default ewald namelist vars

   def __init__(self, program = 'sander', verbosity = 1):
      self.program = program
      if self.program == "sander":
         self.cntrl_nml = self.cntrl_obj.sander
         self.ewald_nml = self.ewald_obj.sander
      elif self.program == "pmemd":
         self.cntrl_nml = self.cntrl_obj.pmemd
         self.ewald_nml = self.ewald_obj.pmemd
      else:
         print >> stderr, 'Error: program (%s) unrecognized!' % self.program
         return

      self.ewald_nml_defaults = self.ewald_nml.copy()
      self.cntrl_nml_defaults = self.cntrl_nml.copy()

   def write(self, filename = 'mdin'):

      if not self.check():
         print >> stderr, 'Mdin file not written!'
         return

      file = open(filename,'w')
      file.write('Mdin prepared by mdin.py\n')
      file.write('&cntrl\n')
      line = ' '
      for var in self.cntrl_nml.keys():
         if self.cntrl_nml[var] != self.cntrl_nml_defaults[var]:
            line = addOn(line, "%s=%s," % (var, self.cntrl_nml[var]), file)

      if len(line.strip()) != 0:
         file.write(line + '\n')
      file.write('/\n')

      line = ' '
      has_been_printed = False
      for var in self.ewald_nml.keys():
         if self.ewald_nml[var] != self.ewald_nml_defaults[var]:
            if (not has_been_printed):
               file.write('&ewald\n')
               has_been_printed = True
            line = addOn(line, "%s=%s," % (var, self.ewald_nml[var]), file)

      if len(line.strip()) != 0:
         file.write(line + '\n')
      if has_been_printed:
         file.write('/\n')

   
   def read(self, filename = 'mdin'):
      pass

   def change(self, namelist, variable, value): # change the value of a variable without adding new key-pair
      if namelist == "cntrl":
         if variable in self.cntrl_nml.keys():
            self.cntrl_nml[variable] = value
         else:
            print >> stderr, 'Unknown variable (%s) in &cntrl!' % variable
      elif namelist == 'ewald': 
         if variable in self.ewald_nml.keys():
            self.ewald_nml[variable] = value
         else:
            print >> stderr, 'Unknown variable (%s) in &ewald!' % variable
      else:
         print >> stderr, 'Unknown namelist (%s)!' % namelist

   def check(self):
      return True
