# This module will create a sander/pmemd input
from cntrl import cntrl
from ewald import ewald
from sys import stderr, stdout

# =====================================================================================

def addOn(line, string, file):
   if len(line.strip()) == 0:
      return line + string
   elif len(line) + len(string) > 40:
      file.write(line + '\n')
      return ' ' + string
   else:
      return line + string

# =====================================================================================

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
   valid_namelists = []    # array with valid namelists for each program

# =====================================================================================

   def __init__(self, program = 'sander', verbosity = 1):
      self.program = program
      if self.program == "sander":
         self.cntrl_nml = self.cntrl_obj.sander
         self.ewald_nml = self.ewald_obj.sander
         self.valid_namelists = ['cntrl','ewald','qmmm']
      elif self.program == "pmemd":
         self.cntrl_nml = self.cntrl_obj.pmemd
         self.ewald_nml = self.ewald_obj.pmemd
         self.valid_namelists = ['cntrl','ewald']
      else:
         print >> stderr, 'Error: program (%s) unrecognized!' % self.program
         return

      self.ewald_nml_defaults = self.ewald_nml.copy()
      self.cntrl_nml_defaults = self.cntrl_nml.copy()

# =====================================================================================

   def write(self, filename = 'mdin'):

      if not self.check(): # check the validity of the parameters
         print >> stderr, 'Mdin file not written!'
         return -1

      # open the file for writing and write the header and &cntrl namelist
      file = open(filename,'w') 
      file.write('Mdin prepared by mdin.py\n')
      file.write('&cntrl\n')
      # automatic indent of single space
      line = ' '
      # add any variable that is different from the default to the mdin file
      for var in self.cntrl_nml.keys():
         if self.cntrl_nml[var] != self.cntrl_nml_defaults[var]:
            line = addOn(line, "%s=%s, " % (var, self.cntrl_nml[var]), file)

      # flush any remaining items that haven't yet been printed to the mdin file
      if len(line.strip()) != 0:
         file.write(line + '\n')

      # end the namelist
      file.write('/\n')

      # print the ewald namelist if any variables differ from the default
      line = ' '
      has_been_printed = False  # keep track if this namelist has been printed
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

# =====================================================================================
   
   def read(self, filename = 'mdin'):
      try:
         file = open(filename, 'r')
      except IOError:
         print >> stderr, 'File (%s) can\'t be opened for reading...' % filename
         return -1

      lines = file.readlines()
      file.close()

      # split up input file into separate fields according to their comma-delimitations
      blocks = []  # namelists in the order they appear
      block_fields = [] # array of arrays that correspond to entries in namelists found in "blocks" above
      inblock = False
      for i in range(len(lines)):
         if not inblock and not lines[i].strip().startswith('&'):
            continue
         elif not inblock and lines[i].startswith('&'):
            inblock = True
            block = lines[i].strip()[1:]
            blocks.append(block)    # add the name of the namelist to "blocks"
            block_fields.append([]) # add empty array to be filled with entries for given namelist
            if not block in self.valid_namelists:
               print >> stderr, 'Invalid namelist (%s) in input file (%s) for %s' % (lines[i].strip(), filename, self.program)
               return -1
         elif inblock and (lines[i].strip() == '/' or lines[i].strip() == '&end'):
            inblock = False
         elif inblock and lines[i].strip().startswith('&'):
            print >> stderr, 'Invalid input file (%s). Terminate each namelist before another is started' % filename
            return -1
         elif inblock:
            items = lines[i].strip().split(',')
            for j in range(len(items)):
               if len(items[j]) == 0:
                  items.pop(j)
            block_fields[len(block_fields)-1].extend(items)

      # combine any multi-element fields: e.g. rstwt=1,2,3,
      begin_field = -1
      for i in range(len(block_fields)):
         for j in range(len(block_fields[i])):
            if not '=' in block_fields[i][j]:
               if begin_field == -1:
                  print >> stderr, 'Invalid input file (%s).' % filename
                  return -1
               else:
                  block_fields[i][begin_field] += ',' + block_fields[i][j]
            else:
               begin_field = j

      # now parse through the options and add them to the dictionaries
      for i in range(len(block_fields)):
         for j in range(len(block_fields[i])):
            if not '=' in block_fields[i][j]:
               continue
            else:
               var = block_fields[i][j].split('=')
               self.change(blocks[i], var[0].strip(), var[1].strip())

# =====================================================================================

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

# =====================================================================================

   def check(self):
      return True

# =====================================================================================
