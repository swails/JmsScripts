# This module will create a sander/pmemd input
from cntrl import cntrl
from ewald import ewald
from pb import pb
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
   
# =====================================================================================

   def __init__(self, program = 'sander', verbosity = 1):
      # define instance data
      self.program = program        # which program we're creating the input file for
      self.cntrl_obj = cntrl()      # object with cntrl namelist vars in a dictionary
      self.ewald_obj = ewald()      # object with ewald namelist vars in a dictionary
      self.pb_obj = pb()            # object with pb namelist vars in a dictionary
      self.verbosity = 0            # verbosity level: 0 -- print nothing
                                    #                  1 -- print errors
                                    #                  2 -- print errors and warnings
                                    #                  3 -- print errors, warnings, and notes
      self.cards = []               # array that has all of the input cards that come after namelists
      self.cntrl_nml = {}           # dictionary with cntrl namelist vars
      self.cntrl_nml_defaults = {}  # dictionary with default cntrl namelist vars
      self.ewald_nml = {}           # dictionary with ewald namelist vars
      self.ewald_nml_defaults = {}  # dictionary with default ewald namelist vars
      self.pb_nml = {}              # dictionary with pb namelist vars
      self.pb_nml_defaults = {}     # dictionary with default pb namelist vars
      self.valid_namelists = []     # array with valid namelists for each program
      self.title = 'mdin prepared by mdin.py'   # title for the mdin file


      if self.program == "sander":
         self.cntrl_nml = self.cntrl_obj.sander
         self.ewald_nml = self.ewald_obj.sander
         self.pb_nml = self.pb_obj.sander
         self.valid_namelists = ['cntrl','ewald','qmmm','pb']
      if self.program == "sander.APBS":
         self.cntrl_nml = self.cntrl_obj.sander
         self.pb_nml = self.pb_obj.sanderAPBS
         self.valid_namelists = ['cntrl','apbs']
      elif self.program == "pmemd":
         self.cntrl_nml = self.cntrl_obj.pmemd
         self.ewald_nml = self.ewald_obj.pmemd
         self.valid_namelists = ['cntrl','ewald']
      else:
         print >> stderr, 'Error: program (%s) unrecognized!' % self.program
         return

      self.cntrl_nml_defaults = self.cntrl_nml.copy()
      self.ewald_nml_defaults = self.ewald_nml.copy()
      self.pb_nml_defaults = self.pb_nml.copy()

# =====================================================================================

   def write(self, filename = 'mdin'):
      
      if not self.check(): # check the validity of the parameters
         print >> stderr, 'mdin file not written!'
         return -1

      # open the file for writing and write the header and &cntrl namelist
      file = open(filename,'w')
      file.write(self.title + '\n')
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

      # flush any remaining items that haven't been printed to the mdin file
      if len(line.strip()) != 0:
         file.write(line + '\n')

      # end the namelist
      if has_been_printed:
         file.write('/\n')

      # print the pb namelist if any variables differ from the original
      line = ' '
      has_been_printed = False # keep track if this namelist has been printed
      for var in self.pb_nml.keys():
         if self.pb_nml[var] != self.pb_nml_defaults[var]:
            if (not has_been_printed):
               if self.program == 'sander.APBS':
                  file.write('&apbs\n')
               else:
                  file.write('&pb\n')
               has_been_printed = True
            line = addOn(line,'%s=%s' % (var, self.pb_nml[var]), file)

      # flush any remaining items that haven't been printed to the mdin file
      if len(line.strip()) != 0:
         file.write(line + '\n')

      # end the namelist
      if has_been_printed:
         file.write('/\n')

      # Write the cards to the input file
      for i in range(len(self.cards)):
         file.write(self.cards[i].strip() + '\n')
      if len(self.cards) != 0:
         file.write('END\n')
      
      file.close()

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
         elif not inblock and lines[i].strip().startswith('&'):
            inblock = True
            block = lines[i].strip()[1:].lower()
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
      
      variable = variable.lower()
      
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
      elif namelist == 'pb':
         if variable in self.pb_nml.keys():
            self.pb_nml[variable] = value
         else:
            print >> stderr, 'Unknown variable (%s) in &pb!' % variable
      else:
         print >> stderr, 'Unknown namelist (%s)!' % namelist

# =====================================================================================

   def check(self):
      return True

# =====================================================================================

   def SHAKE(self):
      self.change('cntrl','ntf', 2)
      self.change('cntrl','ntc', 2)
      self.change('cntrl','dt', 0.002)

# =====================================================================================

   def constPressure(self, press=1.0, taup=1.0):
      self.change('cntrl','ntb', 2)
      self.change('cntrl','ntp', 1)
      self.change('cntrl','pres0', press)
      self.change('cntrl','taup', taup)

# =====================================================================================

   def constVolume(self):
      self.change('cntrl','ntb', 1)
      self.change('cntrl','ntp', 0)

# =====================================================================================

   def constTemp(self, ntt=3, temp=300.0, gamma_ln=2.0, ig=-1, tautp=1.0):
      self.change('cntrl','ntt', ntt)
      self.change('cntrl','temp0', temp)
      self.change('cntrl','tempi', temp)
      self.change('cntrl','gamma_ln', gamma_ln)
      self.change('cntrl','ig', ig)
      self.change('cntrl','tautp', tautp)

# =====================================================================================

   def constpH(self, solvph=7.0, igb=2, ntcnstph=10):
      self.change('cntrl','icnstph', 1)
      self.change('cntrl','solvph', solvph)
      self.change('cntrl','ntcnstph', ntcnstph)
      self.change('cntrl','igb', igb)
      self.change('cntrl','ntb', 0)
      self.change('cntrl','saltcon', 0.1)

# =====================================================================================

   def restrainHeavyAtoms(self, restraint_wt=0.0):
      self.change('cntrl','ntr', 1)
      self.change('cntrl','restraint_wt', restraint_wt)
      self.change('cntrl','restraintmask', '!@H=')

# =====================================================================================

   def restrainBackbone(self, restraint_wt=0.0):
      self.change('cntrl','ntr', 1)
      self.change('cntrl','restraint_wt', restraint_wt)
      self.change('cntrl','restraintmask', '@N,CA,C')

# =====================================================================================

   def genBorn(self, igb=5, rgbmax=25.0):
      self.change('cntrl','igb', igb)
      self.change('cntrl','ntb', 0)
      self.change('cntrl','ntp', 0)
      self.change('cntrl','rgbmax', rgbmax)

# =====================================================================================

   def time(self, time=1000.0, dt=-1): # time in ps
      if dt == -1:
         if self.cntrl_nml['ntc'] == 2 and self.cntrl_nml['ntf'] == 2:
            dt = 0.002
         else:
            dt = 0.001
      time = int(time / dt)

      self.change('cntrl','dt', dt)
      self.change('cntrl','nstlim', time)
      self.change('cntrl','imin', 0)

# =====================================================================================

   def heat(self, tempi=0.0, temp0=300.0, ntt=3, tautp=5.0, ig=-1, gamma_ln=5.0):
      self.constVolume()
      self.change('cntrl','tempi', tempi)
      self.change('cntrl','temp0', temp0)
      self.change('cntrl','ntt', ntt)
      self.change('cntrl','tautp', tautp)
      self.change('cntrl','ig', ig)
      self.change('cntrl','gamma_ln', gamma_ln)

# =====================================================================================

   def restart(self,ntx=5):
      self.change('cntrl','irest',1)
      self.change('cntrl','ntx',ntx)

# =====================================================================================

   def TI(self, clambda=0.0):
      self.change('cntrl','clambda', clambda)
      self.change('cntrl','icfe',1)

# =====================================================================================

   def softcore_TI(self, scalpha=0.5, scmask='', crgmask='', logdvdl=0):
      self.change('cntrl','icfe',1)
      self.change('cntrl','ifsc',1)
      self.change('cntrl','scalpha',scalpha)
      self.change('cntrl','scmask',scmask)
      self.change('cntrl','crgmask',crgmask)
      self.change('cntrl','logdvdl',logdvdl)

# =====================================================================================

   def minimization(self, imin=1, maxcyc=1, ncyc=10, ntmin=1):
      self.change('cntrl','imin', imin)
      self.change('cntrl','maxcyc', maxcyc)
      self.change('cntrl','ncyc', ncyc)
      self.change('cntrl','ntmin', ntmin)

# =====================================================================================

   def AddCard(self, title='Residues in card', cardString='RES 1'):
      self.cards.append('%s\n%s\nEND' % (title, cardString))

# =====================================================================================
