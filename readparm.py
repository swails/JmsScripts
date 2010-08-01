from sys import stderr, stdout
from datetime import datetime

# Global constants
AMBER_ELECTROSTATIC = 18.2223
POINTER_VARIABLES = """
NATOM  : total number of atoms 
NTYPES : total number of distinct atom types
NBONH  : number of bonds containing hydrogen
MBONA  : number of bonds not containing hydrogen
NTHETH : number of angles containing hydrogen
MTHETA : number of angles not containing hydrogen
NPHIH  : number of dihedrals containing hydrogen
MPHIA  : number of dihedrals not containing hydrogen
NHPARM : currently not used
NPARM  : currently not used
NEXT   : number of excluded atoms
NRES   : number of residues
NBONA  : MBONA + number of constraint bonds
NTHETA : MTHETA + number of constraint angles
NPHIA  : MPHIA + number of constraint dihedrals
NUMBND : number of unique bond types
NUMANG : number of unique angle types
NPTRA  : number of unique dihedral types
NATYP  : number of atom types in parameter file, see SOLTY below
NPHB   : number of distinct 10-12 hydrogen bond pair types
IFPERT : set to 1 if perturbation info is to be read in
NBPER  : number of bonds to be perturbed
NGPER  : number of angles to be perturbed
NDPER  : number of dihedrals to be perturbed
MBPER  : number of bonds with atoms completely in perturbed group
MGPER  : number of angles with atoms completely in perturbed group
MDPER  : number of dihedrals with atoms completely in perturbed groups
IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
NMXRS  : number of atoms in the largest residue
IFCAP  : set to 1 if the CAP option from edit was specified
"""

# ++++  Functions associated with readparm objects...  ++++++++++++++++++++++++++++

def parseFormat(format_string):  # parse a format statement and send back details

   if 'a' in format_string: # this is a string
      format_parts = format_string.split('a')
      return int(format_parts[0]), int(format_parts[1]), 'str'

   elif 'I' in format_string: # this is an integer
      format_parts = format_string.split('I')
      return int(format_parts[0]), int(format_parts[1]), 'int'

   elif 'E' in format_string: # this is a floating point decimal
      format_parts = format_string.split('E')
      decimal_parts = format_parts[1].split('.')
      return int(format_parts[0]), int(decimal_parts[0]), 'dec'

   else:
      print >> stderr, 'Error: Unrecognized format "{0}"!'.format(format_string)
      return 1, 80, 'str'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class amberParm:

   # variables used in amberParm
   formats = {}
   parm_data = {}
   flag_list = []
   version = ''
   prm_name = ''
   overwrite = False
   prm_exists = False
   pointers = {}

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, prm_name): # set up necessary variables

      self.prm_name = prm_name # name of the prmtop file
      self.rdparm()
      try: # try to load all of the pointers into the 
         self.pointers["NATOM"] = self.parm_data["POINTERS"][0]
         self.pointers["NTYPES"] = self.parm_data["POINTERS"][1]
         self.pointers["NBONH"] = self.parm_data["POINTERS"][2]
         self.pointers["MBONA"] = self.parm_data["POINTERS"][3]
         self.pointers["NTHETH"] = self.parm_data["POINTERS"][4]
         self.pointers["MTHETA"] = self.parm_data["POINTERS"][5]
         self.pointers["NPHIH"] = self.parm_data["POINTERS"][6]
         self.pointers["MPHIA"] = self.parm_data["POINTERS"][7]
         self.pointers["NHPARM"] = self.parm_data["POINTERS"][8]
         self.pointers["NPARM"] = self.parm_data["POINTERS"][9]
         self.pointers["NEXT"] = self.parm_data["POINTERS"][10]
         self.pointers["NRES"] = self.parm_data["POINTERS"][11]
         self.pointers["NBONA"] = self.parm_data["POINTERS"][12]
         self.pointers["NTHETA"] = self.parm_data["POINTERS"][13]
         self.pointers["NPHIA"] = self.parm_data["POINTERS"][14]
         self.pointers["NUMBND"] = self.parm_data["POINTERS"][15]
         self.pointers["NUMANG"] = self.parm_data["POINTERS"][16]
         self.pointers["NPTRA"] = self.parm_data["POINTERS"][17]
         self.pointers["NATYP"] = self.parm_data["POINTERS"][18]
         self.pointers["NPHB"] = self.parm_data["POINTERS"][19]
         self.pointers["IFPERT"] = self.parm_data["POINTERS"][20]
         self.pointers["NBPER"] = self.parm_data["POINTERS"][21]
         self.pointers["NGPER"] = self.parm_data["POINTERS"][22]
         self.pointers["NDPER"] = self.parm_data["POINTERS"][23]
         self.pointers["MBPER"] = self.parm_data["POINTERS"][24]
         self.pointers["MGPER"] = self.parm_data["POINTERS"][25]
         self.pointers["MDPER"] = self.parm_data["POINTERS"][26]
         self.pointers["IFBOX"] = self.parm_data["POINTERS"][27]
         self.pointers["NMXRS"] = self.parm_data["POINTERS"][28]
         self.pointers["IFCAP"] = self.parm_data["POINTERS"][29]
      except KeyError:
         print >> stderr, 'Error: POINTERS flag not found! Likely a bad AMBER topology file.'
      except IndexError:
         print >> stderr, 'Error: Fewer integers in POINTERS section than expected! Likely a bad AMBER topology file.'


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ptr(self,pointer):
      global POINTER_VARIABLES
      try:
         return self.pointers[pointer.upper()]
      except KeyError:
         print >> stderr, 'No pointer {0} in prmtop. Pointers are: '.format(pointer.upper()) + POINTER_VARIABLES
         return -1

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rdparm(self):   # read topology file and load all data into arrays/dictionaries
      # global variable(s)
      global AMBER_ELECTROSTATIC

      # variables necessary only to rdparm
      current_flag = ''
      dat_type = ''
      gathering_data = False
      number_items_perline = 0
      size_item = 0
      self.version = '' # reset all topology information each time rdparm is called
      self.formats = {}
      self.parm_data = {}
      self.flag_list = []

      try: # open up the prmtop file, catching if it doesn't exist
         prmtop = open(self.prm_name, 'r')
      except IOError:
         print >> stderr, 'Error: %s does not exist!' % self.prm_name
         return

      self.prm_exists = True
      prmlines = prmtop.readlines() # load all lines into memory
      prmtop.close() # close the file now

      for i in range(len(prmlines)):

         if prmlines[i][0:8] == '%VERSION':
            if self.version != '':
               print >> stderr, 'Warning: %VERSION string defined multiple times.'
            self.version = prmlines[i].strip()

         elif prmlines[i][0:5] == '%FLAG':
            current_flag = prmlines[i][6:].strip()
            self.formats[current_flag] = ''
            self.parm_data[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif prmlines[i][0:7] == '%FORMAT':
            self.formats[current_flag] = prmlines[i][8:len(prmlines[i].strip())-1]
            number_items_perline, size_item, dat_type = parseFormat(self.formats[current_flag])
            gathering_data = True

         elif gathering_data:
            position_inline = 0
            for j in range(number_items_perline):
               data_holder = prmlines[i][position_inline:position_inline+size_item].strip()
               if len(data_holder) == 0: # we've reached the end...
                  break # break out of the loop
               if dat_type == 'int': # if int, put held data item in as integer
                  self.parm_data[current_flag].append(int(data_holder))
               elif dat_type == 'dec': # if dec, put held data item in as float
                  self.parm_data[current_flag].append(float(data_holder))
               else: # otherwise, it stays a string
                  self.parm_data[current_flag].append(data_holder)
               position_inline += size_item # move to the next item

      # eliminate multiplicative constant on charges to reduce to fraction-e charges
      try:
         for i in range(len(self.parm_data["CHARGE"])):
            self.parm_data["CHARGE"][i] /= AMBER_ELECTROSTATIC
      except KeyError:
         print >> stderr, 'Error: CHARGE flag not found in prmtop! Likely a bad AMBER topology file.'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeParm(self, name):   # write a new prmtop with the current prmtop data
      # global variable(s)
      global AMBER_ELECTROSTATIC

      # make sure we want to write the new prmtop file
      if not self.overwrite and name == self.prm_name:
         print >> stderr, 'Error: Object\'s overwrite set to False! Will not overwrite original prmtop.'
         return
      elif self.version == '':
         print >> stderr, 'Error: Version string missing! Load prmtop data before writing new prmtop.'
         return

      # now that we know we will write the new prmtop file, open the new file
      new_prm = open(name, 'w')

      # get current time to put into new prmtop file
      now = datetime.now()
      datestring = "DATE = {0:02d}/{1:02d}/{2:02d}  {3:02d}:{4:02d}:{5:02d}".format(
                  now.month, now.day, now.year % 100, now.hour, now.minute, now.second)
      for i in range(len(self.version)): # replace the date in version string
         if self.version[i:i+2] == "DA":
            self.version = self.version[:i] + datestring
            break

      # convert charges back to amber charges...
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] *= AMBER_ELECTROSTATIC

      # write version to top of prmtop file
      new_prm.write('{0}\n'.format(self.version))

      # write data back to prmtop file, inserting blank line if it's an empty field
      for i in range(len(self.flag_list)):
         flag = self.flag_list[i]
         new_prm.write('%FLAG {0}\n'.format(flag))
         new_prm.write('%FORMAT({0})\n'.format(self.formats[flag]))
         number_items_perline, size_item, dat_type = parseFormat(self.formats[flag])
         if dat_type == 'dec':
            decnum = int(self.formats[flag].split('E')[1].split('.')[1])
         line = ''
         num_items = 0
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         for j in range(len(self.parm_data[flag])): # write data in new_prm
            if dat_type == 'dec':
               line += '{0:{1}.{2}e}'.format(self.parm_data[flag][j],size_item,decnum)
            elif dat_type == 'int':
               line += '{0:{1}d}'.format(self.parm_data[flag][j],size_item)
            else:
               line += '{0}'.format(self.parm_data[flag][j]).ljust(size_item)

            num_items += 1
            if num_items == number_items_perline: # flush line to prmtop
               new_prm.write(line + '\n')
               line = ''
               num_items = 0
         new_prm.write(line + '\n') # flush what's left to the prmtop

      new_prm.close() # close new prmtop

      # eliminate multiplicative constant on charges to reduce to fraction-e charges
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] /= AMBER_ELECTROSTATIC

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
