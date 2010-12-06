########################################################################
#                                                                      #
# This module contains an amber prmtop class that will read in all     #
# parameters and allow users to manipulate that data and write a new   #
# prmtop object. It will also extract parameters and write a frcmod.   #
# See readparm.README for more detailed description                    #
#                                                                      #
#          Last updated: 12/04/2010                                    #
#                                                                      #
########################################################################

########################## GPL LICENSE INFO ############################

#  Copyright (C) 2010  Jason Swails

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
   
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330,
#  Boston, MA 02111-1307, USA.

from sys import stderr, stdout
from datetime import datetime

try: # fsum is only part of python 2.6 or later, I think, so add in a substitute here.
   from math import fsum
except ImportError:
   def fsum(array):
      return sum(array)

# Global constants
AMBER_ELECTROSTATIC = 18.2223
AMBER_POINTERS = """
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
NUMEXTRA: number of extra points
NCOPY  : Number of copies for advanded simulations
"""
# These global variables provide a more natural way of accessing
# the various pointers.  Most useful if they're loaded into the
# top-level namespace.
NATOM  = 0
NTYPES = 1
NBONH  = 2
MBONA  = 3
NTHETH = 4
MTHETA = 5
NPHIH  = 6
MPHIA  = 7
NHPARM = 8
NPARM  = 9
NEXT   = 10
NRES   = 11
NBONA  = 12
NTHETA = 13
NPHIA  = 14
NUMBND = 15
NUMANG = 16
NPTRA  = 17
NATYP  = 18
NPHB   = 19
IFPERT = 20
NBPER  = 21
NGPER  = 22
NDPER  = 23
MBPER  = 24
MGPER  = 25
MDPER  = 26
IFBOX  = 27
NMXRS  = 28
IFCAP  = 29
NUMEXTRA= 30
NCOPY  = 31

# ++++  Functions associated with readparm objects...  ++++++++++++++++++++++++++++

def _parseFormat(format_string):  # parse a format statement and send back details
   """ Parses the fortran format statement. Recognizes ints, exponents, and strings.
       Returns the number of items/line, size of each item, and type of data """

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
      print >> stderr, 'Error: Unrecognized format "%s"!' % format_string
      return 1, 80, 'str'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class amberParm:
   """ Amber Topology (parm7 format) class. Gives low, and some high, level access to
       topology data. """

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, prm_name='prmtop'): # set up necessary variables
      """ Instantiates an amberParm object from data in prm_name and establishes validity
          based on presence of POINTERS and CHARGE sections """

      # instance variables:
      self.prm_name = prm_name # name of the prmtop file
      self.formats = {}        # dictionary of Fortran formats corresponding to each %FLAG
      self.parm_data = {}      # dictionary of all prmtop data referenced by %FLAG *NAME*
      self.flag_list = []      # ordered array of all %FLAGs in prmtop
      self.version = ''        # version string
      self.overwrite = False   # whether writeParm will overwrite filename prm_name
      self.exists = False      # Logical set to true if the prmtop exists
      self.valid = False       # Logical set to true if the prmtop is valid
      self.pointers = {}       # list of all the pointers in the prmtop
      LJ_types = {}            # dictionary in which each atom name pairs with its LJ atom type number
      LJ_radius = []           # ordered array of L-J radii in Angstroms -- indices are elements in LJ_types-1
      LJ_depth = []            # ordered array of L-J depths in kcal/mol analagous to LJ_radius

      self.rdparm() # read the prmtop
      self.valid = self.exists # if it exists, fill the pointers
      if self.exists:
         try: # try to load all of the pointers into the 
            self.pointers["NATOM"] = self.parm_data["POINTERS"][NATOM]
            self.pointers["NTYPES"] = self.parm_data["POINTERS"][NTYPES]
            self.pointers["NBONH"] = self.parm_data["POINTERS"][NBONH]
            self.pointers["MBONA"] = self.parm_data["POINTERS"][MBONA]
            self.pointers["NTHETH"] = self.parm_data["POINTERS"][NTHETH]
            self.pointers["MTHETA"] = self.parm_data["POINTERS"][MTHETA]
            self.pointers["NPHIH"] = self.parm_data["POINTERS"][NPHIH]
            self.pointers["MPHIA"] = self.parm_data["POINTERS"][MPHIA]
            self.pointers["NHPARM"] = self.parm_data["POINTERS"][NHPARM]
            self.pointers["NPARM"] = self.parm_data["POINTERS"][NPARM]
            self.pointers["NEXT"] = self.parm_data["POINTERS"][NEXT]
            self.pointers["NRES"] = self.parm_data["POINTERS"][NRES]
            self.pointers["NBONA"] = self.parm_data["POINTERS"][NBONA]
            self.pointers["NTHETA"] = self.parm_data["POINTERS"][NTHETA]
            self.pointers["NPHIA"] = self.parm_data["POINTERS"][NPHIA]
            self.pointers["NUMBND"] = self.parm_data["POINTERS"][NUMBND]
            self.pointers["NUMANG"] = self.parm_data["POINTERS"][NUMANG]
            self.pointers["NPTRA"] = self.parm_data["POINTERS"][NPTRA]
            self.pointers["NATYP"] = self.parm_data["POINTERS"][NATYP]
            self.pointers["NPHB"] = self.parm_data["POINTERS"][NPHB]
            self.pointers["IFPERT"] = self.parm_data["POINTERS"][IFPERT]
            self.pointers["NBPER"] = self.parm_data["POINTERS"][NBPER]
            self.pointers["NGPER"] = self.parm_data["POINTERS"][NGPER]
            self.pointers["NDPER"] = self.parm_data["POINTERS"][NDPER]
            self.pointers["MBPER"] = self.parm_data["POINTERS"][MBPER]
            self.pointers["MGPER"] = self.parm_data["POINTERS"][MGPER]
            self.pointers["MDPER"] = self.parm_data["POINTERS"][MDPER]
            self.pointers["IFBOX"] = self.parm_data["POINTERS"][IFBOX]
            self.pointers["NMXRS"] = self.parm_data["POINTERS"][NMXRS]
            self.pointers["IFCAP"] = self.parm_data["POINTERS"][IFCAP]
            self.pointers["NUMEXTRA"] = self.parm_data["POINTERS"][NUMEXTRA]
            self.pointers["NCOPY"] = self.parm_data["POINTERS"][NCOPY]
            self.valid = True
         except KeyError:
            print >> stderr, 'Error: POINTERS flag not found! Likely a bad AMBER topology file.'
            self.valid = False
         except IndexError:
            if (len(self.parm_data['POINTERS'])) < 30:
               print >> stderr, 'Error: Fewer integers in POINTERS section than expected! Likely a bad AMBER topology file.'
               self.valid = False

      if self.valid and 'LENNARD_JONES_ACOEF' in self.parm_data.keys() and 'LENNARD_JONES_BCOEF' in self.parm_data.keys():
         try:
            self.fill_LJ() # fill LJ arrays with LJ data for easy manipulations
         except:
            print >> stderr, 'Warning: Problem parsing L-J 6-12 parameters.'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __str__(self):
      """ Returns the name of the topology file as its string representation """
      return self.prm_name

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ptr(self,pointer):
      """ Returns the value of the given pointer, and converts to upper-case so it's case-insensitive.
          A pointer that doesn't exist is met with an error message and a list of valid pointers """

      global AMBER_POINTERS
      return self.pointers[pointer.upper()]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rdparm(self):   # read topology file and load all data into arrays/dictionaries
      """ Reads the topology file and loads data in parm_data array """
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
         self.exists = False
         self.valid = False
         return

      self.exists = True
      prmlines = prmtop.readlines() # load all lines into memory
      prmtop.close() # close the file now

      for i in range(len(prmlines)):

         if prmlines[i][0:8] == '%VERSION':
            if self.version != '':
               print >> stderr, 'Warning: %VERSION string defined multiple times in %s.' % self.prm_name
            self.version = prmlines[i].strip()

         elif prmlines[i][0:5] == '%FLAG':
            current_flag = prmlines[i][6:].strip()
            self.formats[current_flag] = ''
            self.parm_data[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif prmlines[i][0:7] == '%FORMAT':
            self.formats[current_flag] = prmlines[i][8:len(prmlines[i].strip())-1]
            number_items_perline, size_item, dat_type = _parseFormat(self.formats[current_flag])
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
      """ Writes the current data in parm_data into a new topology file with a given name. Will not
          overwrite the original prm_name unless the overwrite variable is set to True. """
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
      datestring = "DATE = %02d/%02d/%02d  %02d:%02d:%02d" % ( now.month, now.day, now.year % 100, 
                                                               now.hour, now.minute, now.second)
      for i in range(len(self.version)): # replace the date in version string
         if self.version[i:i+2] == "DA":
            self.version = self.version[:i] + datestring
            break

      # convert charges back to amber charges...
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] *= AMBER_ELECTROSTATIC

      # write version to top of prmtop file
      new_prm.write('%s\n' % self.version)

      # write data back to prmtop file, inserting blank line if it's an empty field
      for i in range(len(self.flag_list)):
         flag = self.flag_list[i]
         new_prm.write('%%FLAG %s\n' % flag)
         new_prm.write('%%FORMAT(%s)\n' % self.formats[flag])
         number_items_perline, size_item, dat_type = _parseFormat(self.formats[flag])
         if dat_type == 'dec':
            decnum = int(self.formats[flag].split('E')[1].split('.')[1])
         line = ''
         num_items = 0
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         for j in range(len(self.parm_data[flag])): # write data in new_prm
            if dat_type == 'dec':
               line += ('%%%s.%sE' % (size_item, decnum)) % self.parm_data[flag][j] 
            elif dat_type == 'int':
               line += ('%%%sd' % size_item) % self.parm_data[flag][j] 
            else:
               line += ('%s' % self.parm_data[flag][j]).ljust(size_item)

            num_items += 1
            if num_items == number_items_perline: # flush line to prmtop
               new_prm.write(line + '\n')
               line = ''
               num_items = 0
         if len(line.strip()) > 0:
            new_prm.write(line + '\n') # flush what's left to the prmtop

      new_prm.close() # close new prmtop

      # eliminate multiplicative constant on charges to reduce to fraction-e charges
      for i in range(len(self.parm_data["CHARGE"])):
         self.parm_data["CHARGE"][i] /= AMBER_ELECTROSTATIC

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def totMass(self):
      """Returns total mass of the system"""
      return fsum(self.parm_data["MASS"])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def totCharge(self):
      """Returns total charge of the system"""
      return fsum(self.parm_data["CHARGE"])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def frcmod(self, frcmod="frcmod"):
      """Prints an Frcmod file that contains every parameter found in prmtop"""
      from math import pi, pow

      print >> stderr, "Warning: amberParm.Frcmod() does not work for 10-12 non-" \
            + "bonded interactions or variable 1-4 scaling prmtops yet."

      def getMatches(entry, array):
         counter = 0
         for i in range(len(array)):
            if array[i][0:11] == entry:
               counter += 1
         return counter

      file = open(frcmod, 'w')

      found_atomtypes = [] # store all of the atom types that have been used for masses
      atom_type_nums = [] # store the index of which atom type it is
      found_bondtypes = [] # store all of the bond types that have been found
      found_angletypes = [] # store all of the angle types that have been found
      stored_dihedtypes = [] # store all of the dihedral types that have been found
      stored_impropers = [] # storage for all improper dihedral parameters
      unique_diheds = [] # storage for all of the unique dihedral parameters

      # write the title
      file.write("Force field created from parameters in %s\n" % self.prm_name)

      # First we have to write the mass 
      file.write("MASS\n")
      for i in range(self.pointers["NATOM"]):
         # make sure we haven't found this atom type yet
         is_found = False
         for j in range(len(found_atomtypes)):
            if self.parm_data["AMBER_ATOM_TYPE"][i] == found_atomtypes[j]:
               is_found = True

         if is_found:
            continue

         # not found: now print out information and consider it found
         found_atomtypes.append(self.parm_data["AMBER_ATOM_TYPE"][i])
         atom_type_nums.append(self.parm_data["ATOM_TYPE_INDEX"][i])
         file.write("%s%6.3f\n" % (self.parm_data["AMBER_ATOM_TYPE"][i].ljust(6), self.parm_data["MASS"][i]))

      file.write("\n")

      # Now we write the bonds
      file.write("BOND\n")
      # We need to collect terms from 2 different blocks -- BONDS_INC_HYDROGEN and BONDS_WITHOUT_HYDROGEN
      # See http://ambermd.org/formats.html to get the details of how to parse this. The pointers for each
      # of these are NBONH and MBONA. Do not-including H first, then do H-included.
      for i in range(self.pointers["NBONA"]):
         start_index = i * 3
         # This is the bond... see if it's been found before
         bond = "%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["BONDS_WITHOUT_HYDROGEN"][start_index]/3].ljust(2), 
                           self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["BONDS_WITHOUT_HYDROGEN"][start_index+1]/3].ljust(2))
         is_found = False
         for j in range(len(found_bondtypes)):
            if bond == found_bondtypes[j]:
               is_found = True
               break

         if is_found:
            continue

         # not found: now print out information and consider it found
         found_bondtypes.append(bond)
         file.write("%s   %8.3f  %6.3f\n" % (bond, 
                  self.parm_data["BOND_FORCE_CONSTANT"][self.parm_data["BONDS_WITHOUT_HYDROGEN"][start_index+2]-1],
                  self.parm_data["BOND_EQUIL_VALUE"][self.parm_data["BONDS_WITHOUT_HYDROGEN"][start_index+2]-1]     ))

      found_bondtypes = []  # free up this memory now that we're done with bonds without hydrogen

      for i in range(self.pointers["NBONH"]):
         start_index = i * 3
         # This is the bond... see if it's been found before
         bond = "%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["BONDS_INC_HYDROGEN"][start_index]/3].ljust(2),
                           self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["BONDS_INC_HYDROGEN"][start_index+1]/3].ljust(2))
         is_found = False
         for j in range(len(found_bondtypes)):
            if bond == found_bondtypes[j]:
               is_found = True
               break

         if is_found:
            continue

         # not found: now print out information and consider it found
         found_bondtypes.append(bond)
         file.write("%s   %8.3f  %6.3f\n" % (bond, 
                     self.parm_data["BOND_FORCE_CONSTANT"][self.parm_data["BONDS_INC_HYDROGEN"][start_index+2]-1],
                     self.parm_data["BOND_EQUIL_VALUE"][self.parm_data["BONDS_INC_HYDROGEN"][start_index+2]-1]     ))

      del found_bondtypes  # free up this memory now that we're done with all bonds

      file.write('\n')

      # Now we write the angles: same kind of deal as the bonds, but now we have 3 atoms instead of 2 to find
      file.write('ANGLE\n')
      for i in range(self.pointers["NTHETA"]):
         start_index = i * 4
         # This is the angle... see if it's been found before
         angle = "%s-%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_WITHOUT_HYDROGEN"][start_index]/3].ljust(2),
                             self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_WITHOUT_HYDROGEN"][start_index+1]/3].ljust(2),
                             self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_WITHOUT_HYDROGEN"][start_index+2]/3].ljust(2) )
         is_found = False
         for j in range(len(found_angletypes)):
            if angle == found_angletypes[j]:
               is_found = True
               break

         if is_found:
            continue

         # not found: now print out information and consider it found
         found_angletypes.append(angle)
         file.write("%s   %8.3f  %6.3f\n" % (angle,
               self.parm_data["ANGLE_FORCE_CONSTANT"][self.parm_data["ANGLES_WITHOUT_HYDROGEN"][start_index+3]-1],
               self.parm_data["ANGLE_EQUIL_VALUE"][self.parm_data["ANGLES_WITHOUT_HYDROGEN"][start_index+3]-1] * 180 / pi ))

      for i in range(self.pointers["NTHETH"]):
         start_index = i * 4
         # This is the angle... see if it's been found before
         angle = "%s-%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_INC_HYDROGEN"][start_index]/3].ljust(2),
                             self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_INC_HYDROGEN"][start_index+1]/3].ljust(2),
                             self.parm_data["AMBER_ATOM_TYPE"][self.parm_data["ANGLES_INC_HYDROGEN"][start_index+2]/3].ljust(2) )
         is_found = False
         for j in range(len(found_angletypes)):
            if angle == found_angletypes[j]:
               is_found = True
               break

         if is_found:
            continue

         # not found: now print out information and consider it found
         found_angletypes.append(angle)
         file.write("%s   %8.3f  %6.3f\n" % (angle,
               self.parm_data["ANGLE_FORCE_CONSTANT"][self.parm_data["ANGLES_INC_HYDROGEN"][start_index+3]-1],
               self.parm_data["ANGLE_EQUIL_VALUE"][self.parm_data["ANGLES_INC_HYDROGEN"][start_index+3]-1] * 180 / pi ))

      del found_angletypes # done with this, clear the memory

      file.write('\n')
      # now it's time to find the dihedrals

      for i in range(self.pointers["NPHIA"]):
         start_index = i * 5
         # atom1 - atom4 are actually atom# - 1. I only need to check for negative values in atom3 and atom4
         # negative in atom3 means it's multiterm, negative atom4 means it's an improper, so store it
         atom1 = self.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"][start_index]/3
         atom2 = self.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"][start_index+1]/3
         atom3 = self.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"][start_index+2]/3
         atom4 = self.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"][start_index+3]/3
         term  = self.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"][start_index+4]
         dihedral = "%s-%s-%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][atom1].ljust(2),self.parm_data["AMBER_ATOM_TYPE"][atom2].ljust(2), 
                      self.parm_data["AMBER_ATOM_TYPE"][abs(atom3)].ljust(2),self.parm_data["AMBER_ATOM_TYPE"][abs(atom4)].ljust(3))

         if atom4 < 0:
            dihedral = "%s %8.3f %8.3f %5.1f" % (dihedral, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         elif atom3 < 0: # if there's another term in the series
            dihedral = "%s %4i %8.3f %8.3f %5.1f" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, -self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         else:
            dihedral = "%s %4i %8.3f %8.3f %5.1f" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         if atom4 < 0: # if it's a *new* improper, store it if necessary
            is_found = False
            for j in range(len(stored_impropers)):
               if stored_impropers[j] == dihedral:
                  is_found = True
                  break
            if is_found:
               continue
            else:
               stored_impropers.append(dihedral)
         else:
            is_found = False
            for j in range(len(stored_dihedtypes)):
               if stored_dihedtypes[j] == dihedral:
                  is_found = True
                  break
            if is_found:
               continue
            else:
               stored_dihedtypes.append(dihedral)

      for i in range(self.pointers["NPHIH"]):
         start_index = i * 5
         # atom1 - atom4 are actually atom# - 1. I only need to check for negative values in atom3 and atom4
         # negative in atom3 means it's multiterm, negative atom4 means it's an improper, so store it
         atom1 = self.parm_data["DIHEDRALS_INC_HYDROGEN"][start_index]/3
         atom2 = self.parm_data["DIHEDRALS_INC_HYDROGEN"][start_index+1]/3
         atom3 = self.parm_data["DIHEDRALS_INC_HYDROGEN"][start_index+2]/3
         atom4 = self.parm_data["DIHEDRALS_INC_HYDROGEN"][start_index+3]/3
         term  = self.parm_data["DIHEDRALS_INC_HYDROGEN"][start_index+4]
         dihedral = "%s-%s-%s-%s" % (self.parm_data["AMBER_ATOM_TYPE"][atom1].ljust(2),self.parm_data["AMBER_ATOM_TYPE"][atom2].ljust(2), 
                      self.parm_data["AMBER_ATOM_TYPE"][abs(atom3)].ljust(2),self.parm_data["AMBER_ATOM_TYPE"][abs(atom4)].ljust(3))

         if atom4 < 0:
            dihedral = "%s %8.3f %8.3f %5.1f" % (dihedral, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         elif atom3 < 0: # if there's another term in the series
            dihedral = "%s %4i %8.3f %8.3f %5.1f" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, -self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         else:
            dihedral = "%s %4i %8.3f %8.3f %5.1f" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])

         if atom4 < 0: # if it's a *new* improper, store it if necessary
            is_found = False
            for j in range(len(stored_impropers)):
               if stored_impropers[j] == dihedral:
                  is_found = True
                  break
            if is_found:
               continue
            else:
               stored_impropers.append(dihedral)
         else:
            is_found = False
            for j in range(len(stored_dihedtypes)):
               if stored_dihedtypes[j] == dihedral:
                  is_found = True
                  break
            if is_found:
               continue
            else:
               stored_dihedtypes.append(dihedral)

      # Find unique dihedrals -- this part is necessary because of multiterm dihedrals and the fact
      # that the ordering is not necessarily what one would expect
      for i in range(len(stored_dihedtypes)):
         is_found = False
         for j in range(len(unique_diheds)):
            if stored_dihedtypes[i][0:11] == unique_diheds[j]:
               is_found = True
               break
         if is_found:
            continue
         else:
            unique_diheds.append(stored_dihedtypes[i][0:11])
         
      file.write("DIHE\n")
      for i in range(len(unique_diheds)): # now that we have all the unique dihedrals, 
         num_left = getMatches(unique_diheds[i], stored_dihedtypes)
         while num_left > 0:
            if num_left > 1:
               for j in range(len(stored_dihedtypes)):
                  if float(stored_dihedtypes[j][len(stored_dihedtypes[j])-6:]) < 0 and  \
                                          stored_dihedtypes[j][0:11] == unique_diheds[i]:
                     file.write(stored_dihedtypes.pop(j) + '\n')
                     num_left -= 1
                     break
            else:
               for j in range(len(stored_dihedtypes)):
                  if stored_dihedtypes[j][0:11] == unique_diheds[i]:
                     file.write(stored_dihedtypes.pop(j) + '\n')
                     num_left -= 1
                     break

      unique_diheds = []
      del stored_dihedtypes
      # now write impropers

      for i in range(len(stored_impropers)):
         is_found = False
         for j in range(len(unique_diheds)):
            if stored_impropers[i][0:11] == unique_diheds[j]:
               is_found = True
               break
         if is_found:
            continue
         else:
            unique_diheds.append(stored_impropers[i][0:11])

      file.write("\nIMPROPER\n")

      for i in range(len(unique_diheds)): # now that we have all the unique dihedrals, 
         num_left = getMatches(unique_diheds[i], stored_impropers)
         while num_left > 0:
            if num_left > 1:
               for j in range(len(stored_impropers)):
                  if float(stored_impropers[j][len(stored_impropers[j])-6:]) < 0 and  \
                                        stored_impropers[j][0:11] == unique_diheds[i]:
                     file.write(stored_impropers.pop(j) + '\n')
                     num_left -= 1
                     break
            else:
               for j in range(len(stored_impropers)):
                  if stored_impropers[j][0:11] == unique_diheds[i]:
                     file.write(stored_impropers.pop(j) + '\n')
                     num_left -= 1
                     break

      del unique_diheds, stored_impropers
      file.write('\n') # done with dihedrals and improper dihedrals

      # now it's time for the non-bonded terms. 
      file.write("NONB\n")
      for i in range(len(found_atomtypes)):
         file.write("%s  %8.4f %8.4f \n" % (found_atomtypes[i].ljust(2), self.LJ_radius[self.LJ_types[found_atomtypes[i]]-1],
                     self.LJ_depth[self.LJ_types[found_atomtypes[i]]-1]))

      del found_atomtypes # done with these now.

      print >> stdout, "Amber force field modification (%s) finished!" % frcmod
      file.close()
      return 0

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def fill_LJ(self):
      """ Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data from LENNARD_JONES_ACOEF
          and LENNARD_JONES_BCOEF sections of the prmtop files, by undoing the canonical combining rules. """
      self.LJ_radius = []  # empty LJ_radii so it can be re-filled
      self.LJ_depth = []   # empty LJ_depths so it can be re-filled
      self.LJ_types = {}   # empty LJ_types so it can be re-filled
      one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

      for i in range(self.pointers["NATOM"]): # fill the LJ_types array
         self.LJ_types[self.parm_data["AMBER_ATOM_TYPE"][i]] = self.parm_data["ATOM_TYPE_INDEX"][i]
         
      for i in range(self.pointers["NTYPES"]):
         lj_index = (i + 1) * (i + 2) / 2 - 1 # n(n+1)/2 th position adjusted for indexing from 0
         if self.parm_data["LENNARD_JONES_BCOEF"][lj_index] < 1.0e-6:
            self.LJ_radius.append(0)
            self.LJ_depth.append(0)
         else:
            factor = 2 * self.parm_data["LENNARD_JONES_ACOEF"][lj_index] / self.parm_data["LENNARD_JONES_BCOEF"][lj_index]
            self.LJ_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_depth.append(self.parm_data["LENNARD_JONES_BCOEF"][lj_index] / 2 / factor)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_LJ(self):
      """ Takes the values of the LJ_radius and LJ_depth arrays and recalculates the LENNARD_JONES_A/BCOEF
          topology sections from the canonical combining rules. """
      from math import sqrt

      index = 0

      for i in range(self.pointers["NTYPES"]):
         for j in range(i+1):
            rij = self.LJ_radius[i] + self.LJ_radius[j]
            wdij = sqrt(self.LJ_depth[i] * self.LJ_depth[j])
            self.parm_data["LENNARD_JONES_ACOEF"][index] = wdij * pow(rij, 12)
            self.parm_data["LENNARD_JONES_BCOEF"][index] = 2 * wdij * pow(rij, 6)
            index += 1

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
