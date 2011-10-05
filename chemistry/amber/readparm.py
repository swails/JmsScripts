########################################################################
#                                                                      #
# This module contains an amber prmtop class that will read in all     #
# parameters and allow users to manipulate that data and write a new   #
# prmtop object. It will also extract parameters and write a frcmod.   #
# See readparm.README for more detailed description                    #
#                                                                      #
#          Last updated: 07/20/2011                                    #
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
from chemistry import exceptions
from chemistry import periodic_table
from math import ceil

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

   # Get rid of ( and ) specifiers in Fortran format strings. This is a hack, but
   # should work for existing chamber prmtop files

   format_string = format_string.replace('(','').replace(')','')

   # Fix case for E, I, and F

   format_string = format_string.replace('e','E')
   format_string = format_string.replace('i','I')
   format_string = format_string.replace('f','F')

   if 'a' in format_string: # this is a string
      format_parts = format_string.split('a')
      try:
         return int(format_parts[0]), int(format_parts[1]), 'str', None
      except:
         return 1, 80, 'str', None

   elif 'I' in format_string: # this is an integer
      format_parts = format_string.split('I')
      if len(format_parts[0].strip()) == 0: format_parts[0] = 1
      return int(format_parts[0]), int(format_parts[1]), 'int', None

   elif 'E' in format_string: # this is a floating point decimal
      format_parts = format_string.split('E')
      decimal_parts = format_parts[1].split('.')
      if len(format_parts[0].strip()) == 0: format_parts[0] = 1
      return int(format_parts[0]), int(decimal_parts[0]), 'dec', int(decimal_parts[1])
   
   elif 'F' in format_string: # this is also a floating point decimal
      format_parts = format_string.split('F')
      decimal_parts = format_parts[1].split('.')
      if len(format_parts[0].strip()) == 0: format_parts[0] = 1
      return int(format_parts[0]), int(decimal_parts[0]), 'dec', int(decimal_parts[1])

   else:
      print >> stderr, 'Warning: Unrecognized format "%s"! Assuming string.' % format_string
      return 1, 80, 'str', None

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# We now add all of the classes that are necessary to rebuild a topology file's set
# of pointers and stuff if you want to add or remove atoms, bonds, etc.

class _Atom(object):
   """ 
   An atom. Only fill these in _AtomList types, since _AtomList will keep track
   of when indexes and other stuff needs to be updated
   """
   #===================================================

   def __init__(self, atnum, parm=None, index_from=0):
      self.bond_partners = []
      self.angle_partners = []
      self.index_from = int(index_from)
      self.atnum = int(atnum)
      if parm and type(parm).__name__ == 'AmberParm': self.load_from_parm(parm)
   
   #===================================================

   def add_data(self, parm):
      """ 
      Writes this atom's data to the AmberParm object. Don't pitch a fit if we're
      missing some useless (unused) flags.
      """
      # Make sure we're indexing from 0
      atnum = self.atnum - self.index_from
      parm.parm_data['ATOM_NAME'][atnum] = self.atname[:4]
      parm.parm_data['CHARGE'][atnum] = self.charge
      parm.parm_data['MASS'][atnum] = self.mass
      parm.parm_data['ATOM_TYPE_INDEX'][atnum] = self.nb_idx
      parm.parm_data['NUMBER_EXCLUDED_ATOMS'][atnum] = (sum(self.bond_partners) +
                                                        sum(self.angle_partners))
      parm.parm_data['AMBER_ATOM_TYPE'][atnum] = self.atname[:4]
      try: parm.parm_data['JOIN_ARRAY'][atnum] = 0
      except KeyError: pass
      parm.parm_data['TREE_CHAIN_CLASSIFICATION'][atnum] = self.tree[:4]
      try: parm.parm_data['IROTAT'][atnum] = 0
      except KeyError: pass
      parm.parm_data['RADII'][atnum] = self.radii
      parm.parm_data['SCREEN'][atnum] = self.screen

   #===================================================

   def load_from_parm(self, parm):
      """ Load data from the AmberParm class """
      atnum = self.atnum - self.index_from # adjust for where we begin our indexing
      self.atname = parm.parm_data['ATOM_NAME'][atnum]
      self.charge = parm.parm_data['CHARGE'][atnum]
      self.mass = parm.parm_data['MASS'][atnum]
      self.nb_idx = parm.parm_data['ATOM_TYPE_INDEX'][atnum]
      self.atname = parm.parm_data['AMBER_ATOM_TYPE'][atnum]
      self.tree = parm.parm_data['TREE_CHAIN_CLASSIFICATION'][atnum]
      self.radii = parm.parm_data['RADII'][atnum]
      self.screen = parm.parm_data['SCREEN'][atnum]

   #===================================================
      
   def bond_to(self, other):
      """ 
      Log this atom as bonded to another atom. Check if this has already been
      added to the angle list. If so, remove it from there.
      """
      if self == other:
         raise exceptions.BondError("Cannot bond atom to itself!")
      if other in self.angle_partners:
         del self.angle_partners[self.angle_partners.index(other)]
      if other in self.bond_partners: return
      self.bond_partners.append(other)

   #===================================================
      
   def angle_to(self, other):
      """
      Log this atom as angled to another atom. Check if this has already been
      added to the bond list. If so, do nothing
      """
      if self == other:
         raise exceptions.BondError("Cannot angle an atom with itself!")
      if other in self.bond_partners or other in self.angle_partners: return
      self.angle_partners.append(other)
   
   #===================================================
      
   def atom_deleted(self, other):
      """ 
      This should be called for every atom in _AtomList so we know to delete
      it from any of the bond/angle partners
      """
      if other in self.bond_partners:
         del self.bond_partners[self.bond_partners.index(other)]
      if other in self.angle_partners:
         del self.angle_partners[self.angle_partners.index(other)]

   #===================================================

   def __eq__(self, other):
      return id(self) == id(other)
      
   def __ne__(self, other):
      return not _Atom.__eq__(self, other)

   def __gt__(self, other):
      return self.atnum > other.atnum

   def __lt__(self, other):
      return self.atnum < other.atnum

   def __ge__(self, other):
      return not _Atom.__lt__(self, other)

   def __le__(self, other):
      return not _Atom.__gt__(self, other)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _Bond(object):
   """ Bond class. Stores 2 atoms involved and force constant/equil value """

   #===================================================

   def __init__(self, atom1, atom2, bond_type):
      """ Bond constructor """
      # Order the atoms so the lowest atom # is first
      self.atom1 = atom1
      self.atom2 = atom2
      # Register each as bonded to the other
      self.atom1.bond_to(atom2)
      self.atom2.bond_to(atom1)
      # Load the force constant and equilibrium distance
      self.bond_type = bond_type

   #===================================================

   def write_info(self, parm, key, idx):
      """ Writes the bond info to the topology file. idx starts at 0 """
      parm.parm_data[key][3*idx  ] = 3*(self.atom1.atnum - self.atom1.index_from)
      parm.parm_data[key][3*idx+1] = 3*(self.atom2.atnum - self.atom2.index_from)
      parm.parm_data[key][3*idx+2] = self.bond_type.idx + 1
      # Now have this bond type write its info
      self.bond_type.write_info(parm)

   #===================================================
   
   def __contains__(self, thing):
      """ Quick and easy way to see if an _BondType or _Atom is in this _Bond """
      if type(thing).__name__ == '_BondType':
         return self.bond_type == thing
      elif type(thing).__name__ == '_Atom':
         return thing == self.atom1 or thing == self.atom2
      else: return False

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _BondType(object):
   """ A bond type """

   #===================================================

   def __init__(self, k, req, idx):
      """_BondType constructor. idx must start from 0!!! """
      self.idx = idx
      self.k = k
      self.req = req

   #===================================================

   def write_info(self, parm):
      """ Writes the bond parameters in the parameter file """
      parm.parm_data['BOND_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['BOND_EQUIL_VALUE'][self.idx] = self.req

   #===================================================

   def __eq__(self, other):
      return id(self) == id(other)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _Angle(object):
   """ Angle class. Stores 3 atoms involved and force constant/equil value """
      
   #===================================================

   def __init__(self, atom1, atom2, atom3, angle_type):
      """ Angle constructor """
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      # Register each as angled to the others
      self.atom1.angle_to(self.atom2)
      self.atom1.angle_to(self.atom3)
      self.atom2.angle_to(self.atom1)
      self.atom2.angle_to(self.atom3)
      self.atom3.angle_to(self.atom1)
      self.atom3.angle_to(self.atom2)
      # Load the force constant and equilibrium angle
      self.angle_type = angle_type

   #===================================================

   def write_info(self, parm, key, idx):
      """ Write the info to the topology file """
      parm.parm_data[key][4*idx  ] = 3*(self.atom1.atnum - self.atom1.index_from)
      parm.parm_data[key][4*idx+1] = 3*(self.atom2.atnum - self.atom2.index_from)
      parm.parm_data[key][4*idx+2] = 3*(self.atom3.atnum - self.atom3.index_from)
      parm.parm_data[key][4*idx+3] = self.angle_type.idx + 1
      # Now have this bond type write its info
      self.angle_type.write_info(parm)

   #===================================================
   
   def __contains__(self, thing):
      """ Quick and easy way to see if an _AngleType or _Atom is in this _Angle """
      if type(thing).__name__ == '_AngleType':
         return self.angle_type == thing
      elif type(thing).__name__ == '_Atom':
         return thing == self.atom1 or thing == self.atom2 or thing == self.atom3
      else: return False

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _AngleType(object):
   """ An angle type """
   #===================================================

   def __init__(self, k, theteq, idx):
      """ _AngleType constructor. idx must start from 0!!! """
      self.k = k
      self.theteq = theteq
      self.idx = idx

   #===================================================

   def write_info(self, parm):
      """ Writes the bond parameters in the parameter file """
      parm.parm_data['ANGLE_FORCE_CONSTANT'][self.idx] = self.k
      parm.parm_data['ANGLE_EQUIL_VALUE'][self.idx] = self.req

   #===================================================

   def __eq__(self, other):
      return id(self) == id(other)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _Dihedral(object):
   " Dihedral class with 4 atoms involved and force constant/periodicity/phase "
      
   #===================================================

   def __init__(self, atom1, atom2, atom3, atom4, dihed_type, signs):
      """ _Dihedral constructor. idx must start from 0!!! """
      self.atom1 = atom1
      self.atom2 = atom2
      self.atom3 = atom3
      self.atom4 = atom4
      # Load the force constant and equilibrium angle
      self.dihed_type = dihed_type
      self.signs = signs # is our 3rd or 4th term negative?

   #===================================================

   def write_info(self, parm, key, idx):
      """ Write the info to the topology file """
      parm.parm_data[key][5*idx  ] = 3*(self.atom1.atnum - self.atom1.index_from)
      parm.parm_data[key][5*idx+1] = 3*(self.atom2.atnum - self.atom2.index_from)
      parm.parm_data[key][5*idx+2] = 3*(self.atom3.atnum - self.atom3.index_from)\
                                   * self.signs[0]
      parm.parm_data[key][5*idx+3] = 3*(self.atom4.atnum - self.atom4.index_from)\
                                   * self.signs[1]
      parm.parm_data[key][5*idx+4] = self.dihed_type.idx + 1
      # Now have this bond type write its info
      self.dihed_type.write_info(parm)

   #===================================================

   def __contains__(self, thing):
      """ 
      Quick and easy way to find out if either a _DihedralType or _Atom is
      in this _Dihedral
      """
      if type(thing).__name__ == "_DihedralType":
         return self.dihedral_type == thing
      elif type(thing).__name__ == "_Atom":
         return thing == self.atom1 or thing == self.atom2 or \
                thing == self.atom3 or thing == self.atom4
      else: return False

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _DihedralType(object):
   """ A type of dihedral """

   #===================================================
   
   def __init__(self, phi_k, per, phase, idx):
      """ _DihedralType constructor """
      self.phi_k = phi_k
      self.per = per
      self.phase = phase

   #===================================================

   def write_info(self, parm):
      """ Write out the dihedral parameters """
      parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][self.idx] = self.phi_k
      parm.parm_data['DIHEDRAL_PERIODICITY'][self.idx] = self.per
      parm.parm_data['DIHEDRAL_PHASE'][self.idx] = self.phase
   
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _AtomList(list):
   """ Array of _Atoms """
   #===================================================

   def __init__(self, parm):
      self.parm = parm
      list.__init__(self, [_Atom(i,self.parm) for i in range(self.parm.ptr('natom'))])
      self.changed = False

   #===================================================

   def __delitem__(self, idx):
      """ Deletes this atom then re-indexes everybody else """
      for atm in self:
         atm.atom_deleted(self[idx])
      list.__delitem__(self, idx)
      self._reindex()
      self.changed = True
      # Now delete this item from all of our atomic property arrays in parm
      del self.parm.parm_data['ATOM_NAME'][idx]
      del self.parm.parm_data['CHARGE'][idx]
      del self.parm.parm_data['MASS'][idx]
      del self.parm.parm_data['ATOM_TYPE_INDEX'][idx]
      del self.parm.parm_data['NUMBER_EXCLUDED_ATOMS'][idx]
      del self.parm.parm_data['AMBER_ATOM_TYPE'][idx]
      del self.parm.parm_data['JOIN_ARRAY'][idx]
      del self.parm.parm_data['TREE_CHAIN_CLASSIFICATION'][idx]
      del self.parm.parm_data['IROTAT'][idx]
      del self.parm.parm_data['RADII'][idx]
      del self.parm.parm_data['SCREEN'][idx]
      # Adjust the NATOM pointer and broadcast that out
      self.parm.parm_data['POINTERS'][NATOM] -= 1
      self.parm.LoadPointers()

   #===================================================
   
   def _reindex(self):
      """ We have deleted an atom, so now we have to re-index everybody """
      for i in range(len(self)): self[i].atnum = i

   #===================================================

   def append(self, item):
      """ Don't allow this! """
      raise exceptions.AmberParmError("Cannot add to an _AtomList!")

   #===================================================

   def write_to_parm(self):
      """ Writes all of the atom data to the topology file """
      for atm in self: atm.add_data(parm)

   #===================================================

   def determine_exclusions(self):
      """ Figures out the EXCLUDED_ATOMS_LIST """
      self.parm.parm_data['EXCLUDED_ATOMS_LIST'] = []
      for atm in self:
         vals_to_add = list(atm.bond_partners) + list(atm.angle_partners)
         vals_to_add.sort()
         self.parm.parm_data['EXCLUDED_ATOMS_LIST'].extend(vals_to_add)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TypeList(list):
   """ Base class for all type lists """

   #===================================================

   def __init__(self, parm):
      """ Constructs a list of bond types from the topology file """
      self.parm = parm
      self._make_array()
      self.changed = False

   #===================================================

   def __delitem__(self, idx):
      """ Removes a bond from the bond type """
      list.__delitem__(self, idx)
      _TypeList._reindex(self)

   #===================================================

   def _reindex(self):
      """ Re-indexes after a deleted item """
      for i in range(len(self)): self[i].idx = i
   
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _BondTypeList(_TypeList):
   """ Bond type list """

   #===================================================

   def _make_array(self):
      list.__init__(self, [_BondType(self.parm.parm_data['BOND_FORCE_CONSTANT'][i],
                               self.parm.parm_data['BOND_EQUIL_VALUE'][i], i)
                               for i in range(self.parm.ptr('numbnd')) ])
      
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _AngleTypeList(_TypeList):
   """ Angle type list """

   #===================================================

   def _make_array(self):
      list.__init__(self, [_AngleType(self.parm.parm_data['ANGLE_FORCE_CONSTANT'][i],
                                self.parm.parm_data['ANGLE_EQUIL_VALUE'][i], i)
                                for i in range(self.parm.ptr('numang')) ])
      
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _DihedralTypeList(_TypeList):
   """ Dihedral type list """

   #===================================================

   def _make_array(self):
      list.__init__(self, [_DihedralType(self.parm.parm_data['DIHEDRAL_FORCE_CONSTANT'][i],
                                   self.parm.parm_data['DIHEDRAL_PERIODICITY'][i], 
                                   self.parm.parm_data['DIHEDRAL_PHASE'][i], i)
                                   for i in range(self.parm.ptr('nptra')) ])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TrackedList(list):
   """ This creates a list type that allows you to see if anything has changed """
   def __init__(self, arg=[]):
      self.changed = False
      list.__init__(self, arg)

   def __delitem__(self, idx):
      self.changed = True
      list.__delitem__(self, idx)

   def append(self, stuff):
      self.changed = True
      list.append(self, stuff)

   def extend(self, stuff):
      self.changed = True
      list.extend(self, stuff)

   def __delslice__(self, i, j):
      self.changed = True
      list.__delslice__(self, i, j)

   def __setitem__(self, i, y):
      self.changed = True
      list.__setitem__(self, i, y)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberParm(object):
   """ Amber Topology (parm7 format) class. Gives low, and some high, level access to
       topology data.
   """

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, prm_name='prmtop', rst7_name=''): # set up necessary variables
      """ Instantiates an AmberParm object from data in prm_name and establishes validity
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
      self.chamber = False     # Is this a chamber-generated prmtop file?
      self.pointers = {}       # list of all the pointers in the prmtop
      self.LJ_types = {}       # dictionary in which each atom name pairs with its LJ atom type number
      self.LJ_radius = []      # ordered array of L-J radii in Angstroms -- indices are elements in LJ_types-1
      self.LJ_depth = []       # ordered array of L-J depths in kcal/mol analagous to LJ_radius

      self.rdparm() # read the prmtop
      self.valid = self.exists # if it exists, fill the pointers
      if self.exists:
         try: # try to load all of the pointers into the 
            self.LoadPointers()
            self.valid = True
         except KeyError:
            print >> stderr, 'Error: POINTERS flag not found! Likely a bad AMBER topology file.'
            self.valid = False
         except IndexError:
            if (len(self.parm_data['POINTERS'])) < 30:
               print >> stderr, 'Error: Fewer integers in POINTERS section than expected! ' + \
                                'Likely a bad AMBER topology file.'
               self.valid = False

      if 'CTITLE' in self.flag_list: # We've found a chamber prmtop
         self.chamber = True

      if self.valid and 'LENNARD_JONES_ACOEF' in self.parm_data.keys() and 'LENNARD_JONES_BCOEF' in self.parm_data.keys():
         try:
            # fill LJ arrays with LJ data for easy manipulations
            self.fill_LJ()
            if self.chamber:
               self.fill_14_LJ()
         except Exception, err:
            print >> stderr, 'Warning: Problem parsing L-J 6-12 parameters. \n%s' % err

      if rst7_name != '':
         self.LoadRst7(rst7_name)

      # Load the structure arrays
      self._load_structure()
      
      # We now have the following instance arrays: All arrays are dynamic such that
      # removing an item propagates the indices if applicable. bond has angle/dihed analogs.
      # All non-dynamic lists have a check on any modification function to track if they
      # have been changed or not so we know whether we have to reload the data before
      # writing.
      #
      # atom_list          a dynamic list of all _Atom objects
      # bond_type_list     a dynamic list of all _BondType objects
      # bonds_inc_h        list of all bonds including hydrogen
      # bonds_without_h    list of all bonds without hydrogen
      # angle_type_list
      # angles_inc_h
      # angles_without_h
      # dihedral_type_list
      # dihedrals_inc_h
      # dihedrals_without_h

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def LoadPointers(self):
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
      # The next is probably only there for LES-prmtops
      try:
         self.pointers["NCOPY"] = self.parm_data["POINTERS"][NCOPY]
      except: pass

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _load_structure(self):
      """ 
      Loads all of the topology instance variables. This is necessary if we actually
      want to modify the topological layout of our system (like deleting atoms)
      """
      ##### First create our atoms #####
      self.atom_list = _AtomList(self)
      ##### Next create our list of bonds #####
      self.bond_type_list = _BondTypeList(self)
      self.bonds_inc_h, self.bonds_without_h = _TrackedList(), _TrackedList()
      # Array of bonds with hydrogen
      for i in range(self.ptr('nbonh')):
         self.bonds_inc_h.append(
              _Bond(self.atom_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i  ]/3],
                    self.atom_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i+1]/3],
                    self.bond_type_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i+2]-1]))
      # Array of bonds without hydrogen
      for i in range(self.ptr('mbona')):
         self.bonds_without_h.append(
              _Bond(self.atom_list[self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i  ]/3],
                    self.atom_list[self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+1]/3],
                    self.bond_type_list[self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+2]-1]))
      # We haven't changed yet...
      self.bonds_inc_h.changed = self.bonds_without_h.changed = False
      ##### Next create our list of angles #####
      self.angle_type_list = _AngleTypeList(self)
      self.angles_inc_h, self.angles_without_h = _TrackedList(), _TrackedList()
      # Array of angles with hydrogen
      for i in range(self.ptr('ntheth')):
         self.angles_inc_h.append(
              _Angle(self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i  ]/3],
                     self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+1]/3],
                     self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+2]/3],
                     self.angle_type_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+3]-1]))
      # Array of angles without hydrogen
      for i in range(self.ptr('mtheta')):
         self.angles_without_h.append(
              _Angle(self.atom_list[self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i  ]/3],
                     self.atom_list[self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+1]/3],
                     self.atom_list[self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+2]/3],
                     self.angle_type_list[self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+3]-1]))
      # We haven't changed yet
      self.angles_inc_h.changed = self.angles_without_h.changed = False
      ##### Next create our list of dihedrals #####
      self.dihedral_type_list = _DihedralTypeList(self)
      self.dihedrals_inc_h, self.dihedrals_without_h = _TrackedList(), _TrackedList()
      # Array of dihedrals with hydrogen
      for i in range(self.ptr('nphih')):
         signs = [1,1]
         if self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+2] < 0: signs[0] = -1
         if self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+3] < 0: signs[1] = -1
         self.dihedrals_inc_h.append(
              _Dihedral(self.atom_list[self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i  ]/3],
                        self.atom_list[self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+1]/3],
                        self.atom_list[abs(self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+2]/3)],
                        self.atom_list[abs(self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+3]/3)],
                        self.dihedral_type_list[self.parm_data['DIHEDRALS_INC_HYDROGEN'][5*i+4]-1],
                        signs))
      # Array of dihedrals without hydrogen
      for i in range(self.ptr('mphia')):
         signs = [1,1]
         if self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+2] < 0: signs[0] = -1
         if self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+3] < 0: signs[1] = -1
         self.dihedrals_without_h.append(
              _Dihedral(self.atom_list[self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i  ]/3],
                        self.atom_list[self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+1]/3],
                        self.atom_list[abs(self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+2]/3)],
                        self.atom_list[abs(self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+3]/3)],
                        self.dihedral_type_list[self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][5*i+4]-1],
                        signs))
      # We haven't changed yet
      self.dihedrals_inc_h.changed = self.dihedrals_without_h.changed = False

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
      self.parm_comments = {}
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
            self.parm_comments[current_flag] = []
            self.flag_list.append(current_flag)
            gathering_data = False

         elif prmlines[i][0:8] == '%COMMENT':
            self.parm_comments[current_flag].append(prmlines[i][9:].strip())

         elif prmlines[i][0:7] == '%FORMAT':
            self.formats[current_flag] = prmlines[i][8:len(prmlines[i].strip())-1]
            number_items_perline, size_item, dat_type, junk = _parseFormat(self.formats[current_flag])
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

      # Fill the residue container array
      self.residue_container = []
      for i in range(self.parm_data['POINTERS'][NRES]-1):
         curres = self.parm_data['RESIDUE_POINTER'][i] - 1
         nextres = self.parm_data['RESIDUE_POINTER'][i+1] - 1
         for j in range(curres, nextres):
            self.residue_container.append(i+1)
      for i in range(self.parm_data['RESIDUE_POINTER'][len(self.parm_data['RESIDUE_POINTER'])-1]-1,
                     self.parm_data['POINTERS'][NATOM]):
         self.residue_container.append(self.parm_data['POINTERS'][NRES])

      # Load the bond[] list, so each Atom # references a list of bonded partners. Note
      # that the index into bond[] begins atom indexing at 0, and each atom located in
      # the bonded list also starts from index 0. Start with bonds containing H. Also load
      # all of the angles (for the excluded list)
      self.bonds = [[] for i in range(self.parm_data['POINTERS'][NATOM])]
      self.angles = [[] for i in range(self.parm_data['POINTERS'][NATOM])]
      for i in range(self.parm_data['POINTERS'][NBONH]):
         at1 = self.parm_data['BONDS_INC_HYDROGEN'][3*i  ] / 3
         at2 = self.parm_data['BONDS_INC_HYDROGEN'][3*i+1] / 3
         self.bonds[at1].append(at2)
         self.bonds[at2].append(at1)
      # Now do bonds not including hydrogen
      for i in range(self.parm_data['POINTERS'][NBONA]):
         at1 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i  ] / 3
         at2 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+1] / 3
         self.bonds[at1].append(at2)
         self.bonds[at2].append(at1)
      # Now do angles including hydrogen
      for i in range(self.parm_data['POINTERS'][NTHETH]):
         at1 = self.parm_data['ANGLES_INC_HYDROGEN'][3*i  ] / 3
         at2 = self.parm_data['ANGLES_INC_HYDROGEN'][3*i+1] / 3
         at3 = self.parm_data['ANGLES_INC_HYDROGEN'][3*i+2] / 3
         if not at2 in self.bonds[at1]: self.angles[at1].append(at2)
         if not at3 in self.bonds[at1]: self.angles[at1].append(at3)
         if not at1 in self.bonds[at2]: self.angles[at2].append(at1)
         if not at3 in self.bonds[at2]: self.angles[at2].append(at3)
         if not at1 in self.bonds[at3]: self.angles[at3].append(at1)
         if not at2 in self.bonds[at3]: self.angles[at3].append(at2)
      # Now do angles without hydrogen
      for i in range(self.parm_data['POINTERS'][MTHETA]):
         at1 = self.parm_data['ANGLES_WITHOUT_HYDROGEN'][3*i  ] / 3
         at2 = self.parm_data['ANGLES_WITHOUT_HYDROGEN'][3*i+1] / 3
         at3 = self.parm_data['ANGLES_WITHOUT_HYDROGEN'][3*i+2] / 3
         if not at2 in self.bonds[at1]: self.angles[at1].append(at2)
         if not at3 in self.bonds[at1]: self.angles[at1].append(at3)
         if not at1 in self.bonds[at2]: self.angles[at2].append(at1)
         if not at3 in self.bonds[at2]: self.angles[at2].append(at3)
         if not at1 in self.bonds[at3]: self.angles[at3].append(at1)
         if not at2 in self.bonds[at3]: self.angles[at3].append(at2)

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
         # Insert any comments before the %FORMAT specifier
         for comment in self.parm_comments[flag]:
            new_prm.write('%%COMMENT %s\n' % comment)
         new_prm.write('%%FORMAT(%s)\n' % self.formats[flag])
         number_items_perline, size_item, dat_type, decnum = _parseFormat(self.formats[flag])
         line = ''
         num_items = 0
         if len(self.parm_data[flag]) == 0: # empty field...
            new_prm.write('\n')
            continue
         for j in range(len(self.parm_data[flag])): # write data in new_prm
            if dat_type == 'dec' and 'E' in self.formats[flag].upper():
               line += ('%%%s.%sE' % (size_item, decnum)) % self.parm_data[flag][j] 
            elif dat_type == 'dec' and 'F' in self.formats[flag].upper():
               line += ('%%%s.%sF' % (size_item, decnum)) % self.parm_data[flag][j] 
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

   def remake_parm(self):
      """ Re-fills the topology file arrays if we have changed the underlying structure """
      # First we need to change the pointer arrays
      self.parm_data['POINTERS'][NATOM] = len(self.atom_list)
      self.parm_data['POINTERS'][NBONH] = len(self.bonds_inc_h)
      self.parm_data['POINTERS'][MBONA] = len(self.bonds_without_h)
      self.parm_data['POINTERS'][NBONA] = len(self.bonds_without_h)
      self.parm_data['POINTERS'][NTHETH] = len(self.angles_inc_h)
      self.parm_data['POINTERS'][MTHETA] = len(self.angles_without_h)
      self.parm_data['POINTERS'][NTHETA] = len(self.angles_without_h)
      self.parm_data['POINTERS'][NPHIH] = len(self.dihedrals_inc_h)
      self.parm_data['POINTERS'][MPHIA] = len(self.dihedrals_without_h)
      self.parm_data['POINTERS'][NPHIA] = len(self.dihedrals_without_h)
      self.parm_data['POINTERS'][NUMBND] = len(self.bond_type_list)
      self.parm_data['POINTERS'][NUMANG] = len(self.angle_type_list)
      self.parm_data['POINTERS'][NPTRA] = len(self.dihedral_type_list)

      # Now we have to fill all of the arrays

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def frcmod(self, frcmod="frcmod"):
      """Prints an Frcmod file that contains every parameter found in prmtop"""
      from math import pi

      print >> stderr, "Warning: AmberParm.Frcmod() does not work for 10-12 non-bonded prmtops yet!"

      self.fill_LJ()

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

         # support variable 1-4 scaled prmtops and write out an frcmod that will always have 1-4 scaling info
         try:
            scee_scale = self.parm_data["SCEE_SCALE_FACTOR"][term-1]
            scnb_scale = self.parm_data["SCNB_SCALE_FACTOR"][term-1]
         except:
            scee_scale = 1.2
            scnb_scale = 2.0

         if atom4 < 0:
            dihedral = "%s %8.3f %8.3f %5.1f" % (dihedral, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         elif atom3 < 0: # if there's another term in the series
            dihedral = "%s %4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, -self.parm_data["DIHEDRAL_PERIODICITY"][term-1],
                           scee_scale, scnb_scale)
         else:
            dihedral = "%s %4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1],
                           scee_scale, scnb_scale)
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

         # support variable 1-4 scaled prmtops and write out an frcmod that will always have 1-4 scaling info
         try:
            scee_scale = self.parm_data["SCEE_SCALE_FACTOR"][term-1]
            scnb_scale = self.parm_data["SCNB_SCALE_FACTOR"][term-1]
         except:
            scee_scale = 1.2
            scnb_scale = 2.0

         if atom4 < 0:
            dihedral = "%s %8.3f %8.3f %5.1f" % (dihedral, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1])
         elif atom3 < 0: # if there's another term in the series
            dihedral = "%s %4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, -self.parm_data["DIHEDRAL_PERIODICITY"][term-1],
                           scee_scale, scnb_scale)
         else:
            dihedral = "%s %4i %8.3f %8.3f %5.1f    SCEE=%s SCNB=%s" % (dihedral, 1, self.parm_data["DIHEDRAL_FORCE_CONSTANT"][term-1],
                           self.parm_data["DIHEDRAL_PHASE"][term-1]*180/pi, self.parm_data["DIHEDRAL_PERIODICITY"][term-1],
                           scee_scale, scnb_scale)

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
                  if float(stored_dihedtypes[j][11:].split()[3]) < 0 and stored_dihedtypes[j][0:11] == unique_diheds[i]:
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

   def writeOFF(self, off_file='off.lib'):
      """ Writes an OFF file from all of the residues found in a prmtop """
      from chemistry.amber.residue import ToResidue
   
      file = open(off_file,'w',0)
   
      # keep track of all the residues we have to print to the OFF file
      residues = []
   
      # First create a Molecule object from the prmtop
      mol = self.ToMolecule()
   
      # Now loop through all of the residues in the Molecule object and add
      # unique ones to the list of residues to print
      for i in range(len(mol.residues)):
         res = ToResidue(mol, i)
         present = False
         for compres in residues:
            if res == compres:
               present = True
   
         if not present:
            residues.append(res)
      
      # Now that we have all of the residues that we need to add, put their names
      # in the header of the OFF file
      file.write('!!index array str\n')
      for res in residues:
         file.write(' "%s"\n' % res.name)
   
      # Now write the OFF strings to the file
      for res in residues:
         file.write(res.OFF())

      file.close()

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
         lj_index = self.parm_data["NONBONDED_PARM_INDEX"][
                     self.pointers["NTYPES"] * i + i] - 1
         if self.parm_data["LENNARD_JONES_ACOEF"][lj_index] < 1.0e-6:
            self.LJ_radius.append(0)
            self.LJ_depth.append(0)
         else:
            factor = 2 * self.parm_data["LENNARD_JONES_ACOEF"][lj_index] / self.parm_data["LENNARD_JONES_BCOEF"][lj_index]
            self.LJ_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_depth.append(self.parm_data["LENNARD_JONES_BCOEF"][lj_index] / 2 / factor)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def fill_14_LJ(self):
      """ Fills the LJ_14_radius, LJ_14_depth arrays with data (LJ_types is identical)
          from LENNARD_JONES_14_ACOEF and LENNARD_JONES_14_BCOEF sections of the 
          prmtop files, by undoing the canonical combining rules. """
      if not self.chamber:
         raise TypeError('fill_14_LJ() only valid on a chamber prmtop!')
      self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
      self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
      one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

      for i in range(self.pointers["NTYPES"]):
         lj_index = self.parm_data["NONBONDED_PARM_INDEX"][
                     self.pointers["NTYPES"] * i + i] - 1
         if self.parm_data["LENNARD_JONES_14_ACOEF"][lj_index] < 1.0e-6:
            self.LJ_14_radius.append(0)
            self.LJ_14_depth.append(0)
         else:
            factor = 2 * self.parm_data["LENNARD_JONES_14_ACOEF"][lj_index] / \
                     self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index]
            self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_14_depth.append(self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index] 
                                    / 2 / factor)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_LJ(self):
      """ Takes the values of the LJ_radius and LJ_depth arrays and recalculates the 
          LENNARD_JONES_A/BCOEF topology sections from the canonical combining rules. 
      """
      from math import sqrt

      for i in range(self.pointers["NTYPES"]):
         for j in range(i,self.pointers["NTYPES"]):
            index = self.parm_data['NONBONDED_PARM_INDEX'][self.ptr('ntypes')*i + j] - 1
            rij = self.LJ_radius[i] + self.LJ_radius[j]
            wdij = sqrt(self.LJ_depth[i] * self.LJ_depth[j])
            self.parm_data["LENNARD_JONES_ACOEF"][index] = wdij * pow(rij, 12)
            self.parm_data["LENNARD_JONES_BCOEF"][index] = 2 * wdij * pow(rij, 6)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_14_LJ(self):
      """ Takes the values of the LJ_radius and LJ_depth arrays and recalculates the 
          LENNARD_JONES_A/BCOEF topology sections from the canonical combining rules
          for the 1-4 LJ interactions (CHAMBER only)
      """
      if not self.chamber:
         raise TypeError('recalculate_14_LJ() requires a CHAMBER prmtop!')
      from math import sqrt

      for i in range(self.pointers["NTYPES"]):
         for j in range(i,self.pointers["NTYPES"]):
            index = self.parm_data['NONBONDED_PARM_INDEX'][self.ptr('ntypes')*i + j] - 1
            rij = self.LJ_14_radius[i] + self.LJ_14_radius[j]
            wdij = sqrt(self.LJ_14_depth[i] * self.LJ_14_depth[j])
            self.parm_data["LENNARD_JONES_14_ACOEF"][index] = wdij * pow(rij, 12)
            self.parm_data["LENNARD_JONES_14_BCOEF"][index] = 2 * wdij * pow(rij, 6)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def LoadRst7(self, filename):
      """ Loads coordinates into the AmberParm class """
      self.rst7 = rst7(filename)
      if not self.rst7.valid:
         return -1
      self.coords = self.rst7.coords
      self.hasvels = self.rst7.hasvels
      self.hasbox = self.rst7.hasbox
      if self.hasbox:
         self.box = self.rst7.box
      if self.hasvels:
         self.vels = self.rst7.vels

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def addFlag(self, flag_name, flag_format, num_items=0, comments=[], data=None):
      """ Adds a new flag with the given flag name and Fortran format string
          and initializes the array with the values given, or as an array of 0s
          of length num_items
      """
      self.flag_list.append(flag_name.upper())
      self.formats[flag_name.upper()] = flag_format
      if data:
         self.parm_data[flag_name.upper()] = list(data)
      else:
         if num_items == 0:
            raise exceptions.FlagError("If you do not supply prmtop data, num_items cannot be 0")
         self.parm_data[flag_name.upper()] = [0 for i in range(num_items)]
      if comments:
         if type(comments).__name__ == 'str': comments = [comments]
         elif type(comments).__name__ == 'tuple': comments = list(comments)
         elif type(comments).__name__ == 'list': pass
         else: raise TypeError('Comments are wrong type. Must be string, list, or tuple')
         self.parm_comments[flag_name.upper()] = comments
      else: self.parm_comments[flag_name.upper()] = []

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ToMolecule(self):
      """ Translates an amber system into a molecule format """
      from chemistry.molecule import Molecule
      from copy import copy

      all_bonds = []        # bond array in Molecule format
      residue_pointers = [] # residue pointers adjusted for index starting from 0
      elements = []         # which element each atom is
      radii = []

      # Set up initial, blank, bond array
      for i in range(self.pointers['NATOM']):
         all_bonds.append([])
      
      # Fill up bond arrays with bond partners excluding H atoms
      for i in range(self.pointers['MBONA']):
         atom1 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i  ]/3
         atom2 = self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+1]/3
         all_bonds[atom1].append(atom2)
         all_bonds[atom2].append(atom1)

      # Fill up bond arrays with bond partners including H atoms
      for i in range(self.pointers['NBONH']):
         atom1 = self.parm_data['BONDS_INC_HYDROGEN'][3*i  ]/3
         atom2 = self.parm_data['BONDS_INC_HYDROGEN'][3*i+1]/3
         all_bonds[atom1].append(atom2)
         all_bonds[atom2].append(atom1)

      # Sort bond arrays
      for i in range(len(all_bonds)):
         all_bonds[i].sort()

      # Adjust RESIDUE_POINTER for indexing from 0
      for i in range(len(self.parm_data['RESIDUE_POINTER'])):
         residue_pointers.append(self.parm_data['RESIDUE_POINTER'][i]-1)

      # Determine which element each atom is
      for i in range(self.pointers['NATOM']):
         elements.append(Element(self.parm_data['MASS'][i]))

      # Put together the title
      title = ''
      for i in range(len(self.parm_data['TITLE'])):
         title += self.parm_data['TITLE'][i]

      # Fill the VDW radii array
      self.fill_LJ()
      for i in range(self.pointers['NATOM']):
         radii.append(self.LJ_radius[self.LJ_types[self.parm_data['AMBER_ATOM_TYPE'][i]]-1])

      try:
         if self.valid and self.rst7.valid:
            return Molecule(atoms=copy(self.parm_data['ATOM_NAME']), atom_types=copy(self.parm_data['AMBER_ATOM_TYPE']),
                            charges=copy(self.parm_data['CHARGE']), residues=copy(self.parm_data['RESIDUE_LABEL']), 
                            bonds=all_bonds, residue_pointers=residue_pointers, coords=copy(self.coords),
                            elements=elements, title=title, radii=radii)
      except AttributeError: # in case no coordinates were loaded, use a dummy-list
         if self.valid:
            return Molecule(atoms=copy(self.parm_data['ATOM_NAME']), atom_types=copy(self.parm_data['AMBER_ATOM_TYPE']),
                            charges=copy(self.parm_data['CHARGE']), residues=copy(self.parm_data['RESIDUE_LABEL']), 
                            bonds=all_bonds, residue_pointers=residue_pointers, coords=list(range(self.pointers['NATOM']*3)),
                            elements=elements, title=title, radii=radii)
            

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class rst7:
   """ Amber input coordinate (or restart coordinate) file format """
   
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, filename):
      """ Initialize the inpcrd file """
      self.filename = filename
      self.valid = False

      try:
         self._read()
      except BaseException, err:
         raise(exceptions.ReadError('Error parsing coordinates from %s: %s' % (self.filename, err)))

      self.valid = True

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _read(self):
      """ Read in the coordinates from the file """
      restrt = open(self.filename, 'r')
      lines = restrt.readlines()

      # Load the title, number of atoms, and time
      self.title = lines[0].strip()
      self.natom = int(lines[1].strip().split()[0])
      self.coords = []
      self.vels = []
      try:
         self.time = float(lines[1].strip().split()[1])
      except IndexError:
         self.time = 0.0

      # Check to see if we have velocities or not and box or not
      if len(lines) == int(ceil(self.natom/2.0) + 2):
         self.hasbox = False
         self.hasvels = False
      if len(lines) == int(ceil(self.natom/2.0) + 3):
         self.hasbox = True
         self.hasvels = False
      if len(lines) == int(2*ceil(self.natom/2.0) + 2):
         self.hasbox = False
         self.hasvels = True
      if len(lines) == int(2*ceil(self.natom/2.0) + 3):
         self.hasbox = True
         self.hasvels = True

      startline = 2
      endline = startline + int(ceil(self.natom/2.0))
      # load the coordinates
      for i in range(startline,endline):
         x1 = float(lines[i][0 :12])
         y1 = float(lines[i][12:24])
         z1 = float(lines[i][24:36])
         try:
            x2 = float(lines[i][36:48])
            y2 = float(lines[i][48:60])
            z2 = float(lines[i][60:72])
            self.coords.extend([x1, y1, z1, x2, y2, z2])
         except ValueError:
            self.coords.extend([x1, y1, z1])
         
      startline += int(ceil(self.natom/2.0))
      # load the velocities
      if self.hasvels:
         endline = startline + int(ceil(self.natom/2.0))

         for i in range(startline, endline):
            x1 = float(lines[i][0 :12])
            y1 = float(lines[i][12:24])
            z1 = float(lines[i][24:36])
            try:
               x2 = float(lines[i][36:48])
               y2 = float(lines[i][48:60])
               z2 = float(lines[i][60:72])
               self.vels.extend([x1, y1, z1, x2, y2, z2])
            except ValueError:
               self.vels.extend([x1, y1, z1])

         startline += int(ceil(self.natom/2.0))
      # load the box information
      if self.hasbox:
         endline = startline + 1
         self.box = lines[startline].strip().split()
         self.box[0], self.box[1], self.box[2]  = float(self.box[0]), float(self.box[1]), float(self.box[2])
         self.box[3], self.box[4], self.box[5]  = float(self.box[3]), float(self.box[4]), float(self.box[5])

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Element(mass):
   """ Determines what element the given atom is based on its mass """

   diff = mass
   best_guess = 'EP'

   for element in periodic_table.Element:
      if abs(periodic_table.Mass[element] - mass) < diff:
         best_guess = element
         diff = abs(periodic_table.Mass[element] - mass)

   return best_guess

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# For backwards-compatibility
class amberParm(AmberParm):
   """ 
   This equivalence is made to preserve backwards-compatibility. The standard for
   creating class names is CapitalLettersSeparateWords, with the starting letter
   being capital.
   """
