"""
This module contains an amber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.

Copyright (C) 2010 - 2012  Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

from chemistry import periodic_table, __version__
from chemistry.amber.topologyobjects import (Bond, Angle,
            Dihedral, ResidueList, AtomList, BondTypeList,
            AngleTypeList, DihedralTypeList, TrackedList)
from chemistry.amber.constants import (NATOM, NTYPES, NBONH, MBONA, NTHETH,
            MTHETA, NPHIH, MPHIA, NHPARM, NPARM, NEXT, NRES, NBONA, NTHETA,
            NPHIA, NUMBND, NUMANG, NPTRA, NATYP, NPHB, IFPERT, NBPER, NGPER,
            NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, NUMEXTRA, NCOPY,
            NNB)
from chemistry.amber.amberformat import AmberFormat
from chemistry.exceptions import (AmberParmWarning, AmberParmError, ReadError,
                                  MoleculeError, MoleculeWarning)
from warnings import warn
from math import ceil, sqrt

# Import NetCDF support either from Scientific, netCDF4, or pynetcdf. Set up
# API-independent reading and writing functions
_HAS_NETCDF = True

try:
   from Scientific.IO.NetCDF import NetCDFFile as NetCDFFile
   open_netcdf = lambda name, mode: NetCDFFile(name, mode)
   get_int_dimension = lambda obj, name: obj.dimensions[name]
   get_float_array = lambda obj, name: obj.variables[name].getValue()
   get_float = lambda obj, name: obj.variables[name].getValue()
except ImportError:
   try:
      from netCDF4 import Dataset as NetCDFFile
      open_netcdf = lambda name, mode: NetCDFFile(name, mode, format='NETCDF3_64BIT')
      get_int_dimension = lambda obj, name: len(obj.dimensions[name])
      get_float_array = lambda obj, name: obj.variables[name][:]
      get_float = lambda obj, name: obj.variables[name].getValue()[0]
   except ImportError:
      try:
         from pynetcdf import NetCDFFile
         open_netcdf = lambda name, mode: NetCDFFile(name, mode)
         get_int_dimension = lambda obj, name: obj.dimensions[name]
         get_float_array = lambda obj, name: obj.variables[name].getValue()
         get_float = lambda obj, name: obj.variables[name].getValue()
      except ImportError:
         _HAS_NETCDF = False

# Since numpy is a prereq for all NetCDF packages, this is OK
if _HAS_NETCDF:
   import numpy as np
else:
   np = None # to avoid NameError's

class AmberParm(AmberFormat):
   """
   Amber Topology (parm7 format) class. Gives low, and some high, level access
   to topology data.
   """

   solvent_residues = ['WAT', 'HOH']

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, prm_name=None, rst7_name=None):
      """
      Instantiates an AmberParm object from data in prm_name and establishes
      validity based on presence of POINTERS and CHARGE sections. During
      instantiation, the following functions are called:
         o  rdparm()  -- Reads the top file (calls rdparm_old() if necessary)
                         rdparm() scales the charges to electron-charges
         o  LoadPointers() -- Loads the POINTERS section into the pointers dict
         o  fill_LJ() -- Fills the Lennard Jones radius/depth arrays from
                         diagonal terms
         o  fill_14_LJ() -- Fills the 1-4 Lennard Jones radius/depth arrays for
                            chamber topology files
         o  _load_structure() -- Fills data structures (atom/bond/angle/dihedral
                                 lists, etc.) to allow rebuilding the topology
      """

      AmberFormat.__init__(self, prm_name)

      # instance variables other than those in AmberFormat
      self.overwrite = False   # whether writeParm will overwrite existing file
      self.chamber = False     # Is this a chamber-generated prmtop file?
      self.pointers = {}       # list of all the pointers in the prmtop
      self.LJ_types = {}       # dict pairing atom name with its LJ atom type #
      self.LJ_radius = []      # ordered array of L-J radii in Ang -- indices
                               # are elements in LJ_types-1
      self.LJ_depth = []       # similarly ordered array of L-J depths

      if self.prm_name is None: return

      # If we were given a prmtop, read it in
      if self.valid:
         try:
            self.LoadPointers()
            self.valid = True
         except KeyError:
            warn('POINTERS flag not found! Likely a bad AMBER topology file.',
                 AmberParmWarning)
            self.valid = False
         except IndexError:
            if (len(self.parm_data['POINTERS'])) < 30:
               warn('Fewer integers in POINTERS section than expected! Likely '
                    'a bad AMBER topology file.', AmberParmWarning)
               self.valid = False

      if 'CTITLE' in self.flag_list: # We've found a chamber prmtop
         self.chamber = True

      if self.valid and 'LENNARD_JONES_ACOEF' in self.parm_data.keys() and \
                              'LENNARD_JONES_BCOEF' in self.parm_data.keys():
         try:
            # fill LJ arrays with LJ data for easy manipulations
            self.fill_LJ()
            if self.chamber:
               self.fill_14_LJ()
         except Exception, err:
            warn('Problem parsing L-J 6-12 parameters. \n%s' % err,
                 AmberParmWarning)

      # Load the structure arrays
      if self.valid and not self.chamber: 
         try:
            self._load_structure()
         except (KeyError, IndexError, AttributeError):
            warn('Error loading molecule topology. Cannot use delete_mask',
                 AmberParmWarning)
            self.valid = False
      
      # We now have the following instance arrays: All arrays are dynamic such
      # that removing an item propagates the indices if applicable. bond has
      # angle/dihed analogs. All non-dynamic lists have a check on any
      # modification function to track if they have been changed or not so we
      # know whether we have to reload the data before writing.
      #
      # atom_list          a dynamic list of all _Atom objects
      # residue_list       a dynamic list of all _Residue objects
      # bond_type_list     a dynamic list of all _BondType objects
      # bonds_inc_h        list of all bonds including hydrogen
      # bonds_without_h    list of all bonds without hydrogen
      # angle_type_list
      # angles_inc_h
      # angles_without_h
      # dihedral_type_list
      # dihedrals_inc_h
      # dihedrals_without_h

      if rst7_name is not None:
         self.LoadRst7(rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def LoadPointers(self):
      """
      Loads the data in POINTERS section into a pointers dictionary with each
      key being the pointer name according to http://ambermd.org/formats.html
      """
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _load_structure(self):
      """ 
      Loads all of the topology instance variables. This is necessary if we
      actually want to modify the topological layout of our system
      (like deleting atoms)
      """
      ##### First create our atoms #####
      self.atom_list = AtomList(self)
      ##### Next, load our residues #####
      self.residue_list = ResidueList(self)
      for i in range(self.ptr('natom')):
         residx = self.residue_container[i] - 1
         self.atom_list[i].residue = self.residue_list[residx]
      ##### Next create our list of bonds #####
      self.bond_type_list = BondTypeList(self)
      self.bonds_inc_h, self.bonds_without_h = TrackedList(), TrackedList()
      # Array of bonds with hydrogen
      for i in range(self.ptr('nbonh')):
         self.bonds_inc_h.append(
           Bond(self.atom_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i  ]/3],
            self.atom_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i+1]/3],
            self.bond_type_list[self.parm_data['BONDS_INC_HYDROGEN'][3*i+2]-1]))
      # Array of bonds without hydrogen
      for i in range(self.ptr('mbona')):
         self.bonds_without_h.append(
           Bond(self.atom_list[
                        self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i  ]/3],
                 self.atom_list[
                        self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+1]/3],
                 self.bond_type_list[
                        self.parm_data['BONDS_WITHOUT_HYDROGEN'][3*i+2]-1]))
      # We haven't changed yet...
      self.bonds_inc_h.changed = self.bonds_without_h.changed = False
      ##### Next create our list of angles #####
      self.angle_type_list = AngleTypeList(self)
      self.angles_inc_h, self.angles_without_h = TrackedList(), TrackedList()
      # Array of angles with hydrogen
      for i in range(self.ptr('ntheth')):
         self.angles_inc_h.append(
          Angle(self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i  ]/3],
          self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+1]/3],
          self.atom_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+2]/3],
          self.angle_type_list[self.parm_data['ANGLES_INC_HYDROGEN'][4*i+3]-1]))
      # Array of angles without hydrogen
      for i in range(self.ptr('mtheta')):
         self.angles_without_h.append(
          Angle(self.atom_list[
                      self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i  ]/3],
          self.atom_list[
                      self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+1]/3],
          self.atom_list[
                      self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+2]/3],
          self.angle_type_list[
                      self.parm_data['ANGLES_WITHOUT_HYDROGEN'][4*i+3]-1]))
      # We haven't changed yet
      self.angles_inc_h.changed = self.angles_without_h.changed = False
      ##### Next create our list of dihedrals #####
      self.dihedral_type_list = DihedralTypeList(self)
      self.dihedrals_inc_h = TrackedList()
      self.dihedrals_without_h = TrackedList()
      # Array of dihedrals with hydrogen
      flag = 'DIHEDRALS_INC_HYDROGEN'
      for i in range(self.ptr('nphih')):
         signs = [1,1]
         if self.parm_data[flag][5*i+2] < 0: signs[0] = -1
         if self.parm_data[flag][5*i+3] < 0: signs[1] = -1
         self.dihedrals_inc_h.append(
          Dihedral(self.atom_list[self.parm_data[flag][5*i  ]/3],
                    self.atom_list[self.parm_data[flag][5*i+1]/3],
                    self.atom_list[abs(self.parm_data[flag][5*i+2]/3)],
                    self.atom_list[abs(self.parm_data[flag][5*i+3]/3)],
                    self.dihedral_type_list[self.parm_data[flag][5*i+4]-1],
                    signs))
      # Array of dihedrals without hydrogen
      flag = 'DIHEDRALS_WITHOUT_HYDROGEN'
      for i in range(self.ptr('mphia')):
         signs = [1,1]
         if self.parm_data[flag][5*i+2] < 0: signs[0] = -1
         if self.parm_data[flag][5*i+3] < 0: signs[1] = -1
         self.dihedrals_without_h.append(
              Dihedral(self.atom_list[self.parm_data[flag][5*i  ]/3],
                        self.atom_list[self.parm_data[flag][5*i+1]/3],
                        self.atom_list[abs(self.parm_data[flag][5*i+2]/3)],
                        self.atom_list[abs(self.parm_data[flag][5*i+3]/3)],
                        self.dihedral_type_list[self.parm_data[flag][5*i+4]-1],
                        signs))
      # We haven't changed yet
      self.dihedrals_inc_h.changed = self.dihedrals_without_h.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __str__(self):
      """ Returns the name of the topology file as its string representation """
      return self.prm_name

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ptr(self,pointer):
      """
      Returns the value of the given pointer, and converts to upper-case so it's
      case-insensitive. A non-existent pointer meets with a KeyError
      """
      global AMBER_POINTERS
      return self.pointers[pointer.upper()]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rdparm(self, fname=None):
      """ Reads the topology file and loads data in parm_data array """
      
      AmberFormat.rdparm(self, fname)

      # If we don't have a version, then read in an old-file topology
      if self.version is None:
         self.rdparm_old(open(self.prm_name, 'r').readlines())

      # Fill the residue container array
      self._fill_res_container()

      # Load the bond[] list, so each Atom # references a list of bonded
      # partners. Note that the index into bond[] begins atom indexing at 0,
      # and each atom located in the bonded list also starts from index 0.
      # Start with bonds containing H. Also load all of the angles
      # (for the excluded list)
      self.bonds = [[] for i in range(self.parm_data['POINTERS'][NATOM])]
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rdparm_old(self, prmtop_lines):
      """
      This reads an old-style topology file and stores the results in the same
      data structures as a new-style topology file
      """
      def read_integer(line_idx, lines, num_items):
         # line_idx should be the line _before_ the first line you
         # want data from.
         i, tmp_data = 0, []
         while i < num_items:
            idx = i % 12
            if idx == 0:
               line_idx += 1
            try:
               tmp_data.append(int(lines[line_idx][idx*6:idx*6+6]))
            except ValueError, err:
               print 'Error parsing line (', line_idx, '):', \
                        lines[line_idx], idx*6, idx*6+6, num_items, i
               raise err
            i += 1
         # If we had no items, we need to jump a line:
         if num_items == 0: line_idx += 1
         return tmp_data, line_idx

      def read_string(line_idx, lines, num_items):
         # line_idx should be the line _before_ the first line you
         # want data from.
         i, tmp_data = 0, []
         while i < num_items:
            idx = i % 20
            if idx == 0:
               line_idx += 1
            tmp_data.append(lines[line_idx][idx*4:idx*4+4])
            i += 1
         # If we had no items, we need to jump a line:
         if num_items == 0: line_idx += 1
         return tmp_data, line_idx

      def read_float(line_idx, lines, num_items):
         # line_idx should be the line _before_ the first line you
         # want data from.
         i, tmp_data = 0, []
         while i < num_items:
            idx = i % 5
            if idx == 0:
               line_idx += 1
            try:
               tmp_data.append(float(lines[line_idx][idx*16:idx*16+16]))
            except ValueError, err:
               print 'Error parsing line (', line_idx, '):', \
                        lines[line_idx], idx*6, idx*6+6, num_items, i
               raise err
            i += 1
         # If we had no items, we need to jump a line:
         if num_items == 0: line_idx += 1
         return tmp_data, line_idx

      # First add a title
      self.addFlag('TITLE', '20a4', data=['| Converted old-style topology'])

      # Next, read in the pointers
      line_idx = 0
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 30)
      # Add a final pointer of 0, which corresponds to NUMEXTRA
      tmp_data.append(0)
      self.addFlag('POINTERS', '10I8', data=tmp_data)

      # Next read in the atom names
      tmp_data, line_idx = read_string(line_idx, prmtop_lines,
                                       self.parm_data['POINTERS'][NATOM])
      self.addFlag('ATOM_NAME', '20a4', data=tmp_data)

      # Next read the charges
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NATOM])
      self.addFlag('CHARGE', '5E16.8', data=tmp_data)

      # Next read the masses
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NATOM])
      self.addFlag('MASS', '5E16.8', data=tmp_data)

      # Next read atom type index
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NATOM])
      self.addFlag('ATOM_TYPE_INDEX', '10I8', data=tmp_data)

      # Next read number excluded atoms
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NATOM])
      self.addFlag('NUMBER_EXCLUDED_ATOMS', '10I8', data=tmp_data)
      
      # Next read nonbonded parm index
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NTYPES]**2)
      self.addFlag('NONBONDED_PARM_INDEX', '10I8', data=tmp_data)

      # Next read residue label
      tmp_data, line_idx = read_string(line_idx, prmtop_lines,
                                       self.parm_data['POINTERS'][NRES])
      self.addFlag('RESIDUE_LABEL', '20a4', data=tmp_data)

      # Next read residue pointer
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NRES])
      self.addFlag('RESIDUE_POINTER', '10I8', data=tmp_data)
   
      # Next read bond force constant
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NUMBND])
      self.addFlag('BOND_FORCE_CONSTANT', '5E16.8', data=tmp_data)

      # Next read bond equil value
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NUMBND])
      self.addFlag('BOND_EQUIL_VALUE', '5E16.8', data=tmp_data)

      # Next read angle force constant
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NUMANG])
      self.addFlag('ANGLE_FORCE_CONSTANT', '5E16.8', data=tmp_data)

      # Next read the angle equilibrium value
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NUMANG])
      self.addFlag('ANGLE_EQUIL_VALUE', '5E16.8', data=tmp_data)

      # Next read the dihedral force constant
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPTRA])
      self.addFlag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', data=tmp_data)

      # Next read dihedral periodicity
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPTRA])
      self.addFlag('DIHEDRAL_PERIODICITY', '5E16.8', data=tmp_data)

      # Next read the dihedral phase
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPTRA])
      self.addFlag('DIHEDRAL_PHASE', '5E16.8', data=tmp_data)

      # Next read SOLTY (?)
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NATYP])
      self.addFlag('SOLTY', '5E16.8', data=tmp_data)

      # Next read lennard jones acoef and bcoef
      numvals = self.parm_data['POINTERS'][NTYPES]
      numvals = numvals * (numvals + 1) / 2
      tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
      self.addFlag('LENNARD_JONES_ACOEF', '5E16.8', data=tmp_data)
      tmp_data, line_idx = read_float(line_idx, prmtop_lines, numvals)
      self.addFlag('LENNARD_JONES_BCOEF', '5E16.8', data=tmp_data)

      # Next read bonds including hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NBONH]*3)
      self.addFlag('BONDS_INC_HYDROGEN', '10I8', data=tmp_data)

      # Next read bonds without hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NBONA]*3)
      self.addFlag('BONDS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

      # Next read angles including hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NTHETH]*4)
      self.addFlag('ANGLES_INC_HYDROGEN', '10I8', data=tmp_data)

      # Next read angles without hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NTHETA]*4)
      self.addFlag('ANGLES_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

      # Next read dihdrals including hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NPHIH]*5)
      self.addFlag('DIHEDRALS_INC_HYDROGEN', '10I8', data=tmp_data)

      # Next read dihedrals without hydrogen
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NPHIA]*5)
      self.addFlag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8', data=tmp_data)

      # Next read the excluded atoms list
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NEXT])
      self.addFlag('EXCLUDED_ATOMS_LIST', '10I8', data=tmp_data)

      # Next read the hbond terms
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPHB])
      self.addFlag('HBOND_ACOEF', '5E16.8', data=tmp_data)
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPHB])
      self.addFlag('HBOND_BCOEF', '5E16.8', data=tmp_data)
      tmp_data, line_idx = read_float(line_idx, prmtop_lines,
                                      self.parm_data['POINTERS'][NPHB])
      self.addFlag('HBCUT', '5E16.8', data=tmp_data)

      # Next read amber atom type
      tmp_data, line_idx = read_string(line_idx, prmtop_lines,
                                       self.parm_data['POINTERS'][NATOM])
      self.addFlag('AMBER_ATOM_TYPE', '20a4', data=tmp_data)

      # Next read tree chain classification
      tmp_data, line_idx = read_string(line_idx, prmtop_lines,
                                       self.parm_data['POINTERS'][NATOM])
      self.addFlag('TREE_CHAIN_CLASSIFICATION', '20a4', data=tmp_data)

      # Next read the join array
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NATOM])
      self.addFlag('JOIN_ARRAY', '10I8', data=tmp_data)

      # Next read the irotat array
      tmp_data, line_idx = read_integer(line_idx, prmtop_lines,
                                        self.parm_data['POINTERS'][NATOM])
      self.addFlag('IROTAT', '10I8', data=tmp_data)

      # Now do PBC stuff
      if self.parm_data['POINTERS'][IFBOX]:
         # Solvent pointers
         tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 3)
         self.addFlag('SOLVENT_POINTERS', '10I8', data=tmp_data)

         # Atoms per molecule
         tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 
                                          self.parm_data['SOLVENT_POINTERS'][1])
         self.addFlag('ATOMS_PER_MOLECULE', '10I8', data=tmp_data)

         # Box dimensions
         tmp_data, line_idx = read_float(line_idx, prmtop_lines, 4)
         self.addFlag('BOX_DIMENSIONS', '5E16.8', data=tmp_data)

      # end if self.parm_data['POINTERS'][IFBOX]

      # Now do CAP stuff
      if self.parm_data['POINTERS'][IFCAP]:
         # CAP_INFO
         tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 1)
         self.addFlag('CAP_INFO', '10I8', data=tmp_data)
         tmp_data, line_idx = read_integer(line_idx, prmtop_lines, 4)
         self.addFlag('CAP_INFO2', '10I8', data=tmp_data)
      # end if self.parm_data['POINTERS'][IFCAP]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeRst7(self, name, netcdf=None):
      """ Writes a restart file with the current coordinates and velocities
          and box info if it's present
      """
      # By default, determine file type by extension (.ncrst is NetCDF)
      if netcdf or (netcdf is None and name.endswith('.ncrst')):
         return self._write_ncrst(name)

      # If we are here, continue with the ASCII format
      restrt = open(name, 'w')
      restrt.write('Restart file written by amberParm\n')
      if len(self.atom_list) < 100000:
         restrt.write('%5d%15.7e\n' % (len(self.atom_list), self.rst7.time))
      elif len(self.atom_list) < 1000000:
         restrt.write('%6d%15.7e\n' % (len(self.atom_list), self.rst7.time))
      elif len(self.atom_list) < 10000000:
         restrt.write('%7d%15.7e\n' % (len(self.atom_list), self.rst7.time))
      else:
         restrt.write('%8d%15e7\n' % (len(self.atom_list), self.rst7.time))
      # Write the coordinates
      num_writ = 0
      for i in range(len(self.atom_list)):
         restrt.write('%12.7f%12.7f%12.7f' %
                     (self.atom_list[i].xx, self.atom_list[i].xy,
                      self.atom_list[i].xz))
         num_writ += 1
         if num_writ % 2 == 0: restrt.write('\n')
      if num_writ % 2 != 0: restrt.write('\n')
      # Write the velocities if they exist
      if self.hasvels:
         num_writ = 0
         for i in range(len(self.atom_list)):
            restrt.write('%12.7f%12.7f%12.7f' %
                        (self.atom_list[i].vx, self.atom_list[i].vy,
                         self.atom_list[i].vz))
            num_writ += 1
            if num_writ % 2 == 0: restrt.write('\n')
         if num_writ % 2 != 0: restrt.write('\n')
         
      if self.hasbox:
         for num in self.box: restrt.write('%12.7f' % num)
         restrt.write('\n')
   
      restrt.close()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _write_ncrst(self, name):
      """ Write a restart file in NetCDF format """
      global np, open_netcdf
      restrt = open_netcdf(name, 'w')
      # Assign the main attributes
      restrt.Conventions = 'AMBERRESTART'
      restrt.ConventionVersion = "1.0"
      restrt.title = self.rst7.title
      restrt.application = "AMBER"
      restrt.program = "ParmEd"
      restrt.programVersion = str(__version__)
      # Make all of the dimensions
      restrt.createDimension('spatial', 3)
      restrt.createDimension('atom', len(self.atom_list))
      if self.rst7.hasbox:
         restrt.createDimension('cell_spatial', 3)
         restrt.createDimension('label', 5)
         restrt.createDimension('cell_angular', 3)
      else:
         restrt.createDimension('label', 5)
      restrt.createDimension('time', 1)
      # Now make the variables
      v = restrt.createVariable('time', 'd', ('time',))
      v.units = 'picosecond'
      v[0] = self.rst7.time
      v = restrt.createVariable('spatial', 'c', ('spatial',))
      v[:] = np.asarray(list('xyz'))
      if self.rst7.hasbox:
         v = restrt.createVariable('cell_angular', 'c', ('cell_angular', 'label'))
         v[0] = np.asarray(list('alpha'))
         v[1] = np.asarray(list('beta '))
         v[2] = np.asarray(list('gamma'))
         v = restrt.createVariable('cell_spatial', 'c', ('cell_spatial',))
         v[0], v[1], v[2] = 'a', 'b', 'c'
         v = restrt.createVariable('cell_lengths', 'd', ('cell_spatial',))
         v.units = 'angstrom'
         v[:] = np.asarray(self.rst7.box[:3])
         v = restrt.createVariable('cell_angles', 'd', ('cell_angular',))
         v.units = 'degree'
         v[:] = np.asarray(self.rst7.box[3:])
      v = restrt.createVariable('coordinates', 'd', ('atom', 'spatial'))

      v.units = 'angstrom'
      v[:] = np.asarray([[a.xx, a.xy, a.xz] for a in self.atom_list])
      if self.rst7.hasvels:
         v = restrt.createVariable('velocities', 'd', ('atom', 'spatial'))
         v.units = 'angstrom/picosecond'
         v.scale_factor = np.float32(20.455)
         v[:] = np.asarray([[a.vx, a.vy, a.vz] for a in self.atom_list])
      # Now we're done
      restrt.close()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeParm(self, name):   # write a new prmtop with the current prmtop data
      """
      Writes the current data in parm_data into a new topology file with a given
      name. Will not overwrite the original prm_name unless the overwrite
      variable is set to True.
      """
      if hasattr(self, 'atom_list') and self._topology_changed(): 
         self.remake_parm()
         # Reload the structure now that we've recalculated it to flush all data
         # structures to what they *should* be.
         self._load_structure()
         # Now we have to redo the ATOMS_PER_MOLECULE/SOLVENT_POINTERS sections
         if self.ptr('ifbox'): self.rediscover_molecules()

      AmberFormat.writeParm(self, name, self.overwrite)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def remake_parm(self):
      """
      Re-fills the topology file arrays if we have changed the underlying
      structure
      """
      # First thing we have to do is load any of our old atom parameters into
      # our atom_list to preserve any changes we've made directly to the data
      self.atom_list.refresh_data()
      # Fill up the atom arrays. This will also adjust NATOM for us if 
      # we've deleted atoms
      self.atom_list.write_to_parm()

      # Recount number of extra points
      nextra = 0
      for atm in self.atom_list:
         if atm.attype[:2] in ['EP', 'LP']: nextra += 1
      self.parm_data['POINTERS'][NUMEXTRA] = nextra

      self.parm_data['POINTERS'][NNB] = len(self.parm_data['EXCLUDED_ATOMS_LIST'])
      # Write the residue arrays
      num_res = 0
      for i in range(len(self.atom_list)):
         atm = self.atom_list[i]
         if atm.residue.idx == -1:
            self.parm_data['RESIDUE_LABEL'][num_res] = atm.residue.resname
            self.parm_data['RESIDUE_POINTER'][num_res] = i+1
            atm.residue.idx = num_res
            num_res += 1
      self.parm_data['POINTERS'][NRES] = num_res
      self.parm_data['RESIDUE_LABEL']=self.parm_data['RESIDUE_LABEL'][:num_res]
      self.parm_data['RESIDUE_POINTER'] = self.parm_data['RESIDUE_POINTER'][
                                                                      :num_res]
      self._fill_res_container()

      # Adjust NMXRS (number of atoms in largest residue) in case that changed
      max_atm = 0
      for i in range(1,self.parm_data['POINTERS'][NRES]):
         atm_in_res = self.parm_data['RESIDUE_POINTER'][i] \
                    - self.parm_data['RESIDUE_POINTER'][i-1]
         if atm_in_res > max_atm: max_atm = atm_in_res

      # Check the last residue
      max_atm = max(self.parm_data['POINTERS'][NATOM] -
        self.parm_data['RESIDUE_POINTER'][self.parm_data['POINTERS'][NRES]-1]+1,
        max_atm)
      
      self.parm_data['POINTERS'][NMXRS] = max_atm

      # Now write all of the bond arrays. We will loop through all of the bonds
      # to make sure that all of their atoms still exist (atm.idx > -1). At the
      # same time, we will start applying indexes to the bond_types so we only
      # print out the bond types that will be used. To do this, we need a couple
      # counters. Different bond types will have an index of -1 until we find
      # out they are needed. Then we assign them an index and write out that
      # bond info. We also have to make sure that every array is at least large
      # enough, so give it enough elements to cover every bond in the list,
      # which will be reduced in size if not every bond is actually added
      bond_num = 0
      bond_type_num = 0
      self.parm_data['BONDS_INC_HYDROGEN'] = [0 for i in 
                                              range(len(self.bonds_inc_h)*3)]
      for i in range(len(self.bonds_inc_h)):
         holder = self.bonds_inc_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx): continue
         if holder.bond_type.idx == -1:
            holder.bond_type.idx = bond_type_num
            bond_type_num += 1
         holder.write_info(self, 'BONDS_INC_HYDROGEN', bond_num)
         bond_num += 1
      self.parm_data['POINTERS'][NBONH] = bond_num
      self.parm_data['BONDS_INC_HYDROGEN'] = \
               self.parm_data['BONDS_INC_HYDROGEN'][:3*bond_num]
      # Now we know how many bonds with hydrogen we have. Note that bond_num is
      # +1 past the last index used, but that last index is -1 from total bond #
      # due to indexing from 0, so it's just right now. So is the bond_type
      # index, but that is applicable for the bonds_without_h as well.
      bond_num = 0
      self.parm_data['BONDS_WITHOUT_HYDROGEN'] = \
               [0 for i in range(len(self.bonds_without_h)*3)]
      for i in range(len(self.bonds_without_h)):
         holder = self.bonds_without_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx): continue
         if holder.bond_type.idx == -1:
            holder.bond_type.idx = bond_type_num
            bond_type_num += 1
         holder.write_info(self, 'BONDS_WITHOUT_HYDROGEN', bond_num)
         bond_num += 1
      # Make sure BOND_FORCE_CONSTANT and BOND_EQUIL_VALUE arrays are big enough
      self.parm_data['BOND_FORCE_CONSTANT'] = [0 for i in range(bond_type_num)]
      self.parm_data['BOND_EQUIL_VALUE'] = [0 for i in range(bond_type_num)]
      # Now we can write all of the bond types out
      self.bond_type_list.write_to_parm()
      # Now we know how many bonds without H we have and our # of bond types
      self.parm_data['POINTERS'][MBONA] = bond_num
      self.parm_data['POINTERS'][NBONA] = bond_num
      self.parm_data['POINTERS'][NUMBND] = bond_type_num
      self.parm_data['BONDS_WITHOUT_HYDROGEN'] = \
               self.parm_data['BONDS_WITHOUT_HYDROGEN'][:3*bond_num]

      # Now do all of the angle arrays
      angle_num = 0
      angle_type_num = 0
      # Make sure we have enough ANGLES_INC_HYDROGEN
      self.parm_data['ANGLES_INC_HYDROGEN'] = [0 
                                       for i in range(len(self.angles_inc_h)*4)]
      for i in range(len(self.angles_inc_h)):
         holder = self.angles_inc_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx, holder.atom3.idx):
            continue
         if holder.angle_type.idx == -1:
            holder.angle_type.idx = angle_type_num
            angle_type_num += 1
         holder.write_info(self, 'ANGLES_INC_HYDROGEN', angle_num)
         angle_num += 1
      self.parm_data['POINTERS'][NTHETH] = angle_num
      self.parm_data['ANGLES_INC_HYDROGEN'] = \
               self.parm_data['ANGLES_INC_HYDROGEN'][:4*angle_num]
      # Time for Angles without H
      angle_num = 0
      self.parm_data['ANGLES_WITHOUT_HYDROGEN'] = \
               [0 for i in range(len(self.angles_without_h)*4)]
      for i in range(len(self.angles_without_h)):
         holder = self.angles_without_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx, holder.atom3.idx):
            continue
         if holder.angle_type.idx == -1:
            holder.angle_type.idx = angle_type_num
            angle_type_num += 1
         holder.write_info(self, 'ANGLES_WITHOUT_HYDROGEN', angle_num)
         angle_num += 1
      # Make sure BOND_FORCE_CONSTANT and BOND_EQUIL_VALUE arrays are big enough
      self.parm_data['ANGLE_FORCE_CONSTANT']=[0 for i in range(angle_type_num)]
      self.parm_data['ANGLE_EQUIL_VALUE']=[0 for i in range(angle_type_num)]
      # Write angle type info to parm
      self.angle_type_list.write_to_parm()
      self.parm_data['POINTERS'][NTHETA] = angle_num
      self.parm_data['POINTERS'][MTHETA] = angle_num
      self.parm_data['POINTERS'][NUMANG] = angle_type_num
      self.parm_data['ANGLES_WITHOUT_HYDROGEN'] = \
               self.parm_data['ANGLES_WITHOUT_HYDROGEN'][:4*angle_num]

      # Now do all of the dihedral arrays
      dihedral_num = 0
      dihedral_type_num = 0
      self.parm_data['DIHEDRALS_INC_HYDROGEN'] = [0 for i in 
                                             range(len(self.dihedrals_inc_h)*5)]
      for i in range(len(self.dihedrals_inc_h)):
         holder = self.dihedrals_inc_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx,
                   holder.atom3.idx, holder.atom4.idx):
            continue
         if holder.dihed_type.idx == -1:
            holder.dihed_type.idx = dihedral_type_num
            dihedral_type_num += 1
         holder.write_info(self, 'DIHEDRALS_INC_HYDROGEN', dihedral_num)
         dihedral_num += 1
      self.parm_data['POINTERS'][NPHIH] = dihedral_num
      self.parm_data['DIHEDRALS_INC_HYDROGEN'] = \
               self.parm_data['DIHEDRALS_INC_HYDROGEN'][:5*dihedral_num]
      # Time for dihedrals without H
      dihedral_num = 0
      self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'] = \
                  [0 for i in range(len(self.dihedrals_without_h)*5)]
      for i in range(len(self.dihedrals_without_h)):
         holder = self.dihedrals_without_h[i]
         if -1 in (holder.atom1.idx, holder.atom2.idx,
                   holder.atom3.idx, holder.atom4.idx):
            continue
         if holder.dihed_type.idx == -1:
            holder.dihed_type.idx = dihedral_type_num
            dihedral_type_num += 1
         holder.write_info(self, 'DIHEDRALS_WITHOUT_HYDROGEN', dihedral_num)
         dihedral_num += 1
      self.parm_data['POINTERS'][NPHIA] = dihedral_num
      self.parm_data['POINTERS'][MPHIA] = dihedral_num
      self.parm_data['POINTERS'][NPTRA] = dihedral_type_num
      self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'] = \
               self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'][:5*dihedral_num]

      # Adjust lengths of the dihedral arrays to make sure they're long enough
      for key in ('DIHEDRAL_FORCE_CONSTANT', 'DIHEDRAL_PERIODICITY',
                  'DIHEDRAL_PHASE', 'SCEE_SCALE_FACTOR', 'SCNB_SCALE_FACTOR'):
         if not key in self.parm_data.keys(): continue
         self.parm_data[key] = [0 for i in range(dihedral_type_num)]
      
      self.dihedral_type_list.write_to_parm()
      
      # Load the pointers now
      self.LoadPointers()
      # Mark atom list as unchanged
      self.atom_list.changed = False
      self.bond_type_list.changed = False
      self.bonds_inc_h.changed = False
      self.bonds_without_h.changed = False
      self.angle_type_list.changed = False
      self.angles_inc_h.changed = False
      self.angles_without_h.changed = False
      self.dihedral_type_list.changed = False
      self.dihedrals_inc_h.changed = False
      self.dihedrals_without_h.changed = False

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
   def _topology_changed(self):
      """ 
      Determines if any of the topological arrays have changed since the
      last upload
      """
      return (self.atom_list.changed or
              self.bond_type_list.changed or self.bonds_inc_h.changed or
              self.bonds_without_h.changed or self.angle_type_list.changed or
              self.angles_inc_h.changed or self.angles_without_h.changed or 
              self.dihedral_type_list.changed or self.dihedrals_inc_h.changed or
              self.dihedrals_without_h.changed)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _fill_res_container(self):
      """ Refills the residue_container array if we changed anything """
      self.residue_container = []
      for i in range(self.parm_data['POINTERS'][NRES]-1):
         curres = self.parm_data['RESIDUE_POINTER'][i] - 1
         nextres = self.parm_data['RESIDUE_POINTER'][i+1] - 1
         for j in range(curres, nextres):
            self.residue_container.append(i+1)
      for i in range(self.parm_data['RESIDUE_POINTER'][
                     self.parm_data['POINTERS'][NRES]-1]-1,
                     self.parm_data['POINTERS'][NATOM]):
         self.residue_container.append(self.parm_data['POINTERS'][NRES])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def delete_mask(self, mask):
      """ Deletes all of the atoms corresponding to an entire mask """
      from chemistry.amber.mask import AmberMask
      # Determine if we were given an AmberMask object or a string. If the
      # latter, turn it into an AmberMask and get the selection
      if isinstance(mask, AmberMask):
         # Make sure the AmberMask's parm is this one!
         if id(self) != id(mask.parm):
            raise AmberParmError('Mask belongs to different prmtop!')
         selection = mask.Selection()
      elif isinstance(mask, str):
         selection = AmberMask(self, mask).Selection()

      # Delete all of the atoms
      for i in reversed(range(len(selection))):
         if selection[i]: del self.atom_list[i]
      
      # Remake the topology file and re-set the molecules if we have periodic
      # boxes (or delete the Molecule info if we removed all solvent)
      self.remake_parm()
      
      # Reconstruct the coordinates and velocities from the remaining atoms
      if hasattr(self, 'coords'):
         self.coords = []
         if self.hasvels: self.vels = []
         for atm in self.atom_list:
            self.coords.extend([atm.xx, atm.xy, atm.xz])
            if self.hasvels: self.vels.extend([atm.vx, atm.vy, atm.vz])
            
      self._load_structure()
      if self.ptr('ifbox'): self.rediscover_molecules()


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def rediscover_molecules(self, solute_ions=True, fix_broken=True):
      """
      This determines the molecularity and sets the ATOMS_PER_MOLECULE and
      SOLVENT_POINTERS sections of the prmtops. Returns the new atom sequence
      in terms of the 'old' atom indexes if re-ordering was necessary to fix the
      tleap bug. Returns None otherwise.
      """
      # Bail out of we are not doing a solvated prmtop
      if not self.ptr('ifbox'): return None

      owner = set_molecules(self)
      ions = ['Br-','Cl-','Cs+','F-','I-','K+','Li+','Mg+','Na+','Rb+','IB',
              'CIO','MG2']
      indices = []
      for res in self.solvent_residues:
         try: indices.append(self.parm_data['RESIDUE_LABEL'].index(res))
         except ValueError: pass
      # Add ions to list of solvent if necessary
      if not solute_ions:
         for ion in ions:
            if ion in self.parm_data['RESIDUE_LABEL']:
               indices.append(self.parm_data['RESIDUE_LABEL'].index(ion))
      # If we have no water, we do not have a molecules section!
      if not indices:
         self.parm_data['POINTERS'][IFBOX] = 0
         self.LoadPointers()
         self.deleteFlag('SOLVENT_POINTERS')
         self.deleteFlag('ATOMS_PER_MOLECULE')
         self.deleteFlag('BOX_DIMENSIONS')
         self.hasbox = False
         try: 
            del self.box
            del self.rst7.box
         except AttributeError:
            # So we don't have box information... doesn't matter :)
            pass
         return None
      # Now remake our SOLVENT_POINTERS and ATOMS_PER_MOLECULE section
      self.parm_data['SOLVENT_POINTERS'] = [min(indices), len(owner), 0]
      first_solvent = self.parm_data['RESIDUE_POINTER'][min(indices)]
      # Find the first solvent molecule
      for i, mol in enumerate(owner):
         if first_solvent-1 == mol[0]:
            self.parm_data['SOLVENT_POINTERS'][2] = i + 1
            break
      else: # this else belongs to 'for', not 'if'
         raise MoleculeError('Could not find first solvent atom!')

      # Now set up ATOMS_PER_MOLECULE and catch any errors
      self.parm_data['ATOMS_PER_MOLECULE'] = [len(mol) for mol in owner]

      # Check that all of our molecules are contiguous, because we have to
      # re-order atoms if they're not
      try:
         for mol in owner:
            for i in range(1, len(mol)):
               if mol[i] != mol[i-1] + 1:
                  raise StopIteration()
      except StopIteration:
         if not fix_broken:
            raise MoleculeError('Molecule atoms are not contiguous!')
         # Non-contiguous molecules detected... time to fix (ugh!)
         warn('Molecule atoms are not contiguous! I am '
              'attempting to fix this, but it may take a while.',
              MoleculeWarning)
         new_atoms = AtomList(self, fill_from=self.atom_list)
         i = 0
         for mol in owner:
            for atm in mol:
               new_atoms[i] = self.atom_list[atm]
               i += 1
         self.atom_list = new_atoms
         return owner

      return None

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def writeOFF(self, off_file='off.lib'):
      """ Writes an OFF file from all of the residues found in a prmtop """
      from chemistry.amber.residue import ToResidue
   
      off_file = open(off_file,'w',0)
   
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
      off_file.write('!!index array str\n')
      for res in residues:
         off_file.write(' "%s"\n' % res.name)
   
      # Now write the OFF strings to the file
      for res in residues:
         off_file.write(res.OFF())

      off_file.close()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def fill_LJ(self):
      """
      Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
      from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
      files, by undoing the canonical combining rules.
      """
      self.LJ_radius = []  # empty LJ_radii so it can be re-filled
      self.LJ_depth = []   # empty LJ_depths so it can be re-filled
      self.LJ_types = {}   # empty LJ_types so it can be re-filled
      one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

      for i in range(self.pointers["NATOM"]): # fill the LJ_types array
         self.LJ_types[self.parm_data["AMBER_ATOM_TYPE"][i]] = \
                                 self.parm_data["ATOM_TYPE_INDEX"][i]
         
      for i in range(self.pointers["NTYPES"]):
         lj_index = self.parm_data["NONBONDED_PARM_INDEX"][
                     self.pointers["NTYPES"] * i + i] - 1
         if self.parm_data["LENNARD_JONES_ACOEF"][lj_index] < 1.0e-10:
            self.LJ_radius.append(0)
            self.LJ_depth.append(0)
         else:
            factor = (2 * self.parm_data["LENNARD_JONES_ACOEF"][lj_index] / 
                        self.parm_data["LENNARD_JONES_BCOEF"][lj_index])
            self.LJ_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_depth.append(self.parm_data["LENNARD_JONES_BCOEF"][lj_index]
                                    / 2 / factor)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def fill_14_LJ(self):
      """
      Fills the LJ_14_radius, LJ_14_depth arrays with data (LJ_types is
      identical) from LENNARD_JONES_14_ACOEF and LENNARD_JONES_14_BCOEF sections
      of the prmtop files, by undoing the canonical combining rules.
      """
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
            factor = (2 * self.parm_data["LENNARD_JONES_14_ACOEF"][lj_index] /
                          self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index])
            self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
            self.LJ_14_depth.append(
                self.parm_data["LENNARD_JONES_14_BCOEF"][lj_index] / 2 / factor)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_LJ(self):
      """
      Takes the values of the LJ_radius and LJ_depth arrays and recalculates the
      LENNARD_JONES_A/BCOEF topology sections from the canonical combining
      rules.
      """
      for i in range(self.pointers["NTYPES"]):
         for j in range(i,self.pointers["NTYPES"]):
            index = self.parm_data['NONBONDED_PARM_INDEX'][
                                          self.ptr('ntypes')*i + j
                                                          ] - 1
            rij = self.combine_rmin(self.LJ_radius[i], self.LJ_radius[j])
            wdij = self.combine_epsilon(self.LJ_depth[i], self.LJ_depth[j])
            self.parm_data["LENNARD_JONES_ACOEF"][index] = wdij * rij ** 12
            self.parm_data["LENNARD_JONES_BCOEF"][index] = 2 * wdij * rij ** 6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def recalculate_14_LJ(self):
      """
      Takes the values of the LJ_radius and LJ_depth arrays and recalculates the
      LENNARD_JONES_A/BCOEF topology sections from the canonical combining rules
      for the 1-4 LJ interactions (CHAMBER only)
      """
      if not self.chamber:
         raise TypeError('recalculate_14_LJ() requires a CHAMBER prmtop!')

      for i in range(self.pointers["NTYPES"]):
         for j in range(i,self.pointers["NTYPES"]):
            index = self.parm_data['NONBONDED_PARM_INDEX'][
                                                   self.ptr('ntypes')*i + j
                                                          ] - 1
            rij = self.combine_rmin(self.LJ_14_radius[i],self.LJ_14_radius[j])
            wdij = self.combine_epsilon(self.LJ_14_depth[i],self.LJ_14_depth[j])
            self.parm_data["LENNARD_JONES_14_ACOEF"][index] = wdij * rij ** 12
            self.parm_data["LENNARD_JONES_14_BCOEF"][index] = 2 * wdij * rij**6

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def combine_rmin(rmin1, rmin2):
      """ Define the combining rule for Rmin """
      return rmin1 + rmin2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def combine_epsilon(eps1, eps2):
      """ Define the combining rule for Epsilon """
      return sqrt(eps1 * eps2)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def LoadRst7(self, filename):
      """ Loads coordinates into the AmberParm class """
      self.rst7 = rst7(filename)
      if not self.rst7.valid:
         raise AmberParmError("Invalid restart file!")
      self.coords = self.rst7.coords
      self.hasvels = self.rst7.hasvels
      self.hasbox = self.rst7.hasbox
      if self.hasbox:
         self.box = self.rst7.box
      if self.hasvels:
         self.vels = self.rst7.vels
      # Load all of the coordinates and velocities into the atoms
      for i in range(len(self.atom_list)):
         self.atom_list[i].xx = self.coords[3*i  ]
         self.atom_list[i].xy = self.coords[3*i+1]
         self.atom_list[i].xz = self.coords[3*i+2]
         if self.hasvels:
            self.atom_list[i].vx = self.vels[3*i  ]
            self.atom_list[i].vy = self.vels[3*i+1]
            self.atom_list[i].vz = self.vels[3*i+2]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def ToMolecule(self):
      """ Translates an amber system into a molecule format """
      from chemistry.molecule import Molecule
      from copy import copy

      # Remake the topology file if it's changed
      if self._topology_changed():
         self.remake_parm()
         if self.ptr('ifbox'): self.rediscover_molecules()
         self._load_structure()

      all_bonds = []        # bond array in Molecule format
      residue_pointers = [] # residue pointers adjusted for indexing from 0
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
         radii.append(
           self.LJ_radius[self.LJ_types[self.parm_data['AMBER_ATOM_TYPE'][i]]-1]
                     )
      try:
         if self.valid and self.rst7.valid:
            return Molecule(atoms=copy(self.parm_data['ATOM_NAME']),
                      atom_types=copy(self.parm_data['AMBER_ATOM_TYPE']),
                      charges=copy(self.parm_data['CHARGE']),
                      residues=copy(self.parm_data['RESIDUE_LABEL']),
                      bonds=all_bonds, residue_pointers=residue_pointers,
                      coords=copy(self.coords), elements=elements,
                      title=title, radii=radii)
      except AttributeError: # use dummy list if no coords are loaded
         if self.valid:
            return Molecule(atoms=copy(self.parm_data['ATOM_NAME']),
                      atom_types=copy(self.parm_data['AMBER_ATOM_TYPE']),
                      charges=copy(self.parm_data['CHARGE']),
                      residues=copy(self.parm_data['RESIDUE_LABEL']), 
                      bonds=all_bonds, residue_pointers=residue_pointers,
                      coords=list(range(self.pointers['NATOM']*3)),
                      elements=elements, title=title, radii=radii)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class rst7(object):
   """ Amber input coordinate (or restart coordinate) file format """
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def __init__(self, filename):
      """ Initialize the inpcrd file """
      self.filename = filename
      self.valid = False

      try:
         self._read()
      except BaseException:
         try:
            self._readnc()
         except ImportError:
            raise ReadError(('Error parsing coordinates from %s. If this is a '
                  'NetCDF restart file, you must install ScientificPython or '
                  'pynetcdf') % self.filename)
         except BaseException, err:
            raise(ReadError('Error parsing coordinates from %s: %s' %
                  (self.filename, err)))

      self.valid = True

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
         self.box[0], self.box[1] = float(self.box[0]), float(self.box[1])
         self.box[2], self.box[3]  = float(self.box[2]), float(self.box[3])
         self.box[4], self.box[5]  = float(self.box[4]), float(self.box[5])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   def _readnc(self):
      """ Reads a NetCDF-format restart file """
      global _HAS_NETCDF
      if not _HAS_NETCDF:
         raise ImportError('Could not import ScientificPython or netCDF4!')
      # Get the helper functions
      global get_int_dimension, get_float_array, get_float

      # Open the restart file and get the title, number of atoms, and time
      restrt = NetCDFFile(self.filename, 'r')
      self.title = str(restrt.title)
      self.natom = get_int_dimension(restrt, 'atom')
      self.coords = get_float_array(restrt, 'coordinates').flatten().tolist()
      self.time = get_float(restrt, 'time')
      self.hasvels = 'velocities' in restrt.variables
      self.hasbox = ('cell_lengths' in restrt.variables
                  and 'cell_angles' in restrt.variables)
      if self.hasvels:
         self.vels = get_float_array(restrt, 'velocities').flatten().tolist()
      else:
         self.vels = []
      if 'time' in restrt.variables:
         self.time = get_float(restrt, 'time')
      if self.hasbox:
         self.box = (get_float_array(restrt, 'cell_lengths').tolist() +
                     get_float_array(restrt, 'cell_angles').tolist())

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Element(mass):
   """ Determines what element the given atom is based on its mass """

   diff = mass
   best_guess = 'EP'

   for element in periodic_table.Element:
      if abs(periodic_table.Mass[element] - mass) < diff:
         best_guess = element
         diff = abs(periodic_table.Mass[element] - mass)

   return best_guess

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# For backwards-compatibility
class amberParm(AmberParm):
   """ 
   This equivalence is made to preserve backwards-compatibility. The standard
   for creating class names is CapitalLettersSeparateWords, with the starting
   letter being capital.
   """
   def __init__(self, prm_name=None, rst7_name=None):
      warn('amberParm is deprecated. Use AmberParm.', DeprecationWarning)
      AmberParm.__init__(self, prm_name, rst7_name)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def set_molecules(parm):
   """ Correctly sets the ATOMS_PER_MOLECULE and SOLVENT_POINTERS sections
       of the topology file.
   """
   from sys import setrecursionlimit, getrecursionlimit
   # Since we use a recursive function here, we make sure that the recursion
   # limit is large enough to handle the maximum possible recursion depth we'll
   # need (NATOM). We don't want to shrink it, though, since we use list
   # comprehensions in list constructors in some places that have an implicit
   # (shallow) recursion, therefore, reducing the recursion limit too much here
   # could raise a recursion depth exceeded exception during a _Type/Atom/XList
   # creation. Therefore, set the recursion limit to the greater of the current
   # limit or the number of atoms
   setrecursionlimit( max(parm.ptr('natom'),getrecursionlimit()) )

   # Unmark all atoms so we can track which molecule each goes into
   parm.atom_list.unmark()

   if not parm.ptr('ifbox'):
      raise MoleculeError('Only periodic prmtops can have Molecule definitions')
   # The molecule "ownership" list
   owner = []
   # The way I do this is via a recursive algorithm, in which
   # the "set_owner" method is called for each bonded partner an atom
   # has, which in turn calls set_owner for each of its partners and 
   # so on until everything has been assigned.
   molecule_number = 1 # which molecule number we are on
   for i in range(parm.ptr('natom')):
      # If this atom has not yet been "owned", make it the next molecule
      # However, we only increment which molecule number we're on if 
      # we actually assigned a new molecule (obviously)
      if not parm.atom_list[i].marked:
         tmp = [i]
         _set_owner(parm, tmp, i, molecule_number)
         # Make sure the atom indexes are sorted
         tmp.sort()
         owner.append(tmp)
         molecule_number += 1
   return owner

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _set_owner(parm, owner_array, atm, mol_id):
   """ Recursively sets ownership of given atom and all bonded partners """
   parm.atom_list[atm].marked = mol_id
   for partner in parm.atom_list[atm].bond_partners:
      if not partner.marked:
         owner_array.append(partner.starting_index)
         _set_owner(parm, owner_array, partner.starting_index, mol_id)
      elif partner.marked != mol_id:
         raise MoleculeError('Atom %d in multiple molecules' % 
                             partner.starting_index)
