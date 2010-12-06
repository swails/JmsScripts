# molecule.py:
# This is a molecule class that contains atom names, residue names,
# atom numbers, etc. It contains a generic molecule as well as an
# AMBER molecule class.

# import useful intrinsic python modules
from copy import copy
import sys, os

# import custom modules
from periodic_table import PeriodicTable
from readparm import amberParm

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MoleculeError(Exception):
   """ General class of error for molecules """
   def __init__(self,message=''):
      self.msg = message

   def __str__(self):
      return self.msg

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MissingBox(Exception):
   """ Exception thrown if the AmberMolecule object is missing expected box information """
   pass

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Molecule:
   """ Molecule class that is useful for computational chemistry. Contains
       atoms and residues, along with information about each atom. The idea
       is to be able to print them out in different formats. """

   # Class variables:
   #
   #  BOND_CUTOFF : cut-off value in Angstroms for what is considered a bond
   #                if no bonding information is found/given
   #  BOND_DICT   : dictionary containing element-specific Bond cutoff lengths. The 2 elements
   #                should be alphabetically-ordered in an array as the key. If no key is found
   #                for the atom pair, then the BOND_CUTOFF will be used. The default will have only
   #                those bonds common to Simulations
   #  MAX_BONDS   : Dictionary that matches each element to the maximum number of atoms it can be 
   #                bonded to. This is the typical value in biological systems. It is ignored if
   #                bonding information is provided.
   #  ignore_max_bond : logical for whether or not the maximum bond should be a fatal limit or
   #                    quietly ignored

   BOND_CUTOFF = 1.6  # bond cutoff is 1.6 angstroms
   BOND_DICT   = { ('C' ,'H' ) : 1.2, ('C' ,'Cl') : 1.9, ('Br','C' ) : 2.0, ('C' ,'C' ) : 1.6, 
                   ('C' ,'Fe') : 2.0, ('C' ,'N' ) : 1.6, ('C' ,'O' ) : 1.6, ('C' ,'S' ) : 1.9,
                   ('O' ,'P' ) : 1.7, ('O' ,'S' ) : 1.6, ('H' ,'O' ) : 1.1, ('H', 'S' ) : 1.4,
                   ('H' ,'N' ) : 1.1, ('N' ,'O' ) : 1.4, ('H' ,'P' ) : 1.5 }
   MAX_BONDS   = { 'H'  : 1, 'Be' : 2, 'B'  : 3, 'C'  : 4, 'N'  : 3, 'O'  : 2, 'F'  : 1, 'Ne' : 0,
                   'Na' : 1, 'Mg' : 2, 'Al' : 3, 'Si' : 4, 'P'  : 5, 'S'  : 6, 'Cl' : 1, 'Ar' : 0 }
   ignore_max_bond = False

   # Instance variables:
   #
   #  SystemName  : Name of the system
   #  AtomNames   : array with each atom name
   #  Elements    : array with each atom's element
   #  Coordinates : x,y,z coordinates of each atom
   #  Residues    : Residue names in the system
   #  ResPointers : Index array for which atom each res begins with
   #  Bonds       : Bond; connectivity
   #  Mass        : array with masses (g/mol) of each atom
   #  Charge      : array with partial atomic charges (fraction e)
   #  natom       : number of atoms in the system
   #  valid       : Boolean for if the molecule is valid or not

#-------------------------------------------

   def __init__(self, name, atoms, elements, coords, residues, res_ptr, bonds={}, mass=[], charge=[]):
      """ Instantiate molecule information:
            o  SystemName .. name of the system of this instance
            o  AtomNames .................... names of each atom
            o  Elements ............. which element each atom is
            o  Coordinates ........... x,y,z coords of each atom
            o  Residues .......................... residue names
            o  ResPointers ..... index starting atom of each res
            o  Bonds ..........................  Bond dictionary
            o  Mass .......................... mass of each atom
            o  Charge ...................... charge of each atom
      """
      self.SystemName = name
      self.AtomNames = copy(atoms)
      self.Elements = copy(elements)
      self.Coordinates = copy(coords)
      self.Residues = copy(residues)
      self.ResPointers = copy(res_ptr)
      self.Bonds = copy(bonds)
      self.Mass = copy(mass)
      self.Charge = copy(charge)
      self.natom = len(atoms)


      # fill the mass array if it's not already filled
      if len(self.Mass) == 0:
         self._fillMasses()

      # fill the charge array with zeroes if it's not already filled
      if len(self.Charge) == 0:
         for i in range(len(self.natom)):
            self.Charge.append(0.0)
         
      # validate the molecule: _validate() changes self.valid to True if it's valid
      self.valid = False
      self._validate()

      # Determine the connectivity if self.Bonds is an empty dictionary
      if len(self.Bonds.keys()) == 0 and self.valid:
         self._fillBond()

#-------------------------------------------

   def _validate(self):
      """ Checks to make sure all arrays are of the proper length """
      if len(self.AtomNames) != self.natom:
         raise MoleculeError('AtomNames invalid length (%d %d)' % (len(self.AtomNames), self.natom))

      if len(self.Coordinates) != self.natom * 3:
         raise MoleculeError('Coordinates invalid length (%d %d)' % (len(self.Coordinates), self.natom*3))

      if len(self.Mass) != self.natom:
         raise MoleculeError('Mass invalid length (%d %d)' % (len(self.Mass), self.natom))

      if len(self.Charge) != self.natom:
         raise MoleculeError('Charge invalid length (%d %d)' % (len(self.Charge), self.natom))

      self.valid = True

#-------------------------------------------

   def _fillMasses(self):
      """ Fills the Mass array with atomic masses from the periodic table """
      self.Mass = []
      for i in range(len(natom)):
         self.Mass.append(PeriodicTable.Mass[self.Elements[i]])

#-------------------------------------------

   def _fillBond(self):
      """ Determines the connectivity based on distance cut-off criteria. """
      self.Bonds = {}

      # this is an invalid molecule until proven otherwise
      self.valid = False

      # Loop from atoms i = 1 -- (natom-1) and compare distance with atoms from 
      # j = i + 1 -- natom, so we don't double-count.
      for i in range(1, self.natom):
         self.Bonds[i] = []
         i3 = i * 3 - 3
         for j in range(i+1, self.natom+1):
            try:
               cutdist = Molecule.BOND_DICT[tuple(sorted(self.Elements[i], self.Elements[j]))]
            except KeyError:
               cutdist = Molecule.BOND_CUTOFF
            cutdist *= cutdist
            j3 = j * 3 - 3
            xij = self.Coordinates(i3  ) - self.Coordinates(j3  )
            yij = self.Coordinates(i3+1) - self.Coordinates(j3+1)
            zij = self.Coordinates(i3+2) - self.Coordinates(j3+2)
            dij = xij * xij + yij * yij + zij * zij
            if dij <= cutdist:
               self.Bonds[i].append(j)

      # Now, if we are not allowing over-loading of bonds, check to make sure we haven't gone
      # overboard. Check one atom at a time (have to check all of the atoms before it, since
      # if atom 2 is bonded to atom 1, it only appears in the list for atom 1).
      if not Molecule.ignore_max_bond:
         for i in range(1, self.natom):
            numbonds = len(self.Bonds[i])
            for j in range(1,i):
               if i in self.Bonds[j]:
                  numbonds += 1
            try:
               if numbonds > Molecule.MAX_BONDS[self.Elements[i]]:
                  raise MoleculeError('Too many bonds to %s (%d)' % (self.Elements[i], numbonds))
            except KeyError:
               pass
      # we've made it through, so it must be OK.
      self.valid = True

#-------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmberMolecule:
   """ A specific type of molecule formed from an Amber topology file and restart/coordinate file """

   # Instance variables:
   #
   #  Topology    :  amberParm topology object
   #  Coordinates :  x,y,z coordinates of each atom
   #  Velocities  :  x,y,z velocities of each atom
   #  BoxInfo     :  Box information (sizes, angles) if present
   #
   # External functions:
   #
   #  ToMolecule  :  returns a molecule object of the prmtop in the system

#-------------------------------------------

   def __init__(self, parm7, rst7):
      """ Instantiate the AmberMolecule object with the topology and coordinate files
            o  Topology ........ Amber Topology file; amberParm object
            o  Coordinates ......... x,y,z coords taken from rst7 file
            o  Velocities ...... x,y,z velocities taken from rst7 file
            o  BoxInfo ........ periodic box info taken from rst7 file
      """

      self.Topology = amberParm(parm7)
      self.Coordinates = []
      self.Velocities = []
      self.BoxInfo = []
      
      self.valid = self.Topology.valid

      # Load the coordinates from the rst7 format file
      self._loadCoords(rst7)

#-------------------------------------------

   def _loadCoords(self, rst7):
      """ Loads coordinates from Amber rst7 file """

      # In this function, we will look for coordinates and velocities, and load that data
      # into the respective instance arrays. A couple notes:
      #  o  Do NOT parse rst7 file if self.valid is False, since this means we have a bad prmtop
      #  o  Use the IFBOX pointer from the prmtop to determine if we will look for box info.
      #
      # Variables:
      #  hasbox            : logical variable for if we should look for a box or not
      #  natom             : number of atoms in the system
      #  natom2            : natom * 2
      #  coords_complete   : how many atoms' coordinates have been added to arrays
      #  line              : line holder used for parsing files
      #  crdfile           : file object for rst7 coordinate file
      
      if not self.valid: return

      # if we make it through this subroutine, turn this to True
      self.valid = False

      hasbox = self.Topology.ptr('ifbox') > 0
      natom = self.Topology.ptr('natom')
      natom2 = natom * 2

      # now open the file to read the coordinates
      crdfile = open(rst7, 'r')
      
      # eat title line
      line = crdfile.readline()

      # determine number of atoms and check for consistency
      line = crdfile.readline()
      if int(line.strip()) != natom:
         raise MoleculeError('NATOM mismatch between %s and %s' % (self.Topology, rst7))

      line = crdfile.readline()
      coords_complete = 0
      done = False

      while line != '':
         # skip over all blank lines
         if len(line.strip()) == 0: 
            line = crdfile.readline()
            continue

         # add the 6 coordinates on this line to the self.Coordinates array
         if coords_complete < natom:
            self.Coordinates.append(float(line[ 0:12]))
            self.Coordinates.append(float(line[12:24]))
            self.Coordinates.append(float(line[24:36]))
            coords_complete += 1
            if coords_complete == natom:
               line = crdfile.readline()
               continue
            self.Coordinates.append(float(line[36:48]))
            self.Coordinates.append(float(line[48:60]))
            self.Coordinates.append(float(line[60:72]))
            coords_complete += 1
            line = crdfile.readline()
            continue

         # add the 6 velocities on this line to the self.Velocities array
         if coords_complete < natom2:
            self.Velocities.append(float(line[ 0:12]))
            self.Velocities.append(float(line[12:24]))
            self.Velocities.append(float(line[24:36]))
            coords_complete += 1
            if coords_complete == natom2:
               line = crdfile.readline()
               continue
            self.Velocities.append(float(line[36:48]))
            self.Velocities.append(float(line[48:60]))
            self.Velocities.append(float(line[60:72]))
            coords_complete += 1
            line = crdfile.readline()
            continue

         # if we've made it this far, we're getting box info
         if len(self.BoxInfo) != 0 or not hasbox:
            raise MoleculeError('Too many lines found in %s' % rst7)

         self.BoxInfo = [float(line[ 0:12]), float(line[12:24]), float(line[24:36]),
                         float(line[36:48]), float(line[48:60]), float(line[60:72])]

         line = crdfile.readline() # read the next line
      
      crdfile.close()

      # we're out of the while loop now -- check to see if we caught 6 velocities that really
      # should be BoxInfo. We will try to fill BoxInfo only if hasbox is True. If it is false,
      # and we have 6 velocities and NOT natom=2, then we raise a MoleculeError. Require a box
      # if the prmtop says so.

      if not hasbox and (len(self.Velocities) == 6 and natom != 2):
         raise MoleculeError('Too many lines found in %s' % rst7)

      if hasbox and len(self.Velocities) == 6 and len(self.BoxInfo) == 0:
         self.BoxInfo = copy(self.Velocities)
         self.Velocities = []
      elif len(self.Velocities) == 6 and natom != 2 and len(self.BoxInfo) != 0:
         raise MoleculeError('Too many lines found in %s' % rst7)

      if len(self.Coordinates) != natom * 3:
         raise MoleculeError('Not enough coordinates found')

      if len(self.Velocities) != 0 and len(self.Velocities) != natom * 3:
         raise MoleculeError('Not enough velocities found')

      if hasbox and len(self.BoxInfo) == 0:
         raise MissingBox()

      self.valid = True

#-------------------------------------------

   def ToMolecule(self):
      """ Returns a Molecule object """
      return Molecule(self.Topology.prm_name, self.Topology.parm_data['ATOM_NAME'], self._elements(),
                      self.Coordinates, self.Topology.parm_data['RESIDUE_LABEL'], self.Topology.parm_data['RESIDUE_POINTER'],
                      self._bondArray(), self.Topology.parm_data['MASS'], self.Topology.parm_data['CHARGE'])

#-------------------------------------------

   def _bondArray(self):
      """ Returns the Bond dictionary, in which the atom number is an index for the
          sorted array of atoms that it's bonded to. This array is minimalistic. i.e.,
          it only lists bonded pairs once in the array of the smaller atom number. So
          if 8 is bonded to 10, bonds[8] = [...,10,...], but bonds[10] will NOT
          contain 8 """

      bonds = {}
      natom = self.Topology.ptr('natom')

      for i in range(1,natom):
         bonds[i] = []

      # first take care of bonds including hydrogens
      for i in range(self.Topology.ptr('nbonh')):
         i3 = i * 3
         atom1 = self.Topology.parm_data['BONDS_INC_HYDROGEN'][i3] / 3 + 1
         atom2 = self.Topology.parm_data['BONDS_INC_HYDROGEN'][i3 + 1] / 3 + 1
         bonds[min(atom1,atom2)].append(max(atom1,atom2))

      # next take care of bonds not including hydrogens
      for i in range(self.Topology.ptr('nbona')):
         i3 = i * 3
         atom1 = self.Topology.parm_data['BONDS_WITHOUT_HYDROGEN'][i3] / 3 + 1
         atom2 = self.Topology.parm_data['BONDS_WITHOUT_HYDROGEN'][i3 + 1] / 3 + 1
         bonds[min(atom1,atom2)].append(max(atom1,atom2))

      # now sort all of the arrays for all of the atoms
      for i in bonds.keys():
         bonds[i] = sorted(bonds[i])

      return bonds

#-------------------------------------------

   def _elements(self):
      """ Determines which element each atom is by comparing molecular masses, since this is really
          the only way of telling for *sure*, since weird atom names can be given, i.e. CA can be
          an alpha-Carbon or a calcium. The mass makes it easy to distinguish. """

      # This is going to be done in a very inefficient manner, but it's the best way I can think of
      # to bypass most kinds of errors, and it's quite easy to implement. I'm going to take the difference
      # between the atomic mass of atom "i" with every atom on the periodic table, and whichever difference
      # has the smallest unsigned value will be identified as that element. However, as inefficient as it
      # is, it still scales as O(N), just with a prefactor of 118 :).

      elements = []  # array of elements to be returned
      min_diff = 999 # minimum difference tracked along
      diff_el = '  ' # element corresponding 2 min_diff

      for i in range(self.Topology.ptr('natom')):
         for j in PeriodicTable.Mass.keys():
            diff = abs(self.Topology.parm_data['MASS'][i] - PeriodicTable.Mass[j])
            if diff < min_diff:
               min_diff = diff
               diff_el = j
         elements.append(diff_el)

      return elements

#-------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
