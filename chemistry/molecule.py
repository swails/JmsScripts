""" 
Molecule class for manipulating biomolecular structures and storing
molecular data for use in various transformations.
"""

from . import periodic_table
from . import exceptions

class Molecule:
   """ 
   Molecule class that contains biomolecular data:
      o  atoms       : sequential list of atoms in an array
      o  atom_types  : atom types associated with each atom
      o  elements    : element of each atom
      o  charges     : partial charges on each atom
      o  residues    : array of residue names
      o  bonds       : array in which each element lists the atom indices
                       of the other atoms that atom is bonded to
      o  residue_pointers : index of first atom of each residue
      o  residue_container: which residue each atom belongs to
      o  coords      : cartesian coordinates (x1,y1,z1,x2,y2,z2,...) of
                       each atom in the molecule
   """

   def __init__(self, atoms=[], atom_types=[], charges=[], residues=[], bonds=[], 
                residue_pointers=[], coords=[], elements=[], title='',radii=[]):
      """ Initializing and checking the molecular data """
      self.atoms = atoms
      self.residues = residues
      self.residue_pointers = residue_pointers
      self.coords = coords
      self.bonds = bonds
      self.atom_types = atom_types
      self.charges = charges
      self.elements = elements
      self.title = title
      self.radii = radii
      
      self._check()
      if self.valid:
         self._fillcontainer()

   def DeleteAtom(self, atomno):
      """ Deletes an atom from the molecule, removing all bonds it's involved with
          and translating all bonds by one """
      # Remove atomno atom from all of the arrays
      self.atoms.pop(atomno)
      self.coords.pop(atomno*3 + 2)
      self.coords.pop(atomno*3 + 1)
      self.coords.pop(atomno*3    )
      self.atom_types.pop(atomno)
      self.charges.pop(atomno)
      self.elements.pop(atomno)
      if len(self.radii) > 0:
         self.radii.pop(atomno)

      # Remove atomno from all of the bonds
      self.bonds.pop(atomno)
      for i in range(len(self.bonds)):
         for j in range(len(self.bonds[i])):
            k = len(self.bonds[i]) - 1 - j # go in reverse order
            if self.bonds[i][k] == atomno:
               self.bonds[i].pop(k)
            elif self.bonds[i][k] > atomno:
               self.bonds[i][k] -= 1
      
      # Adjust residue_pointers
      for i in range(self.residue_container[atomno] + 1, len(self.residue_pointers)):
         self.residue_pointers[i] -= 1

      # Get rid of now-vacant residues if we removed the last atom from that residue
      if self.residue_pointers[self.residue_container[atomno]+1] - \
         self.residue_pointers[self.residue_container[atomno]+1] == 0:
         self.residue_pointers.pop(self.residue_container[atomno])
         self.residues.pop(self.residue_container[atomno])

      # Now get rid of atomno from residue_container
      self.residue_container.pop(atomno)


   def _check(self):
      """ Checks for consistency in the molecule """
      self.valid = False
      if len(self.atoms) != len(self.atom_types):
         raise(exceptions.MoleculeError('len(atoms) != len(atom_types)'))

      if len(self.charges) == 0:
         for i in range(len(self.atoms)):
            self.charges.append(0.0)

      if len(self.charges) != len(self.atoms):
         raise(exceptions.MoleculeError('len(atoms) != len(charges)'))

      if len(self.coords) != 3 * len(self.atoms):
         raise(exceptions.MoleculeError('len(atoms) != len(coords) * 3'))

      if len(self.residue_pointers) != len(self.residues):
         raise(exceptions.MoleculeError('len(residue_pointers) != len(residues)'))

      if len(self.elements) != len(self.atoms):
         raise(exceptions.MoleculeError('len(elements) != len(atoms)'))

      if len(self.radii) != 0 and len(self.radii) != len(self.atoms):
         raise(exceptions.MoleculeError('len(radii) != len(atoms)'))

      self.valid = True

   def _fillcontainer(self):
      """ Fills residue_container so we know what residue each atom is in """
      self.residue_container = []
      # Fill residues 1 - nres-1
      for i in range(len(self.residues)-1):
         curres = self.residue_pointers[i]
         nextres = self.residue_pointers[i+1]
         for j in range(curres, nextres):
            self.residue_container.append(i)
      # Fill the last residue
      for i in range(self.residue_pointers[len(self.residue_pointers)-1], len(self.atoms)):
         self.residue_container.append(len(self.residues)-1)

   def ToResidue(self, resnum):
      """ Returns a residue object from the residue passed """
      begin = self.residue_pointers[resnum]
      end = self.residue_pointers[resnum+1] # 1 after the end, as it needs to be for array slices

      # Gather the bonding information for the given residue. Create a new bond array that is specific
      # to this residue (and whose partners are indexed as though this residue starts from atom 0). Also
      # find the "connect" atoms that bond to adjacent residues.
      bonds = []
      connects = []
      head = tail = -1
      for i in range(begin, end):
         bonds.append([])
         for j in self.bonds[i]:
            if j < begin and head != -1:
               connects.append(i - begin)
               connects.append(-1)
               head = i - begin
            elif j < begin:
               connects.append(i - begin)
            elif j >= end:
               if self.residue_container[j] == resnum + 1:
                  # we keep track of index in connects and the atom number so we can
                  # remove it from here and replace it in connects[1] at the end.
                  # we do this so that we only count the last "tail" residue as connect1
                  tail = len(connects)
                  connects.append(i - begin)
               else:
                  connects.append(i - begin)
            else: # otherwise it's an intra-residue bond. adjust j for new numbering and add it to bonds
               bonds[i].append(j - begin)

      # Now find the tail (if there is one) and put it in the 2nd location
      if tail != -1:
         connects[1] = connects.pop(tail)
         tail = connects[1]

      return Residue(self.atoms[begin:end], self.atom_types[begin:end], self.charges[begin:end],
                     bonds, self.coords[3*begin:3*end], head, tail)

class Residue:
   """ This is a defined residue in a biomolecular system """
   def __init__(self, atoms=[], atom_types=[], charges=[], bonds=[], coords=[], head=-1, tail=-1):
      """ initializes the residue """
      self.atoms = atoms
      self.atom_types = atom_types
      self.charges = charges
      self.bonds = bonds
      self.coords = coords
      self.head = head
      self.tail = tail
      self.natom = len(self.atoms)
