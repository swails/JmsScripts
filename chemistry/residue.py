# This module holds the Residue class, which is an organized unit
# of atoms and bonds.

from .atom import Atom
from .bond import Bond

class Residue:
   """ Residue class that contains a collection of atoms and bonds """
   
   #----------------------------------------------------------------
   
   def __init__(self, name, number=0):
      """ Initializes a residue as a collection of atoms and bonds """

      self.Name = name
      self.Number = number
      self.Bonds = []
      self.Atoms = []
      self.StartAtom

   #----------------------------------------------------------------

   def addAtom(self, atom):
      """ Add an atom to this residue """
      if atom.__name__ != "Atom":
         raise ResidueError("addAtom only accepts Atom as a parameter")

      self.Atoms.append(atom)
      if len(self.Atoms) == 1: # if this is the first atom
         self.StartAtom = self.Atoms[0]
      
   #----------------------------------------------------------------

   def addBond(self, atnum1, atnum2, order=1):
      """ Create a bond between the atnum1-th and atnum2-th atoms in the residue """
      # make sure we don't add the same bond twice
      atomnumbers = (self.Atoms[atnum1], self.Atoms[atnum2])
      atomnumbers = sorted(atomnumbers)
      for i in range(len(self.Bonds)):
         if (self.Bonds[i].Atom1, self.Bonds[i].Atom2) == (atomnumbers[0], atomnumbers[1]):
            raise BondError("Bond already exists")

      self.Bonds.append(Bond(self.Atoms[atnum1], self.Atoms[atnum2], order))
      
   #----------------------------------------------------------------

   def AddHydrogens(self):
      """ Add hydrogens to fill up each atom's valence """
