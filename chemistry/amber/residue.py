"""
Holds an Amber residue class for easy writing of mol2 and
OFF files
"""

from ..molecule import Residue

class AmberResidue(Residue):
   """ A residue that has some amber-specific attributes; child of molcule.Residue """

   def __init__(self):
      """ Initialization of AmberResidue object """
      # Adjust all indexes so that it starts from 1 rather than 0
      self.is_nterm = self.head == 0
      self.is_cterm = self.tail == 0
