# This module contains the Bond class that represents a chemical
# bond between two Atom objects.

from .atom import Atom

# Define the average bond lengths
BOND_LENGTHS = { ('C' ,'H' , 1) : 1.09, ('H' ,'O' , 1) : 0.96, ('S' ,'H' , 1) : 1.34,
                 ('H' ,'N' , 1) : 1.04 }

class Bond:

   #----------------------------------------------------------------
   
   def __init__(self, atom1, atom2, order=1):
      """ Initialize a bond """
      # make sure that we're bonding 2 Atom types
      if atom1.__name__ != 'Atom' or atom2.__name__ != 'Atom':
         raise BondError('Bond can only exist between two Atoms')

      try:
         atom1.addBond(order)
         atom2.addBond(order)
      except AtomError as err:
         raise BondError(str(err))

      if atom1.Number < atom2.Number:
         self.Atom1 = atom1
         self.Atom2 = atom2
      elif atom1.Number > atom2.Number:
         self.Atom1 = atom2
         self.Atom2 = atom1
      else:
         raise BondError('An atom cannot bond to itself')

   #----------------------------------------------------------------

   def Destroy(self):
      """ Delete this bond. This should be called before it's deleted """
      # remove this bond from the count of each Atom's bonded partners
      self.Atom1.delBond(order)
      self.Atom2.delBond(order)

   #----------------------------------------------------------------
