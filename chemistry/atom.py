# This module defines the atom class

from . import periodic_table
from ._exceptions import AtomError, ElementError

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _getElement(mass):
   """ Finds the element based on the mass -- looks for the smallest unsigned
       deviation from every element in the periodic_table """
   deviation = 9999
   for element in periodic_table.Element:
      if abs(mass - periodic_table.Mass[element]) > deviation:
         deviation = abs(mass - periodic_table.Mass[element])
         best_guess = element
   return best_guess

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Atom:
   """ Atom class; These are the building blocks of all molecules """

   __name__    = "Atom" # define it as an atom
   strict_bond = True   # enforce strict limits on bonds for small atoms

   #----------------------------------------------------------------

   def __init__(self, number=0, element='UNSPECIFIED', mass=-1, charge=0.0, maxbond=-1, name='UNSPECIFIED'):
      """ Initialize an atom """

      # Attributes: Mass, Element, Charge, Number, BondedTo, MaxBond, Name, Location

      # either the element or the mass must be specified
      if element == 'UNSPECIFIED' and mass == -1:
         raise AtomError('Atom element and mass not specified')

      # Allow the element to be given as an integer (atomic number)
      if element != 'UNSPECIFIED':
         try:
            element = periodic_table.Element[int(element)]
         except ValueError:
            pass

         if not element in periodic_table.Element:
            raise ElementError('Element %s not recognized' % element)

      self.Mass = float(mass)
      self.Element = element
      self.Number = int(number)
      self.Charge = float(charge)
      self.BondedTo = 0
      self.MaxBond = maxbond
      self.Name = name
      self.Location = [0,0,0]

      # If the mass is not specified
      if self.Mass == -1:
         self.Mass = periodic_table.Mass[element]

      # If the element is not specified, find the element that is closest in mass
      if self.Element == 'UNSPECIFIED':
         self.Element = _getElement(mass)

      # Determine how many bonds we can typically have for this atom type
      if self.MaxBond == -1:
         self.MaxBond = periodic_table.MaxBond[self.Element]

      # Default name is just the element name
      if self.Name == 'UNSPECIFIED':
         self.Name = self.Element

   #----------------------------------------------------------------

   def DeletedAtom(self, number):
      """ Shifts the atom number if an "earlier" atom was deleted """
      if number < self.Number:
         self.Number -= 1

   #----------------------------------------------------------------

   def AddedAtom(self, number):
      """ Shifts the atom number if an "earlier" atom was added """
      if number < self.Number:
         self.Number += 1

   #----------------------------------------------------------------

   def addBond(self, order):
      """ Called when this atom becomes part of a bond """
      if self.BondedTo == self.MaxBond:
         raise AtomError('Too many bonds (%d) to %s' % (self.BondedTo + 1, self.Element))
      self.BondedTo += order

   #----------------------------------------------------------------

   def delBond(self, order):
      """ Called when a bond containing this atom is broken/destroyed """
      self.BondedTo -= order

   #----------------------------------------------------------------

   def place(self, x, y, z):
      """ Give a location to this atom """
      self.Location[0] = float(x)
      self.Location[1] = float(y)
      self.Location[2] = float(z)

   #----------------------------------------------------------------
