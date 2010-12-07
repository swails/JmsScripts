"""
Define the exceptions used for Molecules
"""

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomError(Exception):
   """ Errors associated with the Atom class """
   
   #----------------------------------------------------------------

   def __init__(self, msg='Error defining Atom'):
      """ Initialization of AtomError """
      self.message = msg
   
   #----------------------------------------------------------------

   def __str__(self):
      """ String representation of AtomError """
      return self.message

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ElementError(Exception):
   """ Unrecognized elements """
   
   #----------------------------------------------------------------

   def __init__(self, msg='Unrecognized element'):
      """ Initialization of ElementError """
      self.message = msg

   #----------------------------------------------------------------

   def __str__(self):
      """ String representation of ElementError """
      return self.message

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondError(Exception):
   """ Invalid bond """

   #----------------------------------------------------------------

   def __init__(self, msg='Invalid Bond'):
      """ Initialization of BondError """
      self.message = msg

   #----------------------------------------------------------------

   def __str__(self):
      """ String representation of ElementError """
      return self.message

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueError(Exception):
   """ Invalid residue """

   #----------------------------------------------------------------

   def __init__(self, msg='Invalid Residue'):
      """ Initialization of ResidueError """
      self.message = msg

   #----------------------------------------------------------------

   def __str__(self):
      """ String Representation of ResidueError """
      return self.message

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
