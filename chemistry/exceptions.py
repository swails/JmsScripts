""" This module contains all of the exceptions that are used in the chemistry
    package
"""

class ReadError(Exception):
   """ Error when files cannot be properly read """
   def __init__(self, msg='Read error.'):
      """ initializes ReadError """
      self.msg = msg
   def __str__(self):
      """ Returns error message """
      return self.msg

class MoleculeError(Exception):
   """ Error when molecule is illegally defined """
   def __init__(self, msg='Illegal definition of Molecule.'):
      """ Initializes MoleculeError """
      self.msg = msg
   def __str__(self):
      """ Returns error message """
      return self.msg

class FileError(Exception):
   """ Error when something illegal is done with files """
   def __init__(self, msg='Illegal file manipulation'):
      """ Initializes FileError """
      self.msg = msg
   def __str__(self):
      """ Returns error message """
      return self.msg
