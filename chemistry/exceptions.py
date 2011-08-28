""" This module contains all of the exceptions that are used in the chemistry
    package
"""

class ChemError(Exception):
   def __init__(self, msg='Error'):
      self.msg = msg
   def __str__(self):
      return self.msg

class ReadError(ChemError):
   """ Error when files cannot be properly read """
   def __init__(self, msg='Read error.'):
      self.msg = msg

class MoleculeError(ChemError):
   """ Error when molecule is illegally defined """
   def __init__(self, msg='Illegal definition of Molecule.'):
      self.msg = msg

class FileError(ChemError):
   """ Error when something illegal is done with files """
   def __init__(self, msg='Illegal file manipulation'):
      self.msg = msg

class FlagError(ChemError):
   """ Error when a FLAG is manipulated in readparm.amberParm """
   def __init__(self, msg='Bad flag'):
      self.msg = msg

class MaskError(ChemError):
   """ Error when a Mask is poorly formed """
   def __init__(self, msg='Bad mask'):
      self.msg = msg
