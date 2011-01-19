"""
This module interacts with the Molecule class via several different file formats
(i.e. PDB, mol2, PQR, etc.). 
"""

from ..molecule import Molecule
from os.path import exists as file_exists
from .. import exceptions


overwrite = False

def PutPDB(molecule, filename):
   """ Writes a PDB based on the PDB version 3 spec """
   ter_locations = []
   _needTer(molecule, ter_locations)

   atom_counter = 1

   if file_exists(filename) and not overwrite:
      raise(exceptions.FileError('Cannot open %s for writing. It already exists.' % filename))

   file = open(filename,'w',0)
   print >> file, "REMARK Created by Python package chemistry"

   for i in range(len(molecule.atoms)):
      line = "%-6s%4d %4s %3s A%4i    %8.3f%8.3f%8.3f                      %s" % \
             ('ATOM', atom_counter, _format(molecule.atoms[i]), molecule.residues[molecule.residue_container[i]],
              molecule.residue_container[i]+1, molecule.coords[3*i], molecule.coords[3*i+1], molecule.coords[3*i+2], 
              molecule.atom_types[i])
      print >> file, line
      try:
         atom_counter += 1
         if molecule.residue_container[i] != molecule.residue_container[i+1] and \
               molecule.residue_container[i] in ter_locations:
            print >> file, 'TER'
      except IndexError:
         pass

def _format(name):
   """ Formats an atom name for use in PDB files """
   try:
      int(name[len(name)-1])
      return '%4s' % name
   except ValueError:
      return '%3s ' % name.strip()

def _needTer(molecule, ter_locations):
   """ Determines where we need to put TER cards in the PDB file """
   for i in range(len(molecule.residues)-1):
      determined = False
      curres  = molecule.residue_pointers[i]
      nextres = molecule.residue_pointers[i+1]
      for j in range(curres, nextres):
         for k in molecule.bonds[j]:
            if k >= nextres:
               determined = True
               break
         if determined: break
      if not determined: ter_locations.append(i)
