"""
This module interacts with the Molecule class via the PQR file format.
"""

from ..molecule import Molecule
from os.path import exists as file_exists
from .. import exceptions

# Determines if we can overwrite existing files or not with PutPQR
overwrite = False

# Lists the standard residues for which we have topology data
StandardResidues = ['ALA','ARG','ASP','ASN','CYS','GLN','GLU','GLY','HIS','ILE',
                    'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
                    'A'  ,'C'  ,'G'  ,'U'  ,'DA' ,'DC' ,'DG' ,'DT']
# Account for different names with different protonation states in Amber residue naming
AmberResidues    = {'ALA' : 'ALA',
                    'ARG' : 'ARG',
                    'ASP' : 'ASP',
                    'ASH' : 'ASP',
                    'AS4' : 'ASP',
                    'ASN' : 'ASN',
                    'CYS' : 'CYS',
                    'CYM' : 'CYS',
                    'CYX' : 'CYS',
                    'GLN' : 'GLN',
                    'GLU' : 'GLU',
                    'GLH' : 'GLU',
                    'GLY' : 'GLY',
                    'HIS' : 'HIS',
                    'HIE' : 'HIS',
                    'HID' : 'HIS',
                    'HIP' : 'HIS',
                    'ILE' : 'ILE',
                    'LEU' : 'LEU',
                    'LYN' : 'LYS',
                    'LYS' : 'LYS',
                    'MET' : 'MET',
                    'PHE' : 'PHE',
                    'PRO' : 'PRO',
                    'SER' : 'SER',
                    'THR' : 'THR',
                    'TRP' : 'TRP',
                    'TYR' : 'TYR',
                    'VAL' : 'VAL',
                    'A'   : 'A',
                    'A3'  : 'A',
                    'A5'  : 'A',
                    'AN'  : 'A',
                    'C'   : 'C',
                    'C3'  : 'C',
                    'C5'  : 'C',
                    'CN'  : 'C',
                    'G'   : 'G',
                    'G3'  : 'G',
                    'G5'  : 'G',
                    'GN'  : 'G',
                    'U'   : 'U',
                    'U3'  : 'U',
                    'U5'  : 'U',
                    'UN'  : 'U',
                    'DA'  : 'DA',
                    'DA3' : 'DA',
                    'DA5' : 'DA',
                    'DAN' : 'DA',
                    'DC'  : 'DC',
                    'DC3' : 'DC',
                    'DC5' : 'DC',
                    'DCN' : 'DC',
                    'DG'  : 'DG',
                    'DG3' : 'DG',
                    'DG5' : 'DG',
                    'DGN' : 'DG',
                    'RA'  : 'A',
                    'RA3' : 'A',
                    'RA5' : 'A',
                    'RAN' : 'A',
                    'RC'  : 'C',
                    'RC3' : 'C',
                    'RC5' : 'C',
                    'RCN' : 'C',
                    'RG'  : 'G',
                    'RG3' : 'G',
                    'RG5' : 'G',
                    'RGN' : 'G',
                    'RU'  : 'U',
                    'RU3' : 'U',
                    'RU5' : 'U',
                    'RUN' : 'U' }

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def PutPQR(molecule, filename, standard=False, title='Created by Python chemistry package'):
   """ Writes a PQR based on the PQR version 3 spec """

   # define an array with location of TER cards and call the function to initialize it
   ter_locations = []
   _needTer(molecule, ter_locations)

   # Define a counter for which atom we're printing out
   atom_counter = 1

   # Keep track of whether an atom is ATOM or HETATM so we know what to do for CONECTs
   rectypes = []

   # Keep counters for MASTER record
   numters   = 0
   numconect = 0

   # Make sure we aren't overwriting this file if we haven't said it's OK yet
   if file_exists(filename) and not overwrite:
      raise(exceptions.FileError('Cannot open %s for writing. It already exists.' % filename))

   # Open the file for writing and print out the title
   file = open(filename,'w',0)
   print >> file, "REMARK %s" % title
   print >> file, "REMARK %s" % molecule.title

   # Now print out each line for atom
   for i in range(len(molecule.atoms)):
      record, resname = _resname(molecule.residues[molecule.residue_container[i]], standard)
      rectypes.append(record)
      line = "%-6s %4d %4s %3s  %4i    %8.3f %7.3f %7.3f %7.4f %7.4f      %2s" % \
             (record, atom_counter % 10000, _format(molecule.atoms[i]), resname,
              molecule.residue_container[i]+1, molecule.coords[3*i], molecule.coords[3*i+1], 
              molecule.coords[3*i+2], molecule.charges[i], molecule.radii[i], molecule.elements[i])
      print >> file, line
      atom_counter += 1
      try:
         if molecule.residue_container[i] != molecule.residue_container[i+1] and \
               molecule.residue_container[i] in ter_locations:
            print >> file, 'TER'
            numters += 1
      except IndexError:
         pass

   atom_counter -= 1 # decrement atom counter, since it was incremented after the last atom
   # Now print out the CONECT cards if they are needed
   for i in range(len(rectypes)):
      if rectypes[i] == 'HETATM' and molecule.residues[molecule.residue_container[i]] != 'WAT':
         line = 'CONECT%5d' % (i + 1)
         for j in range(len(molecule.bonds[i])):
            line += '%5d' % (molecule.bonds[i][j] + 1)
         print >> file, line
         numconect += 1

   # Now print out the MASTER record
   print >> file, 'MASTER        1    0    0    0    0    0    0    0 %4d %4d %4d    0' % \
            (atom_counter, numters, numconect)

   # End the PQR
   print >> file, 'END'

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _format(name):
   """ Formats an atom name for use in PQR files """
   
   if len(name) == 4: return '%4s' % name
   else: return ' %-3s' % name.strip()

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _needTer(molecule, ter_locations):
   """ Determines where we need to put TER cards in the PQR file """
   # Only check to see if this residue is bonded to a residue that
   # comes later in the sequence. If so, we don't put a TER card
   for i in range(len(molecule.residues)-1):
      bonded  = False
      curres  = molecule.residue_pointers[i]
      nextres = molecule.residue_pointers[i+1]
      for j in range(curres, nextres):
         for k in molecule.bonds[j]:
            if k >= nextres:
               bonded = True
               break
         if bonded: break
      if not bonded: ter_locations.append(i)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _resname(residue, standard):
   """ Determines how a residue should be printed and if the record should be
       ATOM or HETATM """
   res = residue
   if residue in AmberResidues.keys():
      rec = 'ATOM'
      if standard:
         res = AmberResidues[residue]
   else:
      rec = 'HETATM'

   return rec, res

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
