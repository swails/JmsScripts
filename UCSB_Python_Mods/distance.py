#!/usr/bin/env python


Usage = """ Calculates the distance between two residues at the CA

Usage	: distance.py PDB RES1 RES2

PDB	: PDB file
RES1	: residue number of the first residue
RES2	: residue number of the second residue
"""

import os, sys
if len(sys.argv) !=4:
  print Usage
  sys.exit()
import protein, geometry, coords


#get arguments
Pdb = sys.argv[1]
ResNum1 = int(sys.argv[2]) - 1
ResNum2 = int(sys.argv[3]) - 1

#calculates distance between residues
p = protein.ProteinClass(Pdb=Pdb)
CA1 = p.AtomNum(ResNum1, "CA")
CA2 = p.AtomNum(ResNum2, "CA")
Distance = geometry.Length(p.Pos[CA1] - p.Pos[CA2])
print round(Distance,3)



