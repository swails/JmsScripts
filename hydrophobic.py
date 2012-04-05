#!/usr/bin/env python



Usage = """Extracts hydrophobic resiudes that are within the cutoff distance

Usage           : hydrophobic.py InPdb Cutoff 

IndPdb		: input pdb
Cutoff		: cutoff distance to calculate hydrophobic interaction		
"""

import sys, protein, geometry

if len(sys.argv)==1:
  print Usage
  sys.exit()

#get arguments
Pdb = sys.argv[1]
Cutoff = int(sys.argv[2])

p = protein.ProteinClass(Pdb=Pdb)

#check residues for hydrophobic residues 
ResInd=p.ResInd(Hydrophobic=True)
Pos= p.ResPos()
N=len(p)
#print Pos

#Determine which hydrophobic residues are within the cutoff
for i in range(N): 
  if not i in ResInd: continue  
  for j in range(i+1,N):
      if not j in ResInd: continue    
      if geometry.Length(Pos[i] - Pos[j])  <= Cutoff:
        print "%6d %6d %6d" % (i+1, j+1, geometry.Length(Pos[i] - Pos[j]))
  
   
