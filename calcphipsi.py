#!/usr/bin/env python

#LAST MODIFIED: 09-23-09

Usage = """Calculates phi/psi deviations between pdb files or between
a trajectory and a pdb file.

Usage     : calcphipsi.py [OPTIONS] REFPDBFILE PDBFILE1 PDBFILE2 ...
        OR: calcphipsi.py [OPTIONS] REFPDBFILE TRJFILE1 TRJFILE2 ...

OPTIONS   : "--nskip=X" number of configs in trajectory to skip (default is 0)
            "--nread=X" number of configs in trajectory to read; -1 is all (default -1)
            "--nstride=X" read configs every nstride frames (default is 1)

NOTES     : assumes TRJFILES are named "X.mdtrj.crd[.gz]" and looks for
            prmtop files named "X.prmtop.parm7"
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()

from numpy import *  
import copy, os, coords, random
import geometry, sequence, protein, scripttools


#======== COMMAND-LINE RUNNING ========

def Disp(pRef, AvgDiff, VarDiff, NAdd):
  AvgDiff = AvgDiff / (NAdd + 1.e-30)
  VarDiff = VarDiff / (NAdd + 1.e-30)
  VarDiff = VarDiff - AvgDiff**2
  StdDiff = sqrt(VarDiff)
  
  StartRes = 0
  while StartRes < NRes and NAdd[StartRes] == 0:
    StartRes += 1

  StopRes = NRes
  while StopRes > 0 and NAdd[StopRes-1] == 0:
    StopRes -= 1

  if StopRes <= StartRes:
    print "No sequence overlap between reference and other structures."
    exit()

  print "%-8s %-8s %-10s %-10s %-10s" % ("ResNum", "ResName", "AvgAngDiff",
                                         "StdAngDiff", "NumStruct")
  for i in range(StartRes, StopRes):
    print "%-8d %-8s %-10.4f %-10.4f %-10d" % (i+1, pRef.Seq[i], AvgDiff[i],
                                               StdDiff[i], NAdd[i])  


if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv[1:])
  protein.AssignBondsDflt = False

  PdbRef = Args["ARGS"][0]  
  pRef = protein.ProteinClass(Pdb = PdbRef)
  ValRef = [pRef.PhiPsi(i) for i in range(len(pRef))]
  ValRef = array(ValRef, float)
  NRes = len(pRef)
  AvgDiff = zeros(NRes, float)
  VarDiff = zeros(NRes, float)
  NAdd = zeros(NRes, float)
  
  #check for a trajectory
  if any([".crd" in x for x in Args["ARGS"]]):
    NSkip = int(Args.get("nskip", 0))
    NRead = int(Args.get("nread", -1))
    NStride = int(Args.get("nstride", 1))
    TrjFiles = Args["ARGS"][1:]
    for TrjFile in TrjFiles:
      PrmtopFile = TrjFile.replace(".mdtrj.crd.gz", ".prmtop.parm7")
      PrmtopFile = PrmtopFile.replace(".mdtrj.crd", ".prmtop.parm7")
      Trj = coords.TrjClass(TrjFile, PrmtopFile, NSkip = NSkip,
                            NRead = NRead, NStride = NStride)
      if len(Trj) <= 0: continue
      pTrj = protein.ProteinClass()
      pTrj.LinkTrj(Trj)
      
      Map = sequence.SeqMapClass(pRef.Seq, pTrj.Seq)

      PB = scripttools.ProgressBar("Reading %s" % TrjFile, Steps = len(Trj))
      ThisAvgDiff = zeros(NRes, float)
      ThisVarDiff = zeros(NRes, float)
      ThisNAdd = zeros(NRes, float)
      for (j,Pos) in enumerate(Trj):
        PB.Update(j)
        for i in range(len(Map)):
          Phi, Psi = pTrj.PhiPsi(i + Map.c)
          PhiDiff = geometry.NearestAngle(Phi - ValRef[i+Map.a,0], 0)
          PsiDiff = geometry.NearestAngle(Psi - ValRef[i+Map.a,1], 0)
          Diff = sqrt(PhiDiff**2 + PsiDiff**2)
          ThisAvgDiff[i+Map.a] += Diff
          ThisVarDiff[i+Map.a] += (Diff*Diff)
          ThisNAdd[i+Map.a] += 1
      PB.Clear()

      AvgDiff += ThisAvgDiff
      VarDiff += ThisVarDiff
      NAdd += ThisNAdd

      print "AVERAGES FOR %s:" % TrjFile
      Disp(pRef, ThisAvgDiff, ThisVarDiff, ThisNAdd)
      print ""

      DispFinal = len(TrjFiles) > 1
    pTrj.UnlinkTrj()
      
    
  else:
    
    #get pdbs
    Pdbs = [f for f in Args["ARGS"][1:] if os.path.isfile(f)]
    #check for files
    for f in [f for f in Args["ARGS"][1:] if not os.path.isfile(f)]:
      print "Could not find %s." % f
    N = len(Pdbs)
    Filenames = [os.path.basename(x) for x in Pdbs]
    if N <= 1:
      print "Nothing to compare."
      exit()

    PB = scripttools.ProgressBar("Reading PDBs", Steps = len(Pdbs))    
    for (j,Pdb) in enumerate(Pdbs):
      p = protein.ProteinClass(Pdb)
      Map = sequence.SeqMapClass(pRef.Seq, p.Seq)
      PB.Update(j)
      for i in range(len(Map)):
        Phi, Psi = p.PhiPsi(i + Map.c)
        PhiDiff = geometry.NearestAngle(Phi - ValRef[i+Map.a,0], 0)
        PsiDiff = geometry.NearestAngle(Psi - ValRef[i+Map.a,1], 0)
        Diff = sqrt(PhiDiff**2 + PsiDiff**2)
        AvgDiff[i+Map.a] += Diff
        VarDiff[i+Map.a] += (Diff*Diff)
        NAdd[i+Map.a] += 1

    DispFinal = True        

  if DispFinal:
    print "FINAL AVERAGES:"
    Disp(pRef, AvgDiff, VarDiff, NAdd)
    print ""
  

  
