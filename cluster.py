#!/usr/bin/env python

#LAST MODIFIED: 01-26-09

Usage = """Runs clustering on pdb files in a directory or on an amber trajectory.
Produces clustresults.txt, clust####.txt, and clust####.pdb files."

Usage     : cluster.py PDBFILES [OPTIONS]
   OR     : cluster.py TRJFILE PRMTOPFILE [OPTIONS]

PDBFILES  : list of pdb files
TRJFILE   : trajectory CRD file (can be gzipped)
PRMTOPFILE: PARM7 file
OPTIONS   : "--rmsd=X" to set the rmsd tolerance (default 2.0)
            "--nskip=X" number of configs in trajectory to skip (default is 0)
            "--nread=X" number of configs in trajectory to read; -1 is all (default -1)
            "--nstride=X" read configs every nstride frames (default is 1)
            "--allatom" to use all atom rmsd (default is backbone)
            "--alphacarbon" to use just alpha carbon (default is backbone)
            "--maxclust=X" maximum number of clusters; negative values will force all
                           configurations to a cluster (default 10)
            "--maxclustwork=X" maximum number of working clusters; default is no max
            "--maxiter=X" maximum number of iterations (default 10)
            "--prefix=X" to change the output prefix (default "clust")
            "--compres='1,3-5'" to use residues 1 and 3-5 to minimize RMSD
            "--calcres='1-5,9'" to calculate rmsd only for residues 1-5 and 9
"""

#check for instructions
import sys
if len(sys.argv) == 1:
  print Usage
  sys.exit()


import sys, glob, os
import rmsd, coords, scripttools


#======== COMMAND LINE RUNNING ========


if __name__ == '__main__':

  #parse arguments
  Args = scripttools.ParseArgs(sys.argv)
  NSkip = int(Args.get("nskip", 0))
  NRead = int(Args.get("nread", -1))
  NStride = int(Args.get("nstride", 1))
  MaxCluster = int(Args.get("maxclust", 10))
  MaxIter = int(Args.get("maxiter", 10))
  MaxClusterWork = int(Args.get("maxclusterwork", 0))
  if MaxClusterWork == 0: MaxClusterWork = None
  RmsdTol = float(Args.get("rmsd", 2.0))
  Prefix = Args.get("prefix", "clust")
  CompResInd = scripttools.GetNumList(Args.get("compres", None), Offset = -1)
  CalcResInd = scripttools.GetNumList(Args.get("calcres", None), Offset = -1)

  #decide mask
  if "allatom" in Args["FLAGS"]:
    Mask = coords.NoMask
    print "Using all atoms"
  elif "alphacarbon" in Args["FLAGS"]:
    Mask = coords.AlphaCarbonMask
    print "Using alpha carbons"
  else:
    Mask = coords.BackboneMask
    print "Using backbone atoms"

  #decide mode
  if ".pdb" in Args[1]:
    Mode = 0
    PdbFiles = Args["ARGS"][1:]
    cobj = coords.PdbListClass(PdbFiles, Mask = Mask)
  else:
    Mode = 1
    TrjFile, PrmtopFile = Args[1], Args[2]
    cobj = coords.TrjClass(TrjFile, PrmtopFile, Mask = Mask,
                           NRead = NRead, NSkip = NSkip, NStride = NStride)

  #examine any residue specific masks
  AtomInd, CompInd, CalcInd = rmsd.GetCoordsObjMasks(cobj, AtomMask = Mask,
                                                     CompResInd = CompResInd,
                                                     CalcResInd = CalcResInd)

  #run the cluster algorithm
  Pos, ClustNum, ClustWeights, ClustPop, ConfRmsd, ClustRmsd = rmsd.ClusterMSS(cobj,
    RmsdTol, MaxIter = MaxIter, MaxCluster = MaxCluster, MaxClusterWork = MaxClusterWork,
    Method = 0, CompInd = CompInd, CalcInd = CalcInd)
  Indices = cobj.GetIndices()

  if Mode == 0:
    #save output pdbs  
    for (i, Posi) in enumerate(Pos):
      pc = 100. * float(ClustWeights[i]) / float(sum(ClustWeights))
      fn = Prefix + "%02d_%.0fpc.pdb" % (i+1, pc)
      if -i in ClustNum:
        j = ClustNum.index(-i)
      else:
        j = 0
      coords.SavePdbCoordsPdb(Posi, fn, PdbFiles[j])
    #save results
    rmsd.SaveClustResults(Pos, ClustNum, ClustWeights, ClustPop, ConfRmsd, ClustRmsd,
      ConfIndices = Indices, Verbose = False, Prefix = Prefix)
    
  elif Mode == 1:
    #save output configs
    for (i, Posi) in enumerate(Pos):
      pc = 100. * float(ClustWeights[i]) / float(sum(ClustWeights))
      fn = Prefix + "%02d_%.0fpc.pdb" % (i+1, pc) 
      coords.SavePdbCoordsAmb(Posi, PrmtopFile, fn)
    #save results
    rmsd.SaveClustResults(Pos, ClustNum, ClustWeights, ClustPop, ConfRmsd, ClustRmsd,
      ConfIndices = Indices, Verbose = False, Prefix = Prefix)

