#!/usr/bin/env python

#LAST MODIFIED: 12-19-08

Usage = """Calculates rmsd between pdb files.

Usage     : rmsd.py [OPTIONS] PDBFILE1 PDBFILE2 ...
        OR: rmsd.py [OPTIONS] --traj [--avg=n] [--frac=X] REFPDBFILE TRJFILE PRMTOPTILE 

OPTIONS   : "--first" will only compare the first to all other pdb files
            "--align" to align the sequences before rmsd comparison
            "--compres='1,3-5'" to use residues 1 and 3-5 to minimize RMSD
            "--calcres='1-5,9'" to calculate rmsd only for residues 1-5 and 9
            "--frac=X" to compute fraction of structures within X angstroms
            "--nskip=X" number of configs in trajectory to skip (default is 0)
            "--nread=X" number of configs in trajectory to read; -1 is all (default -1)
            "--nstride=X" read configs every nstride frames (default is 1)
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()

from numpy import *  
import copy, os, coords, random
import geometry, sequence, protein, scripttools


def RMSD(Pos1, Pos2, Align = False, Center = True,
         CompInd = None, CalcInd = None, Verbose = False,
         RetAlignment = False):
  """Calculates the RMSD between two conformations.
* Pos1: array of dimensions [N,3] for conformation 1
* Pos2: array of dimensions [N,3] for conformation 2
* Align: True = modify Pos2 to be aligned to Pos1
  (default is False)
* CompInd: indices in [0,N) for positions in Pos to perform alignment
* CalcInd: indices in [0,N) for positions in Pos to compute RMSD
* Verbose: true to report shape/size mismatches
* RetAlignment: True to return the translation vectors and rotation matrix"""
  #clean indices
  AllInd = arange(len(Pos1), dtype = int)
  if all(CompInd == AllInd): CompInd = None
  if all(CalcInd == AllInd): CalcInd = None
  #get indices
  if CompInd is None:
    p1, p2, n = Pos1, Pos2, len(Pos1)
  else:
    p1, p2, n = Pos1.take(CompInd, axis=0), Pos2.take(CompInd, axis=0), len(CompInd)
  #check for correct shapes
  if not shape(p1) == shape(p2):
    if Verbose: print "Position vectors are not the same size."
    return
  elif not len(shape(p1)) == 2:
    if Verbose: print "Position vectors are not the correct rank."
    return
  #get alignment
  Pos1Vec, Pos2Vec, RotMat, Resid = geometry.AlignmentRMSD(p1, p2, Center = Center)
  #compute rmsd
  if not all(CompInd == CalcInd):
    if CalcInd is None:
      p1, p2, n = Pos1, Pos2, len(Pos1)
    else:
      p1, p2, n = Pos1.take(CalcInd, axis=0), Pos2.take(CalcInd, axis=0), len(CalcInd)
    p1 = p1 + Pos1Vec
    p2 = dot(p2 + Pos2Vec, RotMat)
    Resid = sum((p1 - p2)**2, axis=None)
  r = sqrt(Resid / float(n))
  #align Pos2 to Pos1
  if Align: Pos2[:,:] = dot(Pos2 + Pos2Vec, RotMat) - Pos1Vec
  if RetAlignment:
    return r, Pos1Vec, Pos2Vec, RotMat
  else:
    return r


def GetProteinClassMasks(p, AtomMask = None, CompResInd = None,
                         CalcResInd = None):
  if AtomMask is None:
    AtomInd = array(range(len(p.Atoms)), int)
  else:
    AtomInd = p.AtomInd(AtomName = AtomMask)
  CompInd, CalcInd = None, None
  if not CompResInd is None:
    CompInd = p.AtomInd(AtomName = AtomMask, ResNum = CompResInd)
    CompInd = [i for (i,j) in enumerate(AtomInd) if j in CompInd]
  if not CalcResInd is None:
    CalcInd = p.AtomInd(AtomName = AtomMask, ResNum = CalcResInd)
    CalcInd = [i for (i,j) in enumerate(AtomInd) if j in CalcInd]
  return AtomInd, CompInd, CalcInd

def GetProteinClassMasks2(p1, p2, AtomMask = None,
                          CompResInd = None, CalcResInd = None):
  AtomInd1, AtomInd2, CompInd, CalcInd = [], [], [], []
  count = 0
  for (i, a1) in enumerate(p1.Atoms):
    ResNum = p1.AtomResNum[i]
    InComp, InCalc = True, True
    if not CompResInd is None: InComp = ResNum in CompResInd
    if not CalcResInd is None: InCalc = ResNum in CalcResInd
    j = p2.AtomNum(ResNum, a1.Name, NotFoundError = False)
    if j >= 0:
      AtomInd1.append(i)
      AtomInd2.append(j)
      if InComp: CompInd.append(count)
      if InCalc: CalcInd.append(count)
      count += 1
  AtomInd1 = array(AtomInd1, int)
  AtomInd2 = array(AtomInd2, int)
  CompInd = array(CompInd, int)
  CalcInd = array(CalcInd, int)
  return AtomInd1, AtomInd2, CompInd, CalcInd


def GetCoordsObjMasks(c, AtomMask = None, CompResInd = None,
                      CalcResInd = None):
  if AtomMask is None:
    AtomInd = array(range(len(c.AtomNames)), int)
  else:
    AtomInd = [i for (i, an) in enumerate(c.AtomNames)
               if an.strip() in AtomMask]
    AtomInd = array(AtomInd, int)
  CompInd, CalcInd = None, None
  if not CompResInd is None:
    CompInd = [i for (i, an) in enumerate(c.AtomNames)
               if an.strip() in AtomMask and c.AtomRes[i] in CompResInd]
    CompInd = [i for (i,j) in enumerate(AtomInd) if j in CompInd]
  if not CalcResInd is None:
    CalcInd = [i for (i, an) in enumerate(c.AtomNames)
               if an.strip() in AtomMask and c.AtomRes[i] in CalcResInd]
    CalcInd = [i for (i,j) in enumerate(AtomInd) if j in CalcInd]
  return AtomInd, CompInd, CalcInd


def RMSDProteinClass(p1, p2, Center = True, Backbone = True, AlignSeq = False,
                     CompResInd = None, CalcResInd = None, UpdateBFactors = False,
                     AlignAtoms = False):
  "Returns (RMSD,NRes)"
  #align the sequences
  p2orig = p2
  if AlignSeq:
    Map = sequence.SeqMapClass(p1.Seq, p2.Seq)
    p1 = p1[Map.a:Map.b]
    p2 = p2[Map.c:Map.d]
  #filter res indices
  n = min(len(p1), len(p2))
  if not CompResInd is None: CompResInd = [i for i in CompResInd if i < n]
  if not CalcResInd is None: CalcResInd = [i for i in CalcResInd if i < n]
  #get the right atoms
  if Backbone:
    AtomMask = ["N", "CA", "C"]
  else:
    AtomMask = None
  if AlignAtoms:
    AtomInd1, AtomInd2, CompInd, CalcInd = GetProteinClassMasks2(p1, p2,
      AtomMask = AtomMask, CompResInd = CompResInd, CalcResInd = CalcResInd)

  else:
    AtomInd1, CompInd, CalcInd = GetProteinClassMasks(p1, AtomMask = AtomMask,
                                                      CompResInd = CompResInd,
                                                      CalcResInd = CalcResInd)
    AtomInd2, CompInd, CalcInd = GetProteinClassMasks(p2, AtomMask = AtomMask,
                                                      CompResInd = CompResInd,
                                                      CalcResInd = CalcResInd)
  #filter positions
  Pos1, Pos2 = p1.Pos.take(AtomInd1, axis=0), p2.Pos.take(AtomInd2, axis=0).copy()
  #compute rmsd
  r = RMSD(Pos1, Pos2, Align = True, Center = Center,
           CompInd = CompInd, CalcInd = CalcInd)
  #determine num of residues used in rmsd calculation
  if CalcResInd is None:
    NRes = len(p1)
  else:
    NRes = len(CalcResInd)
  #see if we need to update bfactors
  if UpdateBFactors:
    if CalcInd is None: CalcInd = range(len(AtomInd2))
    Shift = p2orig.Res[Map.c].StartAtom
    for i in CalcInd:
      an = AtomInd2[i] + Shift
      p2orig.Atoms[an].BFactor = sqrt(sum((Pos1[i] - Pos2[i])**2))   
  return r, NRes


def RMSDPdb(PdbFile1, PdbFile2, Center = True, Backbone = True,
            AlignSeq = False, CompResInd = None, CalcResInd = None):
  p1 = protein.ProteinClass(Pdb = PdbFile1)
  p2 = protein.ProteinClass(Pdb = PdbFile2)
  return RMSDProteinClass(p1, p2, Center = Center, Backbone = Backbone,
                          AlignSeq = AlignSeq, CompResInd = CompResInd,
                          CalcResInd = CalcResInd)

def ClusterMSS(CoordsObj, Cutoff, MaxIter = 3, MaxCluster = None,
  MaxClusterWork = None, Method = 0, CompInd = None, CalcInd = None,
  Weights = None, Verbose = True, IterMaxCluster = False,
  IterNormalize = False):
  """Clusters conformations in a trajectory based on RMSD distance.
* CoordsObj: an object exposing the functions GetNextCoords() which
  returns an array object of the next set of coordinates (or None
  if at the end) and has optional list variable Mask, and
  Reset() which moves back to the first set of coordinates
* Cutoff: maximum RMSD distance of a configuration to a cluster
* MaxIter: maximum number of iterations to perform
* MaxCluster: maximum number of clusters; negative values will force
  all configs to be members of a cluster (default is none)
* MaxClusterWork: maximum numeber of working clusters
* Method: 0 to cluster such that each RMSD between a configuration
  and the average cluster configuration is below Cutoff; 1 is
  same except no alignment is performed
* CompInd: indices in [0,N) for positions in Pos to perform alignment
* CalcInd: indices in [0,N) for positions in Pos to compute RMSD
* Weights: weighting factor for each conformation
* IterMaxCluster: True will dump all but MaxCluster configs each iter
* IterNormalize: True will dump previous iter contribs to centroids
"""
  def BasicRMSD(Pos1, Pos2, Ind = None):
    "Calculates the rmsd between two configurations without alignment."
    if Ind is None:
      rmsdsq = sum((Pos1-Pos2)**2) / float(size(Pos1,0))
    else:
      rmsdqs = sum((Pos1[Ind]-Pos2[Ind])**2) / float(size(Pos1[Ind],0))
    rmsdsq = max([rmsdsq,0.])
    return sqrt(rmsdsq)
  Iteration = 0   #iteration number
  WeightSum = []  #total weights of clusters
  PosSum = []        #list of cluster configuration arrays
  FinalIters = 0  #number of iterations without additions/deletions of clusters
  NCoord = len(CoordsObj)
  if Weights is None: Weights = ones(NCoord, float)
  Weights = array(Weights, float)
  if not len(Weights) == NCoord:
    raise IndexError, "Incorrect number of array elements in Weights."
  #filter weights for too low values
  Weights = Weights.copy()
  Weights = Weights / Weights.max()
  Weights[Weights < 1.e-100] = 0.
  NFrameTot = int(sum(Weights > 0))
  StartInd, NewStartInd = 0, -1
  print "Using %d conformations in %d trajectory frames." % (NFrameTot, len(CoordsObj))
  while FinalIters < 2:
    Iteration += 1
    FinalIters += 1
    if Iteration > MaxIter:
      if Verbose: print "Did not converge within maximum number of iterations"
      break
    if Verbose: print "Cluster iteration %d" % Iteration
    if Verbose: print "Starting with %d clusters" % len(PosSum)
    CoordsObj.Reset()
    ClustNum = zeros(NCoord, int)  #cluster number of each configuration, starting at 1
    NAddThis = [0]*len(PosSum)     #number of configs added to each cluster this iteration
    ThisFrame = 0
    PosSumThis = copy.deepcopy(PosSum)
    WeightSumThis = copy.deepcopy(WeightSum)
    #check where to start
    if NewStartInd >= 0: StartInd = NewStartInd
    NewStartInd = -1
    for CurInd in range(StartInd, NCoord) + range(0, StartInd):
      CurPos = CoordsObj[CurInd]
      #check for zero weight
      CurWeight = Weights[CurInd]
      if CurWeight == 0.:
        ClustNum[CurInd] = 0
        continue
      ThisFrame += 1
      ind = -1  #cluster number assigned to this config; -1 means none
      #calculate the rmsd between this configuration and each cluster config,
      #but stop when a rmsd is found which is below the cutoff
      minRMSD = 1.e300
      for (i, PosSumi) in enumerate(PosSum):
        ThisPos = PosSumi / WeightSum[i]
        if Method == 0:
          #rmsd between new config and average cluster config
          r = RMSD(ThisPos, CurPos, Align = True, Center = True,
                   CompInd = CompInd, CalcInd = CalcInd)
        else:
          #rmsd between new config and average cluster config, without alignment
          r = BasicRMSD(ThisPos, CurPos, Ind = CalcInd)
        minRMSD = min(minRMSD, r)
        if r < Cutoff:
          #go with a cluster if rmsd is within the cutoff
          ind = i
          break
      if ind >= 0:
        #add the configuration to the cluster
        PosSum[ind] = PosSum[ind] + CurPos * CurWeight
        WeightSum[ind] = WeightSum[ind] + CurWeight
        NAddThis[ind] = NAddThis[ind] + 1
        ClustNum[CurInd] = ind+1
      elif len(PosSum) < MaxClusterWork or MaxClusterWork is None:
        #create a new cluster with this config, as long as it
        #doesn't exceed the maximum number of working clusters
        if minRMSD == 1.e300: minRMSD = 0.
        if Verbose: print "Adding cluster: config %d (%d/%d) | min RMSD %.1f | %d clusters tot" % (CoordsObj.Index+1,
                          ThisFrame, NFrameTot, minRMSD, len(PosSum)+1)
        PosSum.append(CurPos * CurWeight)
        WeightSum.append(CurWeight)
        NAddThis.append(1)
        ClustNum[CurInd] = len(PosSum)
        FinalIters = 0
      else:
        #cluster is nothing
        ClustNum[CurInd] = 0
        FinalIters = 0
        if NewStartInd < 0:
          NewStartInd = CurInd
          if Verbose: print "Ran out of clusters. Next iteration starting from config %d" % (CoordsObj.Index+1,)
    #remove contribution to centroids from all but this round
    if IterNormalize:
      for i in range(len(PosSumThis)):
        PosSum[i] = PosSum[i] - PosSumThis[i]
        WeightSum[i] = WeightSum[i] - WeightSumThis[i]
    del PosSumThis
    del WeightSumThis
    #loop through clusters
    i = 0
    while i < len(PosSum):
      #remove clusters that have no additions this iteration
      if NAddThis[i] == 0:
        if Verbose: print "Removing cluster %d" % (i+1,)
        del PosSum[i]
        del WeightSum[i]
        del NAddThis[i]
        for (k, cn) in enumerate(ClustNum):
          if cn > i + 1:
            ClustNum[k] -= 1
          elif cn == i + 1:
            ClustNum[k] = -1
        FinalIters = 0
      else:
        i += 1
    #sort clusters and then remove any beyond MaxCluster
    ClustNum = array(ClustNum, int)
    PosSum, ClustNum, WeightSum, NAddThis = __SortClust(PosSum, ClustNum, Weights, WeightSum, NAddThis, Verbose)
    #crop off any extraneous clusters; clusterless configs
    #are assigned a cluster index of 0
    if IterMaxCluster and not MaxCluster is None and len(PosSum) > abs(MaxCluster):
      del PosSum[abs(MaxCluster):]
      WeightSum = WeightSum[:abs(MaxCluster)]
      NAddThis = NAddThis[:abs(MaxCluster)]
      ClustNum[abs(ClustNum) > abs(MaxCluster)] = 0
  #crop off any extraneous clusters; clusterless configs
  #are assigned a cluster index of 0
  if not IterMaxCluster and not MaxCluster is None and len(PosSum) > abs(MaxCluster):
    del PosSum[abs(MaxCluster):]
    WeightSum = WeightSum[:abs(MaxCluster)]
    NAddThis = NAddThis[:abs(MaxCluster)]
    ClustNum[abs(ClustNum) > abs(MaxCluster)] = 0  
  #finalize things
  if Verbose: print "Calculating average structures"
  Pos = [x / y for (x, y) in zip(PosSum, WeightSum)]
  del PosSum
  del WeightSum
  #get cluster populations
  ClustWeights, ClustPop = __CalcClustPop(Pos, ClustNum, Weights)
  #if there is a maximum cluster specification that's negative, force
  #everything to the closest cluster
  if not MaxCluster == None and MaxCluster < 0:
    ClustWeights, ClustPop, Pos, ClustNum = __CalcForceClust(CoordsObj, ClustWeights, ClustPop, Pos, ClustNum, Weights,
                                               CompInd, CalcInd, Verbose)
  #calculate final rmsd values for configs and clusters
  Pos, ClustNum, ConfRmsd, ClustRmsd = __CalcRmsd(CoordsObj, Pos, ClustNum, CompInd, CalcInd, Verbose)
  if Verbose: print "%d configurations sorted into %d clusters" % (NFrameTot, len(Pos))
  return Pos, ClustNum, ClustWeights, ClustPop, ConfRmsd, ClustRmsd


def __SortClust(Pos, ClustNum, Weights, WeightSum, NAddThis, Verbose = True):
  if Verbose: print "Reordering clusters by population"
  Sums = [(sum(Weights[abs(ClustNum) == i+1]), i) for i in range(len(Pos))]
  Sums.sort()
  Sums.reverse()
  #resort the arrays
  Pos = [Pos[j] for (i,j) in Sums]
  WeightSum = [WeightSum[j] for (i,j) in Sums]
  NAddThis = [NAddThis[j] for (i,j) in Sums]
  #create a dictionary which will tell us the new cluster
  #number for a given old cluster number
  Trans = {0:0}
  for i in range(len(Pos)):
    ind = Sums[i][1] + 1
    Trans[ind] = i + 1
    Trans[-ind] = -i - 1
  #update ClustNum with the rearranged cluster numbers
  ClustNum = array([Trans[i] for i in ClustNum], int)
  return Pos, ClustNum, WeightSum, NAddThis
  

def __CalcClustPop(Pos, ClustNum, Weights):
  #update the cluster population
  ClustWeights = array([sum(Weights[abs(ClustNum) == i+1]) for i in range(len(Pos))], float)
  #update the populations
  ClustPop = [float(sum(abs(ClustNum) == i+1)) for i in range(len(Pos))]
  ClustPop = array(ClustPop, float)
  return ClustWeights, ClustPop


def __CalcForceClust(CoordsObj, ClustWeights, ClustPop, Pos, ClustNum, Weights,
                     CompInd = None, CalcInd = None, Verbose = True):
  #count the number of clusterless configurations
  c = sum(ClustNum == 0)
  if Verbose: print "Forcing %d extraneous configurations to existing clusters" % c
  #find the nearest cluster to each clusterless config and assign it
  CoordsObj.Reset()
  for (j, CurPos) in enumerate(CoordsObj):
    if ClustNum[j] == 0:
      ind = -1
      minr = 0.
      for (i, Posi) in enumerate(Pos):
        r = RMSD(Posi, CurPos, Align = False, Center = True,
                 CompInd = CompInd, CalcInd = CalcInd)
        if r < minr or ind < 0: 
          ind = i
          minr = r
      ClustNum[j] = ind + 1
      ClustWeights[ind] = ClustWeights[ind] + Weights[j]
      ClustPop[ind] = ClustPop[ind] + 1.
  return ClustWeights, ClustPop, Pos, ClustNum
  

def __CalcRmsd(CoordsObj, Pos, ClustNum, 
               CompInd = None, CalcInd = None, Verbose = True):
  if Verbose: print "Calculating cluster rmsd values"
  #calculate the pairwise cluster rmsd values
  ClustRmsd = zeros((len(Pos), len(Pos)),float)
  for (i, Posi) in enumerate(Pos):
    for (j, Posj) in enumerate(Pos):
      if j <= i: continue
      ClustRmsd[i,j] = RMSD(Posi, Posj, Align=False, Center=True,
                            CompInd = CompInd, CalcInd = CalcInd)
      ClustRmsd[j,i] = ClustRmsd[i,j]
  if Verbose: print "Calculating final rmsd values"
  #loop through configs and find the one with the lowest
  #rmsd in each cluster
  ConfRmsd = -1. * ones(len(ClustNum), float)
  MinRmsd = [-1]*len(Pos)
  CoordsObj.Reset()
  for (CurInd, CurPos) in enumerate(CoordsObj):
    i = abs(ClustNum[CurInd]) - 1
    if i >= 0:
      ConfRmsd[CurInd] = RMSD(Pos[i], CurPos, Align=False, Center=True,
                              CompInd = CompInd, CalcInd = CalcInd)
      if MinRmsd[i] < 0:
        MinRmsd[i] = CurInd
      elif ConfRmsd[MinRmsd[i]] > ConfRmsd[CurInd]:
        MinRmsd[i] = CurInd
  #loop through the configs again and extract the
  #coords of the minimum-rmsd configs for each clust
  if Verbose: print "Finding nearest cluster structures"
  for (i, ind) in enumerate(MinRmsd):
    ClustNum[ind] = -ClustNum[ind]
    Pos[i] = CoordsObj.Get(ind, coords.NoMask)
  return Pos, ClustNum, ConfRmsd, ClustRmsd


def SaveClustResults(Pos, ClustNum, ClustWeights, ClustPop, ConfRmsd, ClustRmsd,
  Prefix = "clust", ConfIndices = None, Verbose = False):
  #make the indices
  if ConfIndices is None:
    ConfIndices = range(0, len(ConfRmsd))
  #calculate the percent
  x = sum(ClustWeights)
  ClustPct = [100.*w / x for w in ClustWeights]
  #save results file
  s = "CLUSTER POPULATION:\nCluster number, Number of configs, Percent\n"
  s += "\n".join(["%-5d  %-7d  %.2f" % (i+1, ClustPop[i], ClustPct[i])
                  for i in range(0,len(Pos))])
  s += "\n\nCLUSTER-TO-CLUSTER RMSD:\nCluster number, Cluster number, RMSD\n"
  for i in range(0,len(Pos)):
    for j in range(i+1,len(Pos)):
      s += "%-5d  %-5d  %-8.2f\n" % (i+1, j+1, ClustRmsd[i,j])
  s += "\n\nCLUSTER CONFIG NUMBERS:\nCluster number, Config number, RMSD\n"
  ClustConf = [(abs(cn), i) for (i, cn) in enumerate(ClustNum) if cn < 0]
  ClustConf.sort()
  for (cn, i) in ClustConf:
    s += "%-5d  %-5d  %-8.2f\n" % (cn, ConfIndices[i], ConfRmsd[i])
  s += "\n\nCLUSTER MEMBERSHIP:\nConfig number, Cluster number, RMSD\n"
  for (i, r) in enumerate(ConfRmsd):
    if r < 0.:
      s += "%-7d  %-5s  %-8s\n" % (ConfIndices[i], "--", "--")
    else:
      s += "%-7d  %-5d  %-8.2f\n" % (ConfIndices[i], ClustNum[i], ConfRmsd[i])
  fn = Prefix + "results.txt"
  file(fn,"w").write(s)
  #save cluster files
  if Verbose:
    for i in range(0,len(Pos)):
      data = [(ConfRmsd[j], ConfIndices[j]) for j in range(0,len(ConfRmsd))
              if abs(ClustNum[j]) == i + 1]
      data.sort()
      s = "\n".join(["%-7d  %-8.2f" % (j[1],j[0]) for j in data]) + "\n"
      fn = Prefix + "%04d-%dpc.txt" % (i+1, int(ClustPct[i]))
      file(fn, "w").write(s)



#======== COMMAND-LINE RUNNING ========

def GetResList(Arg):
  if Arg is None or Arg == "": return None
  ResList = []
  for s in """'"[]""":
    Arg = Arg.replace(s,"")
  for l in Arg.split(","):
    if "-" in l:
      a, b = [int(x) for x in l.split("-")]
      ResList.extend(range(a-1, b))
    else:
      a = int(l)
      ResList.append(a - 1)
  ResList.sort()
  ResList = [x for (i,x) in enumerate(ResList) if not x in ResList[i+1:]]
  return ResList

def RepRMSD(r):
  if not type(r) is list: r = [r]
  if None in r:
    return 'NA'
  else:
    return "%-8.2f" % mean(r)

def ClipRMSD(r):
  if r is None:
    return 1.e300
  else:
    return r


if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv[1:])
  #get options
  Align = "align" in Args["FLAGS"]
  CompResInd = GetResList(Args.get("compres", None))
  CalcResInd = GetResList(Args.get("calcres", None))
  
  #check for a trajectory
  if "traj" in Args["FLAGS"]:
    PdbRef, TrjFile, PrmtopFile = Args["ARGS"][:3]
    NAvg = int(Args.get("avg", 1))
    Cut = float(Args.get("frac", 0.))
    NSkip = int(Args.get("nskip", 0))
    NRead = int(Args.get("nread", -1))
    NStride = int(Args.get("nstride", 1))
    Trj = coords.TrjClass(TrjFile, PrmtopFile, NSkip = NSkip,
                          NRead = NRead, NStride = NStride)
    pTrj = protein.ProteinClass()
    pTrj.LinkTrj(Trj)
    pRef = protein.ProteinClass(Pdb = PdbRef)
    print "%-10s %-8s %-8s %-5s" % ("Frame","BB_RMSD", "All_RMSD", "NRes")
    i = 0
    y1, y2 = [], []
    z1, z2 = [], []
    for Pos in Trj:
      i += 1
      x1, NRes = RMSDProteinClass(pRef, pTrj, Backbone = True, AlignSeq = Align,
                                  CompResInd = CompResInd, CalcResInd = CalcResInd)
      x2, NRes = RMSDProteinClass(pRef, pTrj, Backbone = False, AlignSeq = Align,
                                  CompResInd = CompResInd, CalcResInd = CalcResInd)
      y1.append(x1)
      y2.append(x2)
      z1.append(ClipRMSD(x1))
      z2.append(ClipRMSD(x2))
      if i % NAvg == 0:
        print "%-10d %-8s %-8s %-5d" % (Trj.Index+1, RepRMSD(y1), RepRMSD(y2), NRes)
        y1, y2 = [], []
    pTrj.UnlinkTrj()
    if Cut > 0.:
      z1 = array(z1, float)
      z2 = array(z2, float)
      frac1 = sum(z1 <= Cut) / float(len(z1))
      frac2 = sum(z2 <= Cut) / float(len(z2))
      print "\nFraction of trajectory within %.2f A:" % Cut
      print "  %.5f  backbone" % frac1
      print "  %.5f  all-atom" % frac2
    
  else:
    
    #get pdbs
    Pdbs = [f for f in Args["ARGS"] if os.path.isfile(f)]
    #check for files
    for f in [f for f in Args["ARGS"] if not os.path.isfile(f)]:
      print "Could not find %s." % f
    N = len(Pdbs)
    Filenames = [os.path.basename(x) for x in Pdbs]
    if N <= 1:
      print "Nothing to compare."
      sys.exit()
    if "first" in Args["FLAGS"]:
      MaxLen1 = len(Filenames[0])
      MaxLen2 = max([len(fn) for fn in Filenames[1:]])
      print "%-*s %-*s %-8s %-8s %-5s" % (MaxLen1, "Pdb1", MaxLen2, "Pdb2",
                                          "BB_RMSD", "All_RMSD", "NRes")
      p1 = protein.ProteinClass(Pdb = Pdbs[0])
      for j in range(1, N):
        p2 = protein.ProteinClass(Pdb = Pdbs[j])
        x1, NRes = RMSDProteinClass(p1, p2, Backbone = True, AlignSeq = Align,
                                    CompResInd = CompResInd, CalcResInd = CalcResInd)
        x2, NRes = RMSDProteinClass(p1, p2, Backbone = False, AlignSeq = Align,
                                    CompResInd = CompResInd, CalcResInd = CalcResInd)
        print "%-*s %-*s %-8s %-8s %-5d" % (MaxLen1, Filenames[0],
                                            MaxLen2, Filenames[j],
                                            RepRMSD(x1), RepRMSD(x2), NRes)

    else:

      MaxLen = max([len(fn) for fn in Filenames])
      print "%-*s %-*s %-8s %-8s %-5s" % (MaxLen, "Pdb1", MaxLen, "Pdb2",
                                          "BB_RMSD", "All_RMSD", "NRes")
      Prots = [protein.ProteinClass(Pdb = f) for f in Pdbs]
      for i in range(0, N):
        for j in range(i+1, N):
          x1, NRes = RMSDProteinClass(Prots[i], Prots[j], Backbone = True, AlignSeq = Align,
                                      CompResInd = CompResInd, CalcResInd = CalcResInd)
          x2, NRes = RMSDProteinClass(Prots[i], Prots[j], Backbone = False, AlignSeq = Align,
                                      CompResInd = CompResInd, CalcResInd = CalcResInd)
          print "%-*s %-*s %-8s %-8s %-5d" % (MaxLen, Filenames[i],
                                              MaxLen, Filenames[j],
                                              RepRMSD(x1), RepRMSD(x2), NRes)
        
        