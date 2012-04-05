#!/usr/bin/env python

#LAST MODIFIED: 07-03-08


Usage = """Reports on consensus information in pdb files.

Usage     : consensus.py [OPTIONS] PDB1 PDB2 ...
            consensus.py [OPTIONS] TRJFILE PRMTOPFILE

PDB*      : input pdb files 
OPTIONS   : "--ssthresh=X" where X is frac pdbs for secondary struct. (dflt 0.8)
            "--contthresh=X" where X is frac pdbs for contacts (dflt 0.9)
            "--dihthresh=X" where X is maximum std dev for dihedrals (dflt 90)
            "--phobics" to only investigate phobic-phobic contacts
            "--runss" to run secondary structure analysis (can be slow)
            "--quiet" to suppress informational messages
            "--zaminputs" to produce inputdihedrals.txt and inputcontacts.txt 
            "--dihfconst=X" to set the dihedral force constant (dflt 4.0)
            "--contfconst=X" to set the contact force constant (dflt 0.5)
            "--scalecont" to scale contact restraints based on consensus
            "--nskip=X" number of configs in trajectory to skip (default 0)
            "--nread=X" number of configs in trajectory to read; -1 is all (default -1)
            "--nstride=X" read configs every nstride frames (default 1)
"""
#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 2:
  print Usage
  sys.exit()


  
from numpy import *
import os, sys
import protein, scripttools, sequence, coords


SSThresh = 0.80
ContThresh = 0.9
DihThresh = 90.
DihTol = 1.25
PhobicsOnly = False


def GetCommon(CObj, RunSS = False, Verbose = False):
  """Returns residues with common secondary structure, common dihedrals,
  and common contacts."""
  if Verbose: print "Getting common secondary structures and contacts..."
  SSList, PhiPsiList, PhiPsiSqList = [], [], []
  CM = None
  NFrame = len(CObj)
  #run through all of the pdbs and tabulate the data
  CObj.Reset()
  p = protein.ProteinClass()
  p.LinkCoordsObj(CObj)
  Seq = p.Seq
  SeqLen = len(Seq)
  for (i, Pos) in enumerate(CObj):
    if Verbose:
      if "PdbFileList" in CObj.__dict__:
        print "Examining pdb %s" % CObj.PdbFileList[i]
      else:
        if (i+1) % 100 == 0: print "Examined %d frames" % (i+1)
    if RunSS: SSList.append(p.SecondaryStructure())
    ThisPhiPsiList = array([p.PhiPsi(i) for i in range(len(p))], float)
    PhiPsiList.append(ThisPhiPsiList)
    if CM is None:
      CM = p.ResContactMap()
    else:
      CM += p.ResContactMap()
  #find the most commond H or E motif at each point along the chain
  if RunSS:
    SSVals = ["H", "E"]
    ResSS, ResSSFrac = "", []
    for i in range(SeqLen):
      s = [SS[i] for SS in SSList]
      nVals = [s.count(x) for x in SSVals]
      ind = argmax(nVals)
      f = float(nVals[ind]) / float(NFrame)
      ResSSFrac.append(f)
      if f >= SSThresh:
        ResSS += SSVals[ind]
      else:
        ResSS += "-"
  else:
    ResSS = "-" * SeqLen
    ResSSFrac = [0.]*SeqLen
  #now find average phi-psi angles for each one
  FixedDih = []
  for i in range(SeqLen):
    ThisResPhiPsi = array([PhiPsi[i] for PhiPsi in PhiPsiList], float)
    s = std(ThisResPhiPsi, axis=0)
    a = mean(ThisResPhiPsi, axis=0)
    if all(s <= DihThresh):
      FixedDih.append((i, a[0], a[1], s[0], s[1]))
  #now compute contacts
  CM = CM / NFrame
  Contacts = []
  SkipMask = [not sequence.Hydrophobic(x) and PhobicsOnly for x in Seq]
  for i in range(SeqLen):
    if SkipMask[i]: continue
    for j in range(i+3, SeqLen):
      if SkipMask[j]: continue
      if CM[i,j] >= ContThresh:
        Contacts.append((i,j,CM[i,j]))
  #summarize
  nums = "1234567890" * (SeqLen/10 + 1)
  nums = nums[:SeqLen]
  if RunSS:
    s = "CONSENSUS SECONDARY STRUCTURE:\n%s\n%s" % (nums, ResSS)
  else:
    s = ""
  s += "\n\nCONSENSUS DIHEDRALS:\nRes PhiAvg PsiAvg PhiStd PsiStd\n"
  s += "\n".join(["%s%d %8.2f %8.2f %8.2f %8.2f" % (Seq[i],i+1,x,y,z,w)
                  for (j, (i,x,y,z,w)) in enumerate(FixedDih)])
  s += "\n\nCONSENSUS CONTACTS:\nResi Resj FracMade\n"
  s += "\n".join(["%s%d %s%d %8.3f" % (Seq[i],i+1,Seq[j],j+1,f) for (i,j,f) in Contacts])
  if Verbose: print ""
  print s
  #get sequence
  Seq = sequence.Standardize(Seq)
  return Seq, Contacts, ResSS, ResSSFrac, FixedDih


def WriteZamFiles(ZamFolder, FixedDih, Contacts, DihFConst, ContFConst, ScaleCont):
    #first output dihedrals
    s = "#Resnum Phi Psi FConst PhiTol PsiTol\n"
    for (i,aPhi,aPsi,sPhi,sPsi) in FixedDih:
      PhiTol, PsiTol = 0.5 * DihTol * sPhi, 0.5 * DihTol * sPsi
      s += "%d %8.2f %8.2f %8.3f %8.3f %8.3f\n" % (i+1,aPhi,aPsi,DihFConst,PhiTol,PsiTol)
    file(os.path.join(ZamFolder, "inputdihedrals.txt"),"w").write(s)
    #now output contacts
    MaxFrac = array([f for (i,j,f) in Contacts], float).max()
    s = "#Resi Resj FConst\n"
    for (i,j,f) in Contacts:
      FConst = ContFConst
      if ScaleCont: FConst *= (f / MaxFrac)
      s += "%d %d %8.3f\n" % (i+1,j+1,FConst)
    file(os.path.join(ZamFolder, "inputcontacts.txt"),"w").write(s)

#command line running
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv, {"contthresh":ContThresh,
                                          "dihthresh":DihThresh,
                                          "ssthresh":SSThresh,
                                          "dihfconst":4.0, "contfconst":0.5,
                                          "nskip":0, "nread":-1, "nstride":1})
  #get options
  ContThresh = Args["contthresh"]
  DihThresh = Args["dihthresh"]
  SSThresh = Args["ssthresh"]
  NSkip, NRead, NStride = Args["nskip"], Args["nread"], Args["nstride"]
  WriteZam = "writezam" in Args["FLAGS"]
  Verbose = not "quiet" in Args["FLAGS"]
  ScaleDih = "scaledih" in Args["FLAGS"]
  ScaleCont = "scalecont" in Args["FLAGS"]
  DihFConst = Args["dihfconst"]
  ContFConst = Args["contfconst"]
  ZamInputs = "zaminputs" in Args["FLAGS"]
  PhobicsOnly = "phobics" in Args["FLAGS"]
  RunSS = "runss" in Args["FLAGS"]
  #check the input files
  Files = [x for x in Args["ARGS"][1:] if os.path.isfile(x)]
  if all([".pdb" in x for x in Files]):
    CObj = coords.PdbListClass(Files)
  else:
    CObj = coords.TrjClass(Files[0], Files[1], NSkip = NSkip, NRead = NRead, NStride = NStride)
  #get common parameters from consensus
  Seq, Contacts, ResSS, ResSSFrac, FixedDih = GetCommon(CObj, RunSS, Verbose = Verbose)
  #write zam output
  if ZamInputs:
    WriteZamFiles(".", FixedDih, Contacts, DihFConst, ContFConst, ScaleCont)


  

        
      
  