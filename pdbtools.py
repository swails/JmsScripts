#!/usr/bin/env python

#last modified: 11-25-08

#CONVENTION: Any command that takes Pdb as an argument can use
#either the pdb filename or the actual pdb data 

#CONVENTION: Residue and atom numbers start as 1 as in the Pdb file



Usage = """Modifies pdb files.

Usage     : pdbtools.py cap INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py decap INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py dehydrogen INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py extendc SEQ INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py extendn SEQ INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py extend SEQBEFORE SEQAFTER INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py generate SEQ OUTPDB
   OR     : pdbtools.py renumber INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py center INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py rotate RESNUM PHI PSI INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py showseq INPDB
   OR     : pdbtools.py showoverlaps INPDB [OVERLAPDIST]
   OR     : pdbtools.py splice INPDB1 INPDB2 OUTPDB
   OR     : pdbtools.py spliceopt INPDB1 INPDB2 OUTPDB
   OR     : pdbtools.py trim STARTRES STOPRES INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py standardize INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py alignseq INPDB1 INPDB2
   OR     : pdbtools.py minimize INPDB [--out OUTPDB | --suffix=X] [--steps1=X] [--steps2=X]
   OR     : pdbtools.py longminimize INPDB [--out OUTPDB | --suffix=X]
   OR     : pdbtools.py getenergy INPDB
   OR     : pdbtools.py dssp INPDB
   OR     : pdbtools.py bfactorrmsd REFPDB INPDB [--out OUTPDB | --suffix=X]

INPDB  : input pdb file(s)
OUTPDB : output pdb file(s)
X         : suffix to add to input pdb for output
STARTRES  : starting residue number
STOPRES   : stopping residue number
RESNUM    : residue number
SEQ       : one-letter sequence codes
PHI,PSI   : desired angles after rotation, in degrees
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) == 1:
  print Usage
  sys.exit()
  

import re, string, urllib, os, tempfile, copy, random, gzip
import sequence, protein, scripttools


#global variables
Caps = ["ACE", "NHE", "NME"]
TerminiAtoms = ["OXT", "H2", "H3"]
OverlapDist = 0.65    #cutoff for detecting atomic overlap, in angstrom


def IsFileName(PdbString):
  """Returns True if PdbString is a filename, False if it is the data."""
  return not "\n" in PdbString

def IsData(PdbString):
  """Returns True if PdbString is the Pdb data, False if it is the filename."""
  return "\n" in PdbString

def ReturnPdbData(PdbString):
  """PdbString can be a pdb filename or the actual file contents."""
  if IsData(PdbString):
    return PdbString
  else:
    if os.path.isfile(PdbString):
      #check for gzip
      if PdbString.endswith(".gz"):
        return gzip.GzipFile(PdbString,"r").read()
      else:
        return file(PdbString,"rU").read()
    else:
      raise IOError, "Pdb file %s not found." % PdbString

def SavePdbData(PdbString, PdbFile):
  "Saves a pdb string to a file."
  if PdbFile is not None:
    file(PdbFile, "w").write(PdbString)
    

#-------- INFORMATION --------

def Atoms(Pdb):
  "Returns an array of atom names."
  Pdb = ReturnPdbData(Pdb)
  return [line[12:16] for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomRes(Pdb):
  "Returns an array of residue names for each atom."
  Pdb = ReturnPdbData(Pdb)
  return [line[22:26] for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomNum(Pdb):
  "Returns an array of atom numbers for each atom."
  Pdb = ReturnPdbData(Pdb)
  return [int(line[6:11]) for line in Pdb.split("\n") if line.startswith("ATOM")]

def AtomPos(Pdb):
  "Returns an array of atom coordinates."
  Pdb = ReturnPdbData(Pdb)
  return [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
          for l in Pdb.split("\n") if l.startswith("ATOM")]

def AtomResNums(Pdb):
  "Returns the numbers of the residues for each atom the pdb file data."
  Pdb = ReturnPdbData(Pdb)
  r = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      r.append(ResNum)
  return r 

def Seq(Pdb):
  "Returns an array of residue names for the sequence."
  Pdb = ReturnPdbData(Pdb)
  s = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      s.append(l[17:20])
      ThisResNum = ResNum
  return s

def ResLen(Pdb):
  "Returns the number of residues."
  return len(Seq(Pdb))

def AtomLen(Pdb):
  "Returns the number of atoms."
  return len(Atoms(Pdb))

def ResNums(Pdb):
  "Returns the numbers of the residues in the pdb file data."
  Pdb = ReturnPdbData(Pdb)
  r = []
  ThisResNum = ""
  for l in Pdb.split("\n"):
    if (not l.startswith("ATOM")):
      continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      r.append(ResNum)
  return r  

def GetCoords(Pdb):
  "Returns the Pdb coordinates."
  Pdb = ReturnPdbData(Pdb)
  Coords = [[float(s[30:38]), float(s[38:46]), float(s[46:54])]
            for s in Pdb.split("\n") if s.startswith("ATOM")]
  return Coords

def GetOverlaps(Pdb, MinDist = None, FirstOnly = False):
  """Returns a list of tuples (a,b) of overlaps between atoms a and b.
a and b start at zero and increment consecutively and are not taken
from the pdb data."""
  if MinDist is None: MinDist = OverlapDist
  #load atom positions
  Coords = GetCoords(Pdb)
  N = len(Coords)
  Overlaps = []
  MinDistSq = MinDist*MinDist
  #check atom distances
  for i in range(0,N):
    for j in range(i+1,N):
      DistSq = sum([(Coords[i][k]-Coords[j][k])**2 for k in range(0,3)])
      if DistSq < MinDistSq:
        Overlaps.append((i,j))
        if FirstOnly: return Overlaps
  return Overlaps

def HasOverlap(Pdb, MinDist = OverlapDist):
  return len(GetOverlaps(Pdb, MinDist = MinDist, FirstOnly = True)) > 0

def HasHydrogens(Pdb):
  "Indicates whether or not there are hydrogens."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if s[0:4]=="ATOM" and s[13]=="H"]
  return len(OutPdb) > 0

def HasCaps(Pdb):
  "Indicates if there are cap residues in a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if s[0:4]=="ATOM" and s[17:20] in Caps]
  return len(OutPdb) > 0


#--------PDB MODIFIERS--------

def Renumber(Pdb, StartAtom = 1, StartRes = 1, StartChain = " ",
             OutPdbFile = None):
  """Renumbers pdb residues and atoms starting at StartRes and StartAtoms,
  and starting the chain at no specification."""
  ChainLetters = [" "] + list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  r = StartRes - 1
  a = StartAtom - 1
  c = ChainLetters.index(StartChain)
  ThisResNum = ""
  ThisChain = ""
  #first get all the chains
  Chains = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      if not l[21] in Chains:
        Chains.append(l[21])
  Chains.sort()
  ChainDict = {}
  for i in range(0, len(Chains)):
    ind = (c + i) % len(ChainLetters)
    ChainDict[Chains[i]] = ChainLetters[ind]
  #now renumber
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      a += 1
      ResNum = l[22:29]
      if ThisResNum != ResNum:
        r += 1
        ThisResNum = ResNum
      ChainNum = ChainDict[l[21]]
      s = l[:6] + str(a).rjust(5) + l[11:21] + ChainNum \
          + str(r).rjust(4) + "   " + l[29:]
      OutPdb.append(s)
    else:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb 

def Center(Pdb, OutPdbFile = None):
  "Centers a pdb file."
  p = protein.ProteinClass(Pdb = Pdb)
  p.Center()
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb  

def Dehydrogen(Pdb, OutPdbFile = None):
  "Removes hydrogens from a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = [s for s in Pdb.split("\n") if not (s[0:4]=="ATOM" and s[13]=="H")]
  OutPdb = "\n".join(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Decap(Pdb, OutPdbFile = None):
  "Removes caps and terminal atoms from a pdb file."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      if l[12:16].strip() == "H1":
        OutPdb.append(l[:12] + " H  " + l[16:])
      elif not l[12:16].strip() in TerminiAtoms and not l[17:20] in Caps:
        OutPdb.append(l)
    else:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Cap(Pdb, OutPdbFile = None):
  "Adds caps ACE and NME to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapNRes = protein.ProteinClass(Seq = ">")
  CapCRes = protein.ProteinClass(Seq = "<")
  p = CapNRes + p + CapCRes
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def CapN(Pdb, OutPdbFile = None):
  "Adds n-terminal ACE cap to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapNRes = protein.ProteinClass(Seq = ">")
  p = CapNRes + p
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def CapC(Pdb, OutPdbFile = None):
  "Adds c-terminal NME cap to a pdb file."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  CapCRes = protein.ProteinClass(Seq = "<")
  p = p + CapCRes
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Splice(Pdb1, Pdb2, OutPdbFile = None):
  "Joins two pdb files together."
  #this uses proteinclass
  Pdb1 = ReturnPdbData(Pdb1)
  Pdb2 = ReturnPdbData(Pdb2)
  p1 = protein.ProteinClass(Pdb = Pdb1)
  p2 = protein.ProteinClass(Pdb = Pdb2)
  p3 = p1 + p2
  OutPdb = p3.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Generate(Seq, OutPdbFile = None):
  "Generates a pdb file of extended sequence Seq."
  #uses proteinclass
  p = protein.ProteinClass(Seq = Seq)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Extend(Pdb, SeqBefore, SeqAfter, OutPdbFile = None, Opt = True, N = 5):
  "Extends a pdbfile in either or both directions with given sequences."
  #uses proteinclass
  OutPdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  n = len(p)
  if len(SeqBefore) > 0:
    p = protein.ProteinClass(Seq = SeqBefore) + p
    if Opt and n > 0: p.OptimizePep(len(SeqBefore), N)
  if len(SeqAfter) > 0:
    p = p + protein.ProteinClass(Seq = SeqAfter)
    if Opt and n > 0: p.OptimizePep(n, N)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Extract(Pdb, StartRes, StopRes, OutPdbFile = None):
  "Extracts a portion of a pdb file between specified residue nums."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  ThisResNum = ""
  ResIndex = 0
  for l in Pdb.split("\n"):
    if not l.startswith("ATOM"): continue
    ResNum = l[22:29]
    if ThisResNum != ResNum:
      ThisResNum = ResNum
      ResIndex += 1
    if ResIndex >= StartRes and ResIndex <= StopRes:
      OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb
    
def Trim(Pdb, Chains = [], ResNums = [], OutPdbFile = None):
  """Selects out specific chains or residue numbers from a pdb file.
Here, residue nums are those in the file, and not necessarily in
sequential order starting at 1."""
  Pdb = ReturnPdbData(Pdb)
  pdblines = Pdb.split("\n")
  OutPdb = []
  for line in pdblines:
    if (not line.startswith("ATOM")):
      if line.startswith("TER"):
          line = line[:17] + lastres + line[20:23] + lastresnum + line[26:]
      OutPdb.append(line)
      continue
    chain = line[21]
    resnum = int(line[23:26])
    if (len(Chains) and chain not in Chains):
      continue
    if (len(ResNums) and resnum not in ResNums):
      continue
    lastres = line[17:20]
    lastresnum = line[23:26]
    OutPdb.append(line)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  SavePdbData(OutPdb, OutPdbFile)    
  return OutPdb

def Rotate(Pdb, ResNum, Phi = None, Psi = None, OutPdbFile = None):
  "Rotates a residue to specified phi, psi."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  if not ResNum in range(0, len(p)):
    raise IndexError, "Residue number %d not found" % ResNum
  try:
    p.RotateToPhiPsi(ResNum, Phi, Psi)
  except StandardError:
    print "Could not perform rotation."
    return Pdb
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def SpliceOpt(Pdb1, Pdb2, OutPdbFile = None, N = 5):
  "Splices two pdbs together and rotates to optimize for non-overlap."
  #this uses proteinclass
  Pdb1 = ReturnPdbData(Pdb1)
  Pdb2 = ReturnPdbData(Pdb2)
  p1 = protein.ProteinClass(Pdb = Pdb1)
  p2 = protein.ProteinClass(Pdb = Pdb2)
  p3 = p1 + p2
  if len(p1) > 0 and len(p2) > 0:
    p3.OptimizePep(len(p1), N)
  OutPdb = p3.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def RandDihedrals(Pdb, OutPdbFile = None, DeltaAng = 5.):
  "Adds a random angle to the torsions along the backbone."
  #this uses proteinclass
  Pdb = ReturnPdbData(Pdb)
  p = protein.ProteinClass(Pdb = Pdb)
  for i in range(0,len(p)):
    Phi, Psi = p.PhiPsi(i)
    if not Phi is None and not Psi is None:
      Phi += DeltaAng * (2.*random.random() - 1.)
      Psi += DeltaAng * (2.*random.random() - 1.)
      p.RotateToPhiPsi(i, Phi, Psi)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Standardize(Pdb, OutPdbFile = None):
  "Formats to standard residue and atom names, etc."
  Pdb = ReturnPdbData(Pdb)
  OutPdb = []
  for l in Pdb.split("\n"):
    if l.startswith("ATOM"):
      #fix histidines
      if l[17:20] in ["HID","HIE","HIP"]:
        l = l[:17] + "HIS" + l[20:]
      elif l[17:20] in ["CYX", "CYM"]:
        l = l[:17] + "CYS" + l[20:]
      #fix amide hydrogens
      if l[12:16] == " HN ":
        l = l[:12] + " H  " + l[16:]
    elif l.startswith("HETATM") and l[17:20] == 'MSE':
      #replace selenomethionine
      l = "ATOM  " + l[6:17] + "MET" + l[20:]
    #check alternate location indicator
    if l.startswith("ATOM"):
      if l[16] == "A":
        l = l[:16] + " " + l[17:]
      elif not l[16] == " ":
        continue
    OutPdb.append(l)
  OutPdb = "\n".join(OutPdb)
  OutPdb = Renumber(OutPdb)
  #now generate chains
  p = protein.ProteinClass(Pdb = OutPdb)
  p.GenerateChains()
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb

def Minimize(Pdb, OutPdbFile = None, TempDir = None,
             NSteps1 = 200, NSteps2 = 50):
  "Uses AMBER to minimize the energy, etc."
  import mdsim, shutil
  if TempDir is None:
    TempDir = tempfile.mkdtemp()
  s = mdsim.SimClass(TempDir)
  Pdb = ReturnPdbData(Pdb)
  Pdb = Dehydrogen(Pdb)
  fn = os.path.join(TempDir, "init.pdb")
  file(fn, "w").write(Pdb)
  s.SysInitPdb(fn)
  s.SysBuild()
  s.RunMin(NSteps1 = NSteps1, NSteps2 = NSteps2)
  s.RunEnergy()
  E = s["EPOT2"]
  OutPdb = file(os.path.join(TempDir, "current.pdb")).read()
  shutil.rmtree(TempDir)
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb, E

def GetEnergy(Pdb, TempDir = None):
  "Uses AMBER to get the potential energy."
  import mdsim, shutil
  if TempDir is None:
    TempDir = tempfile.mkdtemp()
  s = mdsim.SimClass(TempDir)
  Pdb = ReturnPdbData(Pdb)
  Pdb = Dehydrogen(Pdb)
  fn = os.path.join(TempDir, "init.pdb")
  file(fn, "w").write(Pdb)
  s.SysInitPdb(fn)
  s.SysBuild()
  s.RunEnergy()
  E = s["EPOT2"]
  shutil.rmtree(TempDir)
  return E

def DSSP(PdbFile):
  "Uses DSSP to get secondary structures."
  import subprocess
  p = subprocess.Popen("dsspcmbi %s" % PdbFile, shell = True,
                       stdout = subprocess.PIPE,
                       stderr = subprocess.STDOUT)
  ret = p.stdout.read()
  p.stdout.close()
  p.wait()
  #check that it ran
  if not "Secondary Structure Definition" in ret:
    raise ValueError, "Could not run dsspcmbi."
  ret = ret[ret.index("#  RESIDUE AA"):].strip()
  lines = ret.split("\n")
  ss = "".join([l[16] for l in lines[1:]])
  ss = ss.replace(" ", "-")
  ss = ss.replace("T", "-")
  ss = ss.replace("S", "-")
  ss = ss.replace("B", "E")
  return ss

def BFactorRMSD(RefPdb, Pdb, OutPdbFile = None):
  "Adds BFactors proportional to atom RMSD."
  import rmsd
  Pdb = ReturnPdbData(Pdb)
  RefPdb = ReturnPdbData(RefPdb)
  p = protein.ProteinClass(Pdb = Pdb)
  pref = protein.ProteinClass(Pdb = RefPdb)  
  r, NRes = rmsd.RMSDProteinClass(pref, p, Backbone = False, AlignSeq = True,
                                  UpdateBFactors = True, AlignAtoms = True)
  OutPdb = p.GetPdb()
  SavePdbData(OutPdb, OutPdbFile)
  return OutPdb
     

#--------PDB DOWNLOAD--------

def Download(Id, OutPdbFile = None):
  "Downloads a pdb file."
  pdbmirror = "ftp://ftp.rcsb.org/pub/pdb"
  if (not re.match("^[1-9][1-9a-z]{3}", Id)):
    print "warning: " + Id + " not a valid PDB code"
    return ""
  Id = string.lower(Id)
  url = pdbmirror + "/data/structures/all/pdb/pdb" + Id + ".ent.Z"
  urlhandle = urllib.urlopen(url)
  pdbz = urlhandle.read()
  urlhandle.close()
  OutPdb = __Uncompress(pdbz)

  SavePdbData(OutPdb, OutPdbFile)  
  return OutPdb

def __Uncompress(Data):
  "Accessory for uncompressing downloaded pdb."
  #uin, uout = os.popen2("uncompress")
  #
  #uin.write(data)
  #uin.close()
  
  compfd, comppath = tempfile.mkstemp(".Z")
  compfile = os.fdopen(compfd, "w")
  compfile.write(Data)
  compfile.close()
  
  uout = os.popen("uncompress -c " + comppath)
  
  uncompdata = uout.read()
  uout.close()
  
  return uncompdata




def ParseArgs(Args, NArg = 0):
  """Parses arguments to get input and output pdb files."""
  #shift to remove command itself
  NArg += 2
  if "suffix" in Args:
    Suff = Args["suffix"]
    return [(x, x.replace(".pdb", "%s.pdb" % Suff)) for x in Args["ARGS"][NArg:]]
  elif "out" in Args["FLAGS"]:
    n = Args["FLAGPOS"][Args["FLAGS"].index("out")] - NArg
    InPdbs = Args["ARGS"][NArg:NArg+n]
    OutPdbs = Args["ARGS"][NArg+n:NArg+2*n]
    if len(InPdbs) != len(OutPdbs):
      print "Number of input and output pdb files does not match!"
      sys.exit()
    return zip(InPdbs, OutPdbs)
  else:
    return [(x, x) for x in Args["ARGS"][NArg:]]


#command-line running
if __name__ == "__main__":
  Args = scripttools.ParseArgs(sys.argv)
  Cmd = Args[1].lower()
  
  if Cmd == "cap":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Cap(InPdb, OutPdbFile = OutPdb)
      print "Capped %s as %s." % (InPdb, OutPdb)

  elif Cmd == "decap":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Decap(InPdb, OutPdbFile = OutPdb)
      print "Decapped %s as %s." % (InPdb, OutPdb)

  elif Cmd == "dehydrogen":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Dehydrogen(InPdb, OutPdbFile = OutPdb)
      print "Dehydrogenated %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extendc":
    Seq = Args[2]
    Pdbs = ParseArgs(Args, 1)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Extend(InPdb, [], Seq, OutPdbFile = OutPdb, Opt = True)
      print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extendn":
    Seq = Args[2]
    Pdbs = ParseArgs(Args, 1)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Extend(InPdb, Seq, [], OutPdbFile = OutPdb, Opt = True)
      print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "extend":
    SeqBefore, SeqAfter = Args[2:4]
    Pdbs = ParseArgs(Args, 2)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Extend(InPdb, SeqBefore = SeqBefore, SeqAfter = SeqAfter,
                   OutPdbFile = OutPdb, Opt = True)
      print "Extended and optimized for overlaps %s as %s." % (InPdb, OutPdb)

  elif Cmd == "generate":
    Seq, OutPdb = Args[2:4]
    pdb = Generate(Seq, OutPdbFile = OutPdb)
    print "Generated as %s." % OutPdb

  elif Cmd == "renumber":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Renumber(InPdb, OutPdbFile = OutPdb)
      print "Renumbered %s as %s." % (InPdb, OutPdb)

  elif Cmd == "center":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Center(InPdb, OutPdbFile = OutPdb)
      print "Centered %s as %s." % (InPdb, OutPdb)

  elif Cmd == "rotate":
    ResNum = int(Args[2]) - 1
    Pdbs = ParseArgs(Args, 3)
    Phi, Psi = [float(x) for x in Args[3:5]]
    for (InPdb, OutPdb) in Pdbs:
      pdb = Rotate(InPdb, ResNum, Phi, Psi, OutPdbFile = OutPdb)
      print "Rotated %s as %s." % (InPdb, OutPdb)    

  elif Cmd == "showseq":
    InPdb = Args[2]
    Seq = sequence.SeqToList(Seq(InPdb))
    Nums = range(1, len(Seq)+1)
    n = 15
    i = 0
    print "\n" + sequence.SeqToAA1(Seq) + "\n"
    while i < len(Seq):
      print " ".join(["%-3d" % x for x in Nums[i:i+n]])
      print " ".join(["%-3s" % x for x in Seq[i:i+n]])
      print ""
      i += n

  elif Cmd == "showoverlaps":
    InPdb = Args[2]
    if len(Args["ARGS"]) > 3:
      OverlapDist = float(Args[3])
    aname = [x.strip() for x in Atoms(InPdb)]
    anum = AtomNum(InPdb)
    Overlaps = GetOverlaps(InPdb)
    if len(Overlaps) == 0:
      print "No overlaps detected."
    else:
      print "Overlaps (less than %.2f A):" % (OverlapDist)
      for (a,b) in Overlaps:
        print "  %s%d-%s%d" % (aname[a], anum[a], aname[b], anum[b])

  elif Cmd == "splice":
    InPdb1, InPdb2, OutPdb = Args[2:5]
    pdb = Splice(InPdb1, InPdb2, OutPdbFile = OutPdb)
    print "Spliced %s and %s as %s." % (InPdb1, InPdb2, OutPdb)

  elif Cmd == "spliceopt":
    InPdb1, InPdb2, OutPdb = Args[2:5]
    pdb = SpliceOpt(InPdb1, InPdb2, OutPdbFile = OutPdb)
    print "Spliced and optimized for overlaps %s and %s as %s." \
          % (InPdb1, InPdb2, OutPdb)

  elif Cmd == "trim":
    StartRes = int(Args[2])
    StopRes = int(Args[3])
    Pdbs = ParseArgs(Args, 2)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Extract(InPdb, StartRes, StopRes, OutPdbFile = OutPdb)
      print "Extracted residues from %s as %s." % (InPdb, OutPdb)

  elif Cmd == "standardize":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      pdb = Standardize(InPdb, OutPdbFile = OutPdb)
      print "Standardized pdb file from %s to %s" % (InPdb, OutPdb)

  elif Cmd == "alignseq":
    InPdb1, InPdb2 = Args[2:]
    Seq1, Seq2 = Seq(InPdb1), Seq(InPdb2)
    Map = sequence.SeqMapClass(Seq1, Seq2)
    if len(Map) <= 1:
      print "Pdb files do not align."
    else:
      print "%s: residues %d to %d (%d of %d)" % (InPdb1, Map.a + 1, Map.b,
                                                  len(Map), len(Seq1))
      print "%s: residues %d to %d (%d of %d)" % (InPdb2, Map.c + 1, Map.d,
                                                  len(Map), len(Seq2))
      Seq = Seq1[Map.a:Map.b]
      Nums1, Nums2 = range(Map.a+1, Map.b+1), range(Map.c+1, Map.d+1)
      n = 15
      i = 0
      print "\nOverlap sequence:"
      print sequence.SeqToAA1(Seq) + "\n"
      print sequence.SeqToAA1(Seq2[Map.c:Map.d]) + "\n"
      while i < len(Seq):
        print " ".join(["%-3d" % x for x in Nums1[i:i+n]])
        print " ".join(["%-3d" % x for x in Nums2[i:i+n]])
        print " ".join(["%-3s" % x for x in Seq[i:i+n]])
        print ""
        i += n

  elif Cmd == "minimize":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      print "Minimizing %s as %s" % (InPdb, OutPdb)
      NSteps1, NSteps2 = None, None
      if "steps1" in Args: NSteps1 = int(Args["steps1"])
      if "steps2" in Args: NSteps2 = int(Args["steps2"])
      pdb, E = Minimize(InPdb, OutPdbFile = OutPdb, TempDir = "_pdbmintmp",
                        NSteps1 = NSteps1, NSteps2 = NSteps2)
      print "Minimized energy is: %.2f\n" % E

  elif Cmd == "longminimize":
    Pdbs = ParseArgs(Args, 0)
    for (InPdb, OutPdb) in Pdbs:
      print "Minimizing %s as %s" % (InPdb, OutPdb)
      pdb, E = Minimize(InPdb, OutPdbFile = OutPdb, TempDir = "_pdbmintmp",
                        NSteps1 = 200, NSteps2 = 2000)
      print "Minimized energy is: %.2f\n" % E

  elif Cmd == "getenergy":
    Pdbs = ParseArgs(Args, 0)
    l = max([len(x) for (x,y) in Pdbs])
    print "%-*s  %s" % (l, "PDB", "Potential_Energy")
    for (InPdb, OutPdb) in Pdbs:
      E = GetEnergy(InPdb, TempDir = "_pdbmintmp")
      print "%-*s  %.3f" % (l, InPdb, E)

  elif Cmd == "dssp":
    Pdbs = ParseArgs(Args, 0)
    l = max([len(x) for (x,y) in Pdbs])
    print "%-*s %-5s %-5s %s" % (l, "PDB", "%H", "%E", "Secondary_structure")
    for (InPdb, OutPdb) in Pdbs:
      SS = DSSP(InPdb)
      PH = float(SS.count("H")) / len(SS) * 100.
      PE = float(SS.count("E")) / len(SS) * 100.
      print "%-*s %-5.1f %-5.1f %s" % (l, InPdb, PH, PE, SS)

  elif Cmd == "bfactorrmsd":
    RefPdb = Args[2]
    Pdbs = ParseArgs(Args, 1)
    for (InPdb, OutPdb) in Pdbs:
      pdb = BFactorRMSD(RefPdb, InPdb, OutPdbFile = OutPdb)
      print "Added atom RMSDs to bfactors for %s as %s" % (InPdb, OutPdb)
      
  else:
    print "Command not recognized."
    