#!/usr/bin/env python

#LAST MODIFIED: 10-29-08

#DESCRIPTION: Provides reading and writing routines for popular
#coordinate file formats.  Currently: pdb and amber trajectory.

from numpy import *
import copy, os, gzip, pdbtools

#Masks for backbone atoms
NoMask = []
PdbBackboneMask = ["CA", "C", "N"]
CrdBackboneMask = ["CA", "C", "N"]
BackboneMask = ["CA", "C", "N"]
AlphaCarbonMask = ["CA"]
Caps = ["NHE","NME","ACE"]


def GetPdbSeq(PdbFile):
  "Gets the sequence from a PdbFile."
  LineList = [l for l in open(PdbFile, "rU").readlines() if l[0:4] == 'ATOM']
  Seq = []
  ThisResNum = 0
  for l in LineList:
    ResNum = int(l[22:29])
    if ThisResNum != ResNum:
      Seq.append(l[17:20].upper())
      ThisResNum = ResNum
  return Seq    

def GetPdbCoords(PdbFile, Mask = NoMask):
  "Gets coordinates from a Pdb file."
  if os.path.isfile(PdbFile):
    f = open(PdbFile, "rU")
    Pos = []
    while True:
      s = f.readline()
      if len(s) == 0: break
      if s[0:4] == "ATOM" and (s[13:15].strip() in Mask or Mask == NoMask):
        Pos.append([float(s[30:38]),float(s[38:46]), float(s[46:54])])
    f.close()
    return array(Pos,float)
  else:
    raise IOError, "Pdb file does not exist."

def SavePdbCoordsPdb(Pos, PdbFile, PdbTmplFile, Mask = NoMask):
  "Saves a coordinate array to a pdb file using a template pdb file."
  if os.path.isfile(PdbTmplFile):
    ft = open(PdbTmplFile, "rU")
    f = open(PdbFile, "w")
    AtomNum = 0
    while True:
      s = ft.readline()
      if len(s) == 0: break
      if s.startswith("TER") or s.startswith("MODEL") or s.startswith("ENDMDL"):
        f.write(s)
      elif s.startswith("ATOM") and (s[13:15].strip() in Mask or Mask == NoMask):
        s = s[0:30] + "%8.3f%8.3f%8.3f" % (Pos[AtomNum,0],
            Pos[AtomNum,1], Pos[AtomNum,2]) + s[54:]
        AtomNum += 1
        f.write(s)
    ft.close()
    f.close()
    if not AtomNum == len(Pos):
      print "Warning: pdb template file %s has %d coordinates, but I have %d." \
            % (PdbTmplFile, AtomNum, len(Pos))
  else:
    raise IOError, "Warning: could not find template pdb file."

def SavePdbCoords(Pos, AtomNames, AtomRes, Seq, PdbFile,
                  Mask = NoMask, Standardize = True):
  "Saves a coordinate array to a pdb file."
  Seq = copy.deepcopy(Seq)
  if Standardize:
    for i in range(0, len(Seq)):
      if Seq[i] in ["HIE","HID","HIP"]:
        Seq[i] = "HIS"
      if Seq[i] in ["CYX"]:
        Seq[i] = "CYS"
  Pdb = "MODEL\n"
  ind = 1
  for (i, (an, arn)) in enumerate(zip(AtomNames, AtomRes)):
    #check mask
    if not (Mask == NoMask or Mask is None):
      if not an.strip() in Mask: continue
    (x, y, z) = Pos[i,:]
    ar = Seq[arn]
    Pdb += "ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f \n" \
      % (ind, an, ar, arn + 1, x, y, z)
    ind += 1
  Pdb += "TER\nENDMDL\n"
  file(PdbFile, "w").write(Pdb)

def SavePdbCoordsObj(Pos, CoordsObj, PdbFile, Mask = NoMask, Standardize = True):
  "Saves a coordinate array to a pdb file using a coords obj."
  SavePdbCoords(Pos, CoordsObj.AtomNames, CoordsObj.AtomRes,
                CoordsObj.Seq, PdbFile, Mask = Mask, Standardize = Standardize)

def ParseCrdString(s, AtomNames = [], Mask = NoMask):
  """Takes a text string of Crd coordinates and parses into an array."""
  #parse into a n by 3 array
  try:
    vals = s.replace("\n","")
    vals = [float(vals[i:i+8]) for i in range(0, len(vals), 8)]
    Pos = array(vals, float)
  except ValueError:
    raise ValueError, "Improper number of coordinates found in Crd string."
  if mod(len(Pos), 3) == 0:
    Pos = reshape(Pos, (-1,3))
    #if using mask remove the extraneous coordinates
    if not Mask == NoMask and len(AtomNames) == len(Pos):
      Pos = compress([a.strip() in Mask for a in AtomNames], Pos, 0)
    return Pos
  else:
    raise ValueError, "Improper number of coordinates found in Crd string."


def ParseRstString(s, AtomNames = [], Mask = NoMask):
  """Takes a text string of Rst coordinates and parses into an array;
Note that this will return velocities as well."""
  #parse into a n by 3 array
  try:
    vals = s.replace("\n","")
    vals = [float(vals[i:i+12]) for i in range(0, len(vals), 12)]
    Pos = array(vals, float)
  except ValueError:
    raise ValueError, "Improper number of coordinates found in Crd string."
  if mod(len(Pos), 3) == 0:
    Pos = reshape(Pos, (-1,3))
    #if using mask remove the extraneous coordinates
    if not Mask == NoMask and len(AtomNames) == len(Pos):
      Pos = compress([a.strip() in Mask for a in AtomNames], Pos, 0)
    return Pos
  else:
    raise IOError, "Improper number of coordinates found in Rst file."


def GetPrmtopBlock(PrmtopFile, Flag, BlockLen = None):
  "Gets Prmtop data for a specified Flag."
  if os.path.isfile(PrmtopFile):
    f = open(PrmtopFile, "rU")
    Dat = []
    line = f.readline()
    while len(line) > 0 and not line.startswith(Flag):
      line = f.readline()
    f.readline()
    line = f.readline()
    while len(line) > 0 and not line[0:1] == "%":
      if BlockLen is None:
        Dat.extend(line.split())
      else:
        Dat.extend([line[i:i+BlockLen] for i in xrange(0,len(line)-1,BlockLen)])
      line = f.readline()
    f.close()
    return Dat
  else:
    raise IOError, "Could not find Prmtop file."
    return None

def GetPrmtopAtomNames(PrmtopFile):
  "Gets the names of atoms from a Prmtop file."
  Names = GetPrmtopBlock(PrmtopFile, "%FLAG ATOM_NAME", BlockLen = 4)
  return Names

def AmbToPdbAtomNames(Names):
  "Converts Amber-style atom names to Pdb-style."
  return [an[-1] + an[:-1] for an in Names]

def GetPrmtopSeq(PrmtopFile):
  "Gets the names of residues from a Prmtop file."
  Seq = GetPrmtopBlock(PrmtopFile, "%FLAG RESIDUE_LABEL", BlockLen = 4)
  Seq = [s.strip() for s in Seq]
  return Seq

def GetPrmtopAtomRes(PrmtopFile):
  "Gets the atom's residue numbers from a Prmtop file."
  ResPtr = GetPrmtopBlock(PrmtopFile, "%FLAG RESIDUE_POINTER", BlockLen = 8)
  if ResPtr is None: return None
  AtomNames = GetPrmtopAtomNames(PrmtopFile)
  ResPtr = [int(x) for x in ResPtr] + [len(AtomNames) + 1]
  AtomRes = []
  for i in range(0, len(ResPtr)-1):
    AtomRes.extend([i] * (ResPtr[i+1] - ResPtr[i]))
  return AtomRes


def GetCrdCoords(CrdFile, PrmtopFile = "", Mask = NoMask):
  "Gets coordinates from a Crd and Prmtop set of files."
  if len(Mask) > 0:
    AtomNames = GetPrmtopAtomNames(PrmtopFile)
    if AtomNames is None: return
  else:
    AtomNames = []
  if os.path.isfile(CrdFile):
    f = open(CrdFile, "rU")
    f.readline()  #bypass the header lines
    f.readline()
    s = f.read()
    f.close()
    return ParseCrdString(s, AtomNames, Mask)
  else:
    raise IOError, "Crd file not found."


def SaveCrdCoords(Pos, CrdFile, Mode = "wb", Head = True):
  "Saves coordinates to a Crd file."
  NPerLine = 10
  Fmt = "%8.3f"
  f = open(CrdFile, Mode)
  if Head: f.write("ACE".ljust(80) + "\n")
  p = [Fmt % x for x in Pos.flat]
  r = mod(len(p),NPerLine)
  if r > 0: p.extend([""]*(NPerLine-r))
  for i in xrange(0, len(p), NPerLine):
    f.write("".join(p[i:i+NPerLine]) + "\n")
  f.close()


def GetRstCoords(RstFile, PrmtopFile = "", Mask = NoMask):
  "Gets coordinates from an amber restart file."
  if len(Mask) > 0:
    AtomNames = GetPrmtopAtomNames(PrmtopFile)
    if AtomNames is None: return
  else:
    AtomNames = []
  if os.path.isfile(RstFile):
    f = open(RstFile, "rU")
    f.readline()  #bypass the header lines
    f.readline()
    s = f.read()
    f.close()
    return ParseRstString(s, AtomNames, Mask)
  else:
    raise IOError, "Rst file not found."


def SaveRstCoords(Pos, RstFile):
  "Saves coordinates to an amber restart file."
  NPerLine = 6
  Fmt = "%12.7f"
  f = open(RstFile, "w")
  f.write("ACE".ljust(80) + "\n")
  f.write("%5d  0.0000000E+00\n" % len(Pos))
  p = [Fmt % x for x in Pos.flat]
  r = mod(len(p),NPerLine)
  if r > 0: p.extend([""]*(NPerLine-r))
  for i in xrange(0, len(p), NPerLine):
    f.write("".join(p[i:i+NPerLine]) + "\n")
  f.close()  
  

def AmbToPdb(RstFile, PrmtopFile, PdbFile, AAtm = False, BRes = False):
  """Covert parmtop and coordinate files into a PDB string
* RstFile: file path to amber restart file file
* PrmtopFile: file path to parmtop file
* PdbFile: pdb file to write
* AAtm: boolean, use Amber atom names in the PDB file
* BRes: boolean, use PDB-standard residue names"""
  ambcmd = "ambpdb -p " + PrmtopFile
  if (AAtm):
    ambcmd += " -aatm"
  if (BRes):
    ambcmd += " -bres"
  ambpdb = os.popen(ambcmd + " < " + RstFile + " 2> /dev/null", "r")
  pdb = ambpdb.read()
  ambpdb.close()
  file(PdbFile,"w").write(pdb)


def SavePdbCoordsAmb(Pos, PrmtopFile, PdbFile, Standardize = True):
  """Saves a coordinate array to a pdb file using an amber prmtop file."""
  AtomNames = GetPrmtopAtomNames(PrmtopFile)
  AtomNames = AmbToPdbAtomNames(AtomNames)
  AtomRes = GetPrmtopAtomRes(PrmtopFile)
  Seq = GetPrmtopSeq(PrmtopFile)
  SavePdbCoords(Pos, AtomNames, AtomRes, Seq, PdbFile,
                Standardize = Standardize)

def SavePdbCoordsAmb2(Pos, PrmtopFile, PdbFile, AAtm = False, BRes = False):
  """Saves a coordinate array to a pdb file using an amber prmtop file."""
  SaveRstCoords(Pos, "rst.tmp")
  AmbToPdb("rst.tmp", PrmtopFile, PdbFile, AAtm, BRes)
  os.remove("rst.tmp")

def PdbsToCrd(PdbFileList, CrdFile):
  """Saves multiple pdb files as a crd file"""
  NPerLine = 10
  Fmt = "%8.3f"
  NAtom = -1
  f = open(CrdFile, "w")
  f.write("ACE".ljust(80) + "\n")
  for fn in PdbFileList:
    Pos = GetPdbCoords(fn)
    if NAtom < 0:
      NAtom = len(Pos)
    elif not len(Pos) == NAtom:
      f.close()
      raise IOError, "Pdb files do not have same number of atoms."      
    p = [Fmt % x for x in Pos.flat]
    r = mod(len(p),NPerLine)
    if r > 0: p.extend([""]*(NPerLine-r))
    for i in xrange(0, len(p), NPerLine):
      f.write("".join(p[i:i+NPerLine]) + "\n")
  f.close()
    
  

def GetTrjLenOld(TrjFile, PrmtopFile):
  "Gets the number of frames in a trajectory."
  #get the number of atoms using the prmtop file
  Names = GetPrmtopAtomNames(PrmtopFile)
  if Names is None: return 0
  NAtom = len(Names)
  #calculate the number of lines to read per coordinate set
  BlockLen = (3 * NAtom) / 10
  if (3 * NAtom) % 10 > 0: BlockLen += 1
  #use gzip?
  if TrjFile.split(".")[-1].strip().lower() == "gz":
    FileMthd = gzip.GzipFile
  else:
    FileMthd = file
  #open the file and count the lines
  NLine = -1
  f = FileMthd(TrjFile, "r")
  for lines in f: NLine += 1
  f.close()
  return NLine/BlockLen

def GetTrjLen(TrjFile, PrmtopFile):
  "Gets the number of frames in a trajectory."
  Names = GetPrmtopAtomNames(PrmtopFile)
  if Names is None: return 0
  NAtom = len(Names)
  #calculate the number of lines to read per coordinate set
  r = mod(3 * NAtom, 10)
  NLine = (3 * NAtom) / 10
  if r > 0: NLine += 1
  #use gzip?
  if TrjFile.strip().lower().endswith("gz"):
    FileMthd = gzip.GzipFile
  else:
    FileMthd = file
  #open the trj file
  try:
    Trj = FileMthd(TrjFile, "r")
  except IOError:
    raise IOError, "There was an error opening the trajectory file."
    return
  #read the header
  BytesHead = len(Trj.readline())
  #read the first config
  BytesCoords = 0
  for i in range(NLine):
    BytesCoords += len(Trj.readline())
  #unfortunately for gzip files, we can't seek to end, so we need to
  #do this bit to get the total size:
  Trj.seek(0)
  BytesTot = 0
  #try to catch errors in the gzipped file; stride 1024 bytes
  try:
    while len(Trj.read(1024)) > 0: pass
  except IOError:
    #find the exact byte number
    try:
      while len(Trj.read(1)) > 0: pass
    except IOError:
      pass
  BytesTot = max(Trj.tell(), 0)
  #return the total number of configs
  Trj.close()
  del Trj
  return int((BytesTot - BytesHead) / BytesCoords)


class TrjClass:
  "Provides a class for reading successive sets of coordinates from Amber trajectory files."
  
  def __init__(self, TrjFile, PrmtopFile, Mask = NoMask,
               NSkip = 0, NRead = None, NStride = 1,
               LinkPos = None):
    """Initializes the class and opens the trajectory file for reading.
* TrjFile: string name of trj file
* PrmtopFile: string name of prmtop file
* Mask: list of strings; filter for atom names (default is no mask/empty list)
* NSkip: number of configurations to skip
* NRead: maximum number of configurations to read (default is all)
* NStride: stride between configuration frames (default is 1)
* LinkPos: an outside array that is updated automatically as Pos are read"""
    IsFile1, IsFile2 = os.path.isfile(PrmtopFile), os.path.isfile(TrjFile)
    if IsFile1 and IsFile2:
      #set the filenames
      self.TrjFile = TrjFile
      self.PrmtopFile = PrmtopFile
      #set the frames to skip, read, and stride
      if NSkip < 0:
        raise ValueError, "NSkip is less than zero."
      self.NSkip = NSkip
      self.NRead = NRead
      self.NStride = max(NStride, 1)
      if self.NRead < 0: self.NRead = None
      #set the byte fields
      self.BytesHead = 0
      self.BytesCoords = 0
      self.BytesTot = 0
      self.NCoords = 0
      self.SliceNCoords = 0
      #set the counters
      self.Count = 0
      self.Index = -1
      self.SliceIndex = -1
      #set the mask option
      self.Mask = Mask
      if self.Mask is None: self.Mask = NoMask
      #set the linked pos
      self.LinkPos = LinkPos
      #initialize everything
      self.__Init()
    else:
      if not IsFile2:
        raise IOError, "Could not find %s." % TrjFile
      if not IsFile1:
        raise IOError, "Could not find %s." % PrmtopFile

  def Reset(self):
    "Resets current configuration to trajectory start."
    if not self.__Trj is None:
      #skip to the right configuration
      self.__Trj.seek(self.BytesHead + self.NSkip * self.BytesCoords)
    #reset the indices
    self.Index, self.SliceIndex = -1, -1
    #set the coordinate set number, corresponding to last set read in
    self.Count = 0

  def __Open(self):
    """Opens the trajectory file for reading"""
    if self.__Trj is None:
      try:
        self.__Trj = self.__FileMthd(self.TrjFile, "r")
      except IOError:
        raise IOError, "There was an error opening the trajectory file."
      
  def Close(self):
    "Closes any open files."
    if self.__Trj is not None:
      self.__Trj.close()
      self.__Trj = None

  def __Init(self):
    """Initializes internal variables from disk data."""
    #get the atom names, seq, and residue nums
    self.AtomNames = GetPrmtopAtomNames(self.PrmtopFile)
    self.AtomNames = AmbToPdbAtomNames(self.AtomNames)
    self.NAtom = len(self.AtomNames)
    if self.NAtom == 0:
      raise IOError, "There was an error reading the Prmtop file."
      return
    self.AtomRes = GetPrmtopAtomRes(self.PrmtopFile)
    self.Seq = GetPrmtopSeq(self.PrmtopFile)
    #set the file method
    if self.TrjFile.split(".")[-1].strip().lower() == "gz":
      self.__FileMthd = gzip.GzipFile
      self.Gzip = True
    else:
      self.__FileMthd = open
      self.Gzip = False
    self.__Trj = None
    #count and reset
    self.__Open()
    self.__CountBytes()
    self.Reset()
    #get an initial set of positions
    if len(self) > 0:
      self.Pos = self[0]
    else:
      self.Pos = None
    #close for now
    self.Close()

  def __CountBytes(self):
    """Counts the number of bytes per configuration by reading the first set."""
    #calculate the number of lines to read per coordinate set
    r = mod(3 * self.NAtom, 10)
    self.__NLine = (3 * self.NAtom) / 10
    if r > 0: self.__NLine += 1
    #close the trj file if it's open
    self.Close()
    #open the trj file
    self.__Open()
    #read the header
    self.BytesHead = len(self.__Trj.readline())
    #read the first config
    self.BytesCoords = 0
    for i in range(0, self.__NLine):
      self.BytesCoords += len(self.__Trj.readline())
    if self.Gzip:
      #unfortunately for gzip files, we need to do this bit to get the total size:
      self.__Trj.seek(0)
      self.BytesTot = 0
      Incr = os.path.getsize(self.TrjFile)
      while self.BytesTot == self.__Trj.tell():
        self.BytesTot += Incr
        self.__Trj.seek(self.BytesTot)
      self.BytesTot = max(self.__Trj.tell(), 0)
    else:
      self.BytesTot = os.path.getsize(self.TrjFile)
    if 1==2: 
      #try to catch errors in the gzipped file; stride 1024 bytes
      try:
        while len(self.__Trj.read(1024)) > 0: pass
      except IOError:
        #find the exact byte number
        try:
          while len(self.__Trj.read(1)) > 0: pass
        except IOError:
          pass
      self.BytesTot = max(self.__Trj.tell(), 0)
    #set the total number of configs
    self.NCoords = int((self.BytesTot - self.BytesHead) / self.BytesCoords)
    self.NSkip = min(self.NSkip, self.NCoords)
    a = max(0, self.BytesTot - self.BytesHead - self.BytesCoords*self.NSkip)
    b = self.NStride * self.BytesCoords
    self.SliceNCoords = int(a / b)  
    if a % b > 0: self.SliceNCoords += 1
    #set the limit of how many to read in
    if self.NRead is None:
      self.NRead = self.SliceNCoords
    else:
      self.SliceNCoords = min(self.NRead, self.SliceNCoords)

  def Get(self, ind, Mask = None):
    if Mask is None: Mask = self.Mask
    #this index is relative to the sliced version
    #check for reverse notation
    if ind < 0:
      ind += self.SliceNCoords
    #check bounds
    if ind < 0 or ind >= self.SliceNCoords:
      raise IndexError, "Index out of bounds for trj class."
    else:
      #make sure we're open
      self.__Open()
      #calculate the absolute index
      self.SliceIndex = ind
      self.Index = self.NSkip + self.NStride * ind
      #seek to the right position
      self.__Trj.seek(self.BytesHead + self.BytesCoords*self.Index)
      #read the data
      s = self.__Trj.read(self.BytesCoords)
      #check to see if we ran out of coordinates
      if len(s) < self.BytesCoords:
        raise IOError
        return
      #parse the crd string
      self.Pos = ParseCrdString(s, self.AtomNames, Mask)
      #update a linked coord array
      if not self.LinkPos is None: self.LinkPos[:,:] = self.Pos
      return self.Pos

  def __getitem__(self, ind):
    return self.Get(ind)

  def __len__(self):
    return self.SliceNCoords

  def __iter__(self):
    self.Reset()
    return self

  def next(self):
    ind = self.SliceIndex + 1
    if ind < self.SliceNCoords:
      self.Count += 1
      return self[ind]
    else:
      self.Reset()
      raise StopIteration

  def GetNextCoords(self, Mask = None):
    "Returns the next set of coordinates in the Trj file, or None if the end is reached."
    ind = self.SliceIndex + 1
    if ind < self.SliceNCoords:
      self.Count += 1
      result = self.Get(ind, Mask)
    else:
      result = None
    return result

  def GetIndices(self):
    "Returns the configuration indices of all read in so far."
    return range(self.NSkip, self.NCoords, self.NStride)


class PdbListClass:
  "Provides a class for reading successive sets of coordinates from pdb files."
  
  def __init__(self, PdbFileList, Mask = NoMask,
               LinkPos = None):
    """Initializes the class and checks for pdb file existence.
* PdbFileList: list of string names of pdb files
* Mask: list of strings; filter for atom names (default is no mask/empty list)
* LinkPos: an outside array that is updated automatically as coords are read
"""
    #check for file existence
    self.PdbFileList = [f for f in PdbFileList if os.path.isfile(f)]
    Diff = len(PdbFileList) - len(self.PdbFileList)
    if Diff > 0:
      raise IOError, "Could not find %d files; removing from list" % Diff
    #set the number of positions
    self.LastLen = -1
    #set the total count
    self.TotalCount = len(self.PdbFileList)
    #set the mask option
    self.Mask = Mask
    if self.Mask is None: self.Mask = NoMask
    #set the linked pos
    self.LinkPos = LinkPos
    #get the sequence, atom names, and atom residues
    f = self.PdbFileList[0]
    self.AtomNames = pdbtools.Atoms(f)
    self.AtomRes = [int(x)-1 for x in pdbtools.AtomResNums(pdbtools.Renumber(f))]
    self.Seq = pdbtools.Seq(f)
    #reset
    self.Reset()
    #get initial positions
    if len(self) > 0:
      self.Pos = self[0]
      self.Reset()
    else:
      self.Pos = None

  def __len__(self):
    "Returns the number of configurations."
    return len(self.PdbFileList)
    
  def Reset(self):
    "Resets current configuration to list start."
    #set the coordinate set number, corresponding to last set read in
    self.Count = 0
    self.LastLen = -1
    #set the index
    self.Index = -1

  def GetNextCoords(self, Mask = None):
    "Returns the next set of pdb coordinates, or None if list end is reached."
    ind = self.Index + 1
    if ind < len(self.PdbFileList):
      self.Count += 1
      result = self.Get(ind, Mask)
    else:
      result = None
    return result

  def GetIndices(self):
    "Returns the configuration indices of all read in so far."
    return range(len(self.PdbFileList))

  def Get(self, ind, Mask = None):
    if Mask is None: Mask = self.Mask
    if ind < 0: ind += len(self.PdbFileList)
    if ind < 0 or ind >= len(self.PdbFileList):
      raise IndexError, "Index out of bounds for pdb coords class."
    self.Index = ind
    self.Pos = GetPdbCoords(self.PdbFileList[ind], Mask)
    if self.LastLen > 0 and not self.LastLen == len(self.Pos):
      raise ValueError, "Configuration read with different number of atoms from last read."
    self.LastLen = len(self.Pos)
    #update the linked pos
    if not self.LinkPos is None: self.LinkPos[:,:] = self.Pos
    return self.Pos

  def __getitem__(self, ind):
    return self.Get(ind)

  def __iter__(self):
    self.Reset()
    return self

  def next(self):
    ind = self.Index + 1
    if ind < len(self.PdbFileList):
      self.Count += 1
      return self[ind]
    else:
      self.Reset()
      raise StopIteration

    
class MultiCoordClass:
  "Provides a class for linking multiple coordinate files together."
  
  def __init__(self, CoordObjList, LinkPos = None):
    """Initializes the class and opens the trajectory file for reading.
* CoordObjList gives a list of the coordinate objects
* LinkPos: an outside array that is updated automatically as coords are read"""
    def CmpArrays(a,b):
      if not len(a) == len(b):
        return False
      return all(a==b)
    nc = len(CoordObjList)
    self.CoordObjList = CoordObjList
    self.NCoordsList = [len(x) for x in self.CoordObjList]
    self.NCoords = sum(self.NCoordsList)
    self.CoordStartInd = [int(sum(self.NCoordsList[:i])) for i in range(nc + 1)]
    self.__OpenObjNum = None
    self.AtomNames = self.CoordObjList[0].AtomNames
    self.NAtom = len(self.AtomNames)
    self.Seq = self.CoordObjList[0].Seq
    self.AtomRes = self.CoordObjList[0].AtomRes
    #check that we have the same sequences and atoms
    for Obj in self.CoordObjList[1:]:
      if not CmpArrays(self.AtomNames, Obj.AtomNames):
        raise ValueError, "Multiple coordinate objects have different atom names."
      elif not CmpArrays(self.Seq, Obj.Seq):
        raise ValueError, "Multiple coordinate objects have different sequences."
      elif not CmpArrays(self.AtomRes, Obj.AtomRes):
        raise ValueError, "Multiple coordinate objects have different residue names."
    #set the linked pos
    self.LinkPos = LinkPos
    #reset
    self.Reset()
    #initial positions
    self.Pos = self[0]
    self.Close()

  def Close(self):
    for Obj in self.CoordObjList:
      Obj.Close()
    self.__OpenObjNum = None

  def Reset(self):
    self.Index = -1
    self.Count = 0
    self.__OpenObjNum = None
    
  def Get(self, ind, Mask = None):
    #check for reverse notation
    if ind < 0:
      ind += self.NCoords
    #check bounds
    if ind < 0 or ind >= self.NCoords:
      raise IndexError, "Index out of bounds for multicoords class."
    else:
      #calculate the absolute index
      self.Index = ind
      #find the obj in which this one belongs
      for (ObjNum, StartInd) in enumerate(self.CoordStartInd[:-1]):
        if ind < self.CoordStartInd[ObjNum+1]: break
      if not self.__OpenObjNum is None and not ObjNum == self.__OpenObjNum:
        self.CoordObjList[self.__OpenObjNum].Close()
      self.__OpenObjNum  = ObjNum
      j = ind - StartInd
      self.Pos = self.CoordObjList[self.__OpenObjNum].Get(j, Mask)
      #update a linked coord array
      if not self.LinkPos is None: self.LinkPos[:,:] = self.Pos
      return self.Pos

  def __getitem__(self, ind):
    return self.Get(ind)

  def __len__(self):
    return self.NCoords

  def __iter__(self):
    self.Reset()
    return self

  def next(self):
    ind = self.Index + 1
    if ind < self.NCoords:
      self.Count += 1
      return self[ind]
    else:
      self.Reset()
      raise StopIteration

  def GetNextCoords(self):
    "Returns the next set of coordinates in the Trj file, or None if the end is reached."
    ind = self.Index + 1
    if ind < self.NCoords:
      self.Count += 1
      result = self[ind]
    else:
      result = None
    return result

  def GetIndices(self):
    "Returns the configuration indices of all read in so far."
    return range(self.NCoords)
  
