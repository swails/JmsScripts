#!/usr/bin/env python

#LAST MODIFIED: 08-21-08

import copy
import numpy as np

class AminoAcidClass:
  "Class containing information about an amino acid."
  def __init__(self, Name = "", AA1 = "", AA3 = "", Hydrophobic = True,
    Charge = 0, Natural = False, KDHydro = 0., ESGHydro = 0.,
    Cap = False, AA3Aliases = None, PosAtoms = None, NegAtoms = None,
    CapType = "",):
    self.Name = Name
    self.AA1 = AA1
    self.AA3 = AA3
    self.Hydrophobic = Hydrophobic
    self.Charge = Charge
    self.Natural = Natural
    self.Cap = Cap
    self.CapType = CapType
    self.AA3Aliases = AA3Aliases
    if self.AA3Aliases is None: self.AA3Aliases = []
    self.KDHydro = KDHydro
    self.ESGHydro = ESGHydro
    self.PosAtoms = PosAtoms
    self.NegAtoms = NegAtoms
    if self.PosAtoms is None: self.PosAtoms = []
    if self.NegAtoms is None: self.NegAtoms = []
  def Copy(self):
    return copy.deepcopy(self)

class SeqMapClass:
  """Class that maps one sequence to another, finding maximum
     contiguous alignment.  Class values a and b indicate the start
     and stopping residues in Seq1, whereas c and d indicate the start
     and stopping residues in Seq2.  Seq1 can also be a target sequence.
     Note: stopping residue is not inclusive, per Python convention."""
  def __init__(self, Seq1 = None, Seq2 = None,
               a = 0, b = 0, c = 0, d = 0, ID = None,
               Seq1Start = 0, Seq2Start = 0):
    """There are two ways to initialize.  One is to provide 
       Seq1 and Seq2 and the mapping will automatically be performed.
       The other is to provide the actual starting and stopping
       residues of either or both sequences.  ID is optional."""
    self.ID = ID
    self.a = a      #starting res in seq1
    self.b = b      #stopping res in seq1
    self.c = c      #starting res in seq2
    self.d = d      #stopping res in seq2
    if not Seq1 is None and not Seq2 is None:
      self.Map(Seq1, Seq2, Seq1Start, Seq2Start)
  def __repr__(self):
    return "(" + str(self.ID) + ":%d-%d>%d-%d)" % (self.a, self.b-1, self.c, self.d-1)
  def __len__(self):
    if self.a < 0:
      return 0
    else:
      return max(self.b - self.a, 0)
  def Map(self, Seq1, Seq2, Seq1Start = 0, Seq2Start = 0):
    "Finds the maximum overlap between Seq1 and Seq2."
    Seq1, Seq2 = SeqToList(Seq1), SeqToList(Seq2)
    Seq1, Seq2 = Standardize(Seq1), Standardize(Seq2)
    n1, n2 = len(Seq1), len(Seq2)
    if n1 < n2:
      Switch = True
      Seq1, Seq2 = Seq2, Seq1
      n1, n2 = n2, n1
    else:
      Switch = False
    Seq1Str = SeqToAA3(Seq1)
    for l in range(n2, 0, -1):
      for i in range(0, n2 - l + 1):
        Seq2Str = SeqToAA3(Seq2[i:i+l])
        if Seq2Str in Seq1Str:
          ind = Seq1Str.find(Seq2Str)
          self.a = Seq1Str[:ind].count(" ") + Seq1Start
          self.b = self.a + l 
          self.c = i + Seq2Start
          self.d = self.c + l
          if Switch:
            self.a, self.b, self.c, self.d = self.c, self.d, self.a, self.b
          return
  def SubMap(self, a = None, b = None, c = None, d = None):
    """Returns a new SeqMapClass for a subsection in either Seq1 or Seq2 (but not both)."""
    m = copy.deepcopy(self)
    if not a is None or not b is None:
      if a is None: a = self.a
      if b is None: b = self.b
      m.a, m.b = max(self.a, a), min(self.b, b)
      if m.b <= m.a:
        m.a, m.b, m.c, m.d = 0, 0, 0, 0
      else:
        m.c = self.c + m.a - self.a
        m.d = self.d + m.b - self.b
    elif not b is None or not c is None:
      if c is None: c = self.c
      if d is None: d = self.d
      m.c, m.d = max(self.c, c), min(self.d, d)
      if m.d <= m.c:
        m.a, m.b, m.c, m.d = 0, 0, 0, 0
      else:
        m.a = self.a + m.c - self.c
        m.b = self.b + m.d - self.d    
    return m
  def TrimSeq1(self, Seq1):
    "Returns the overlap span with the Seq1 sequence."
    return SeqToList(Seq1)[self.a:self.b]
  def TrimSeq2(self, Seq2):
    "Returns the overlap span with the Seq2 sequence."
    return SeqToList(Seq2)[self.c:self.d]
  def __getattr__(self, attr):
      if attr == "Take1Ind":
          return np.arange(self.a, self.b).astype(int)
      elif attr == "Take2Ind":
          return np.arange(self.c, self.b).astype(int)
      else:
          raise AttributeError(attr)



class SeqMap2Class:
  """Class that maps one sequence to another, finding maximum
     contiguous alignment.  Allows the second sequence to have missing
     residues/gaps.  Class boolean arrays Take1 and Take2 indicate
     residues in the sequence to use that are aligned."""
  def __init__(self, Seq1, Seq2, MinContig = 7):
    """Will map Seq1 and Seq2 to each other; Seq2 can have gaps.
MinContig gives the minimum number of contigous residues that
must be found to match."""
    def ToStr(Seq):
        return "".join(["%-5s" % x for x in Seq])
    Seq1, Seq2 = SeqToList(Seq1), SeqToList(Seq2)
    Seq1, Seq2 = Standardize(Seq1), Standardize(Seq2)
    Seq1Str = ToStr(Seq1)
    Seq2Str = ToStr(Seq2)
    n1, n2 = len(Seq1), len(Seq2) 
    self.Take1 = np.zeros(n1, bool)
    self.Take2 = np.zeros(n2, bool)
    i1 = 0
    i2 = 0
    while i1 < n1 and i2 < n2:
        l = MinContig
        r1a, r1b = 0, 0
        r2a, r2b = 0, 0
        lastind = None
        while i2 + l <= n2:
            ind = Seq1Str.find(Seq2Str[:5*l])
            if ind < 0: 
                break
            if not (lastind is None or ind==lastind):
                break
            lastind = ind
            r1a = ind / 5
            r1b = r1a + l
            r2a = 0
            r2b = l
            l += 1
        l = r1b - r1a
        if l == 0:
            i2 += 1
            Seq2Str = Seq2Str[5:]
        else:
            self.Take1[i1+r1a:i1+r1b] = 1
            Seq1Str = Seq1Str[5*r1b:]
            i1 += r1b
            self.Take2[i2+r2a:i2+r2b] = 1
            Seq2Str = Seq2Str[5*r2b:]
            i2 += r2b
    self.n = int(np.sum(self.Take1))
    #make a string representation
    self.Str = []
    for (i1, i2) in zip(self.Take1Ind, self.Take2Ind):
        self.Str += ["%s%d -> %s%d" % (Seq1[i1], i1, Seq2[i2], i2)]
    self.Str = "\n".join(self.Str)
    
  def __repr__(self):
    return "(Map2:%d res)" % self.n
  def __len__(self):
    return self.n
  def __getattr__(self, attr):
      if attr == "Take1Ind":
          return np.nonzero(self.Take1)[0]
      elif attr == "Take2Ind":
          return np.nonzero(self.Take2)[0]
      else:
          raise AttributeError(attr)
  def Switch(self):
      """Switches identities of Seq1, Seq2"""
      self.Take1, self.Take2 = self.Take2, self.Take1
      self.Str = [x.split(" -> ") for x in self.Str.split("\n")]
      self.Str = "\n".join(["%s -> %s" % tuple(x) for x in self.Str])
          
def GetSeqMap2Class(Seq1, Seq2):
    """Returns a SeqMap2Class for whichever has the highest match."""
    Map1 = SeqMap2Class(Seq1, Seq2)
    Map2 = SeqMap2Class(Seq2, Seq1)
    if len(Map1) >= len(Map2):
        return Map1
    else:
        Map2.Switch()
        return Map2
      


def GetAllAtomAlign(p1, p2, Backbone = False, ResInd = None):
    """Returns an alignment between two protein classes at the all-atom level."""
    #align the sequences
    Map = SeqMapClass(p1.Seq, p2.Seq)
    p1 = p1[Map.a:Map.b]
    p2 = p2[Map.c:Map.d]
    #filter res indices
    n1, n2 = len(p1), len(p2)
    if not ResInd is None: 
        if not n1 == n2:
            raise ValueError("ProteinClasses have different number of residues.")
        ResInd = [i for i in ResInd if i < n1]
    #get the right atoms
    if Backbone:
        AtomMask = ["N", "CA", "C"]
    else:
        AtomMask = None
    #find overlap atoms
    AtomInd1, AtomInd2 = [], []
    for (i, a1) in enumerate(p1.Atoms):
        ResNum = p1.AtomResNum[i]
        InResInd = True
        if not ResInd is None: InResInd = ResNum in ResInd
        j = p2.AtomNum(ResNum, a1.Name, NotFoundError = False)
        if j >= 0:
            if InResInd:
                AtomInd1.append(i)
                AtomInd2.append(j)
    AtomInd1 = np.array(AtomInd1, int)
    AtomInd2 = np.array(AtomInd2, int)
    return AtomInd1, AtomInd2


def GetAllMaps(Seq1, Seq2, MinLen = 3, ID = None, Seq1Start = 0):
  """Returns a list of all disjoint mappings of Seq2 to Seq1"""
  if len(Seq1) == 0: return []
  Seq1, Seq2 = SeqToList(Seq1), SeqToList(Seq2)
  Map = SeqMapClass(Seq1 = Seq1, Seq2 = Seq2, ID = ID, Seq1Start = Seq1Start)
  if len(Map) < MinLen:
    return []
  else:
    return GetAllMaps(Seq1[:Map.a - Seq1Start], Seq2, MinLen, ID, Seq1Start) + \
           [Map] + \
           GetAllMaps(Seq1[Map.b - Seq1Start:], Seq2, MinLen, ID, Map.b)

def RecurseAlign(StartRes, StopRes, SeqMaps, MinLen = 3, MaxTrim = None,
                 MaximizeMaps = False):
  """Produces a list of lists of mappings of sequences in SeqMaps;
   StartRes is inclusive, StopRes is not."""
  #check that target sequence is zero length
  if StopRes <= StartRes:
    return [[]]
  #if target is too small, return non matched sequence
  if StopRes - StartRes < MinLen:
    return [[SeqMapClass(a = StartRes, b = StopRes)]]  
  #start with the base alignment of non matched
  AlignList = [[SeqMapClass(a = StartRes, b = StopRes)]]
  #check all sequence possibilities for alignments
  for SeqMap in SeqMaps:
    Map = SeqMap.SubMap(a = StartRes, b = StopRes)
    #if the alignment took...
    if len(Map) >= MinLen:
      #check to see if we didn't trim too much 
      if not MaxTrim is None:
        NTrim = (SeqMap.b - SeqMap.a) - (Map.b - Map.a)
        if NTrim > MaxTrim: continue
      #search for more
      LeftList = RecurseAlign(StartRes, Map.a, SeqMaps, MinLen, MaxTrim, MaximizeMaps)
      RightList = RecurseAlign(Map.b, StopRes, SeqMaps, MinLen, MaxTrim, MaximizeMaps)
      for l in LeftList:
        for r in RightList:
          AlignList.append(l + [Map] + r)
  #check if we only want the maximum number of maps
  if MaximizeMaps and len(AlignList) > 1: AlignList = AlignList[1:]
  #remove duplicates
  i = 0
  StrAlignList = [str(x) for x in AlignList]
  while i < len(AlignList) - 1:
    StrAlign = str(AlignList[i])
    if StrAlign in StrAlignList[i+1:]:
      ind = StrAlignList.index(StrAlign, i+1)
      del AlignList[ind]
      del StrAlignList[ind]
    else:
      i += 1
  return AlignList

def GetAllAlign(TarSeq, Seqs, MultiMap = False, MinLen = 3,
                MaxTrim = None, MaximizeMaps = False):
  """Returns all alignments of Seqs to TarSeq.
* MultiMap = True will allow each Seq to be mapped to TarSeq
  in multiple locations.
* MinLen is the minimum number of residues that will be mapped.
* MaxTrim is the maximum number of residues that will be trimmed
  from the ends of each seq in Seqs to do the alignment
* MaximizeMaps = True will only return alignments with the
  maximum number of alignments (i.e., minimal unaligned regions)
  """
  #make maps
  SeqMaps = []
  for (i,Seq) in enumerate(Seqs):
    if MultiMap:
      Maps = GetAllMaps(TargetSeq, Seq, MinLen = MinLen, ID = i)
      SeqMaps.extend(Maps)
    else:
      Map = SeqMapClass(Seq1 = TarSeq, Seq2 = Seq, ID = i)
      if len(Map) >= MinLen: SeqMaps.append(Map)
  #perform alignments
  AlignList = RecurseAlign(0, len(TarSeq), SeqMaps, MinLen = MinLen,
                           MaxTrim = MaxTrim, MaximizeMaps = MaximizeMaps)
  #remove the fully unmatched alignment
  if not MaximizeMaps: AlignList = AlignList[1:]
  return AlignList



#List of standard amino acids; lower index means higher priority
#for aliases, searching, etc.
AminoAcids = [
  AminoAcidClass("glycine", "G", "GLY", True, 0, True, 4.5, 1.0),
  AminoAcidClass("alanine", "A", "ALA", True, 0, True, 1.8, 1.6),
  AminoAcidClass("valine", "V", "VAL", True, 0, True, 4.2, 2.6),
  AminoAcidClass("leucine", "L", "LEU", True, 0, True, 3.8, 2.8),
  AminoAcidClass("isoleucine", "I", "ILE", True, 0, True, 4.5, 3.1),
  AminoAcidClass("methionine", "M", "MET", True, 0, True, 1.9, 3.4),
  AminoAcidClass("phenylalanine", "F", "PHE", True, 0, True, 2.8, 3.7),
  AminoAcidClass("tryptophan", "W", "TRP", True, 0, True, -0.9, 1.9),
  AminoAcidClass("proline", "P", "PRO", True, 0, True, -1.6, -0.2),
  AminoAcidClass("cysteine", "C", "CYS", True, 0, True, 2.5, 2.0,
                 AA3Aliases = ["CYM", "CYX"]),
  AminoAcidClass("threonine", "T", "THR", False, 0, True, -0.7, 1.2),
  AminoAcidClass("serine", "S", "SER", False, 0, True, -0.8, 0.6),
  AminoAcidClass("tyrosine", "Y", "TYR", True, 0, True, -1.3, -0.7),
  AminoAcidClass("asparagine", "N", "ASN", False, 0, True, -3.5, -4.8),
  AminoAcidClass("glutamine", "Q", "GLN", False, 0, True, -3.5, -4.1),
  AminoAcidClass("aspartic acid", "D", "ASP", False, -1, True, -3.5, -9.2,
                 NegAtoms = ["OD1", "OD2"], AA3Aliases = ["ASH"]),
  AminoAcidClass("glutamic acid", "E", "GLU", False, -1, True, -3.5, -8.2,
                 NegAtoms = ["OE1", "OE2"], AA3Aliases = ["GLH"]),
  AminoAcidClass("lysine", "K", "LYS", False, 1, True, -3.9, -8.8,
                 PosAtoms = ["NZ"]),
  AminoAcidClass("arginine", "R", "ARG", False, 1, True, -4.5, -12.3,
                 PosAtoms = ["NH1", "NH2"]),
  AminoAcidClass("histidine", "H", "HIS", False, 0, True, -3.2, -3.0,
                 AA3Aliases = ["HIE", "HID", "HIP"]),
  AminoAcidClass("acetyl", ">", "ACE", False, 0, True, Cap = True, CapType = "N"),
  AminoAcidClass("amine", "{", "NHE", False, 0, True, Cap = True, CapType = "C"),
  AminoAcidClass("n-methylamine", "<", "NME", False, 0, True, Cap = True, CapType = "C"),
  AminoAcidClass("aspartic acid", "D", "ASH", False, 0, True, -3.5, -9.2),
  AminoAcidClass("cysteine", "C", "CYM", True, 0, True, 2.5, 2.0),
  AminoAcidClass("cysteine", "C", "CYX", True, 0, True, 2.5, 2.0),
  AminoAcidClass("glutamic acid", "E", "GLH", False, 0, True, -3.5, -8.2),
  AminoAcidClass("histidine", "H", "HIE", False, 0, True, -3.2, -3.0),
  AminoAcidClass("histidine", "H", "HID", False, 0, True, -3.2, -3.0),
  AminoAcidClass("histidine", "H", "HIP", False, 1, True, -3.2, -3.0,
                 PosAtoms = ["ND1", "NE2"])  
  ]
#Add N and C-terminal aliases and salt atoms
for aa in AminoAcids:
  AllNames = [aa.AA3] + aa.AA3Aliases
  for n in AllNames:
    aa.AA3Aliases.extend(["N" + n, "C" + n])


def AAInst(a):
  """Returns the amino acid class instance for a, which can be 1 or 3 letters.
  First checks all names, then checks aliases."""
  if len(a) == 1:
    l = [aa for aa in AminoAcids if aa.AA1 == a]
    if len(l) > 0:
      return l[0]
  else:
    l = [aa for aa in AminoAcids if aa.AA3 == a]
    if len(l) > 0:
      return l[0]
    else:
      l = [aa for aa in AminoAcids if a in aa.AA3Aliases]
      if len(l) > 0:
        return l[0]      

def AAInstAlias(a):
  """Returns the amino acid class instance for a, which can be 1 or 3 letters.
  Checks names and aliases concurrently."""
  if len(a) == 1:
    l = [aa for aa in AminoAcids if aa.AA1 == a]
    if len(l) > 0:
      return l[0]
  else:
    l = [aa for aa in AminoAcids if aa.AA3 == a or a in aa.AA3Aliases]
    if len(l) > 0:
      return l[0]

def AreAliases(a, b):
  """Returns True if a is an alias of b or vice versa."""
  if a == b: return True
  a = AAInst(a)
  b = AAInst(b)
  if a is None or b is None: return False
  return a.AA3 == b.AA3 or a.AA3 in b.AA3Aliases or b.AA3 in a.AA3Aliases

def Natural(a):
  "Returns True if a is a natural amino acid."
  aa = AAInst(a)
  if aa is None:
    return False
  else:
    return aa.Natural

def Exists(a):
  "Returns True if a is found in the AA database."
  return (not AAInst(a) is None)

def Charge(a):
  aa = AAInst(a)
  if aa is None:
    return 0
  else:
    return aa.Charge

def Charged(a):
  aa = AAInst(a)
  if aa is None:
    return False
  else:
    return (aa.Charge != 0)

def Hydrophobic(a):
  "Returns True if a is hydrophobic or not found."
  aa = AAInst(a)
  if aa is None:
    return True
  else:
    return aa.Hydrophobic

def Cap(a):
  "Returns True if a is a cap."
  aa = AAInst(a)
  if aa is None:
    return False
  else:
    return aa.Cap

def CCap(a):
  "Returns True if a is a C - cap."
  aa = AAInst(a)
  if aa is None:
    return False
  else:
    return aa.Cap and aa.CapType == "C"

def NCap(a):
  "Returns True if a is a N - cap."
  aa = AAInst(a)
  if aa is None:
    return False
  else:
    return aa.Cap and aa.CapType == "N"   

def Terminal(a):
  "Returns True if a is a terminal residue."
  a = a.strip()
  return len(a) > 3 and a[0] in "NCnc"

def NTerminal(a):
  "Returns True if a is a terminal residue."
  a = a.strip()
  return len(a) > 3 and a[0] in "Nn"

def CTerminal(a):
  "Returns True if a is a terminal residue."
  a = a.strip()
  return len(a) > 3 and a[0] in "Cc"

def NTerm(a):
  "Converts a to a N-terminal residue."
  if Terminal(a):
    return copy.copy(a)
  elif len(a) >= 3:
    if a[0].isupper():
      return "N" + a
    else:
      return "n" + a
  else:
    return a

def CTerm(a):
  "Converts a to a C-terminal residue."
  if Terminal(a):
    return copy.copy(a)
  elif len(a) >= 3:
    if a[0].isupper():
      return "C" + a
    else:
      return "c" + a
  else:
    return a  

def NonTerm(a):
  "Converts a to a normal residue, if terminal."
  if Terminal(a) and len(a) >= 3:
    return a[1:]
  else:
    return a
  
def AttrPair(a, b):
  """Returns True if a and b are likely to be mutually attractive
     residues or if they are not found."""
  aa1 = AAInst(a)
  aa2 = AAInst(b)
  if aa1 is None or aa2 is None:
    return True
  else:
    return ((aa1.Charge * aa2.Charge == -1) or (aa1.Hydrophobic and aa2.Hydrophobic))

def HydrophobicPair(a, b):
  """Returns True if a and b are likely to be a hydrophobic pair
     or if they are not found."""
  aa1 = AAInst(a)
  aa2 = AAInst(b)
  if aa1 is None or aa2 is None:
    return True
  else:
    return (aa1.Hydrophobic and aa2.Hydrophobic)

def SaltPair(a, b):
  """Returns True if a and b can form a salt bridge."""
  aa1 = AAInst(a)
  aa2 = AAInst(b)
  if aa1 is None or aa2 is None:
    return False
  else:
    return (aa1.Charge * aa2.Charge < 0)  


def ParseList(Seq, MatchAliases = False):
  "Parses a list of 3-letter codes and matches any aliases."
  if MatchAliases:
    Seq1 = []
    for r in Seq:
      aa = AAInst(r)
      if aa is None:
        Seq1.append(r)
      else:
        Seq1.append(aa.AA3)
    return Seq1
  else:
    return copy.copy(Seq)

def AA3ToAA1(Seq):
  "Converts a string of 3-letter codes to a string of 1-letter codes."
  Seq1 = ""
  for r in Seq.split():
    aa = AAInst(r)
    if aa is None:
      Seq1 = Seq1 + "?"
    else:
      Seq1 = Seq1 + aa.AA1
  return Seq1

def AA1ToAA3(Seq):
  "Converts a string of 1-letter codes to a string of 3-letter codes."
  Seq3 = ""
  for r in Seq:
    aa = AAInst(r)
    if aa is None:
      Seq3 = Seq3 + "??? "
    else: 
      Seq3 = Seq3 + aa.AA3 + " " 
  return Seq3.strip()

def SeqToList(Seq, MatchAliases = False):
  "Converts an arbitrary sequence type to a list of 3-letter codes."
  #check if Seq is in list format
  if type(Seq) is list:
    return ParseList(Seq, MatchAliases)
  #check if Seq is in AA3 format
  elif " " in Seq:
    return ParseList(Seq.split(), MatchAliases)
  #Seq must be AA1 format
  else:
    return AA1ToAA3(Seq).split()

def SeqToAA3(Seq, MatchAliases = False):
  "Converts an arbitrary sequence type to a string of 3-letter codes."
  #check if Seq is in list format
  if type(Seq) is list:
    return " ".join(ParseList(Seq, MatchAliases))
  #check if Seq is in AA3 format
  elif " " in Seq:
    return " ".join(ParseList(Seq.split(" "), MatchAliases))
  #Seq must be AA1 format
  else:
    return AA1ToAA3(Seq) 

def SeqToAA1(Seq):
  "Converts an arbitrary sequence type to a string of 1-letter codes."
  #check if Seq is in list format
  if type(Seq) is list:
    return AA3ToAA1(" ".join(Seq))
  #check if Seq is in AA3 format
  elif " " in Seq:
    return AA3ToAA1(Seq)
  #Seq must be AA1 format
  else:
    return copy.copy(Seq)


def AllNatural(Seq):
  "Returns True if all amino acids in Seq are canonical."
  return ([Natural(a) for a in SeqToList(Seq)].count(False) == 0)

def AllExist(Seq):
  "Returns True if all amino acids in Seq exist in the AA database."
  return ([Exists(a) for a in SeqToList(Seq)].count(False) == 0)


def Standardize(Seq):
  "Returns a list with aliases replaced by standard names."
  Seq1 = []
  for r in SeqToList(Seq):
    aa = AAInstAlias(r)
    if aa is None:
      Seq1.append(r)
    else:
      Seq1.append(aa.AA3)
  return Seq1
  
    