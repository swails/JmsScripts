#!/usr/bin/env python

#LAST MODIFIED: 06-25-08


Usage = """\nExtracts pdb snapshots from a trajectory that maintain
           a specific residue-residue distance.

Usage     : extractconfs.py -y      {IN_TRJFILE} 
                            -prm    {PRMTOP} 
                            -prefix {PREFIX} 
                            -res1   {RES1} 
                            -res2   {RES2} 
                            -min    {MINDIST} 
                            -max    {MAXDIST}
                            -x      {OUT_TRJFILE}
                   {OPTIONS -skip   {initial frame}   
                            -int    {interval}
IN_TRJFILE    : trajectory CRD file (can be gzipped)
PRMTOPFILE    : PARM7 file [OPTIONS]
PREFIX        : prefix of output pdb files (will add 001.pdb, 002.pdb, etc)
RES1          : residue number of first residue (starting at 1)
RES1          : residue number of second residue (starting at 1)
MINDIST       : minimum pairwise distance between specified residues
MAXDIST       : maximum pairwise distance between specified residues
OUT_TRJFILE   : trajectory of just the frames which meet the criterion
initial frame : frame number to begin reading [default =0]
interval      : number of frames to be read in between previous frame [default=1]
"""

#check for instructions
import sys
if __name__ == "__main__" and len(sys.argv) < 3:
  print Usage
  sys.exit()


#import numerics library
from numpy import *
#standard python modules  
import os, sys
#modules developed by the Shell group
import protein, coords, geometry

#______PARSE Arguments_____

#____Defaults______
nskip=0
nstride=1


ARGS=sys.argv
for input in range(len(ARGS)):
  if ARGS[input].startswith('-y'):
    Trjfn=str(ARGS[input+1])
  elif ARGS[input].startswith('-prm'):
    Prmtopfn=str(ARGS[input+1])
  elif ARGS[input].startswith('-pref'):
    Prefix=ARGS[input+1]
  elif ARGS[input].startswith('-res1'):
    ResNum1=int(ARGS[input+1])-1  
  elif ARGS[input].startswith('-res2'):
    ResNum2=int(ARGS[input+1])-1
  elif ARGS[input].startswith('-min'):
    MinDist=float(ARGS[input+1])
  elif ARGS[input].startswith('-max'):
    MaxDist=float(ARGS[input+1])
  elif ARGS[input].startswith('-x'):
    mdcrd=ARGS[input+1]
  elif ARGS[input].startswith('-skip'):
    nskip=int(ARGS[input+1])
  elif ARGS[input].startswith('-int'):
    nstride=int(ARGS[input+1])
  elif ARGS[input].startswith('-'):
    print str(ARGS[input]) + " flag is not recognized\n"
    sys.exit()


#open the trajectory file
Trj = coords.TrjClass(str(Trjfn),str(Prmtopfn), NSkip=nskip, NStride=nstride)
NDigits = int(log(len(Trj) + 0.1) / log(10.) + 1.)
print "Found %d snapshots in trajectory file" % len(Trj)

#link a proteinclass object
p = protein.ProteinClass()
p.LinkTrj(Trj)
pdbfiles=[]

#check residues
if ResNum1 > len(p):
  print "Error: first residue number exceeds protein length."
  sys.exit()
if ResNum2 > len(p):
  print "Error: second residue number exceeds protein length."
  sys.exit()
#get CA indices
CA1 = p.AtomNum(ResNum1, "CA")
CA2 = p.AtomNum(ResNum2, "CA")

#print out results
print "First atom is : CA%d in residue %s%d" % (CA1+1, p.Seq[ResNum1], ResNum1+1)
print "Second atom is: CA%d in residue %s%d" % (CA2+1, p.Seq[ResNum2], ResNum2+1)
print "Distance range is from %8.3f to %8.3f A" % (MinDist, MaxDist)

#now sort through all of the configurations
Ind = 0
for (i, Pos) in enumerate(Trj):
  if (i+1) % 500 == 0: print "...examined %d snapshots" % (i+1)
  Dist = geometry.Length(Pos[CA1] - Pos[CA2])
  if Dist >= MinDist and Dist < MaxDist:
    fn = "%s%0*d.pdb" % (str(Prefix), NDigits, Ind)
    p.WritePdb(fn)
    Ind += 1
    pdbfiles.append(fn)
#close down
p.UnlinkTrj()
Trj.Close()
#Creates the mdcrd of the frames the meet the criterion
print " Writing " + str(len(pdbfiles)) + " frames to " + str(mdcrd)

coords.PdbsToCrd(pdbfiles,str(mdcrd))
    

