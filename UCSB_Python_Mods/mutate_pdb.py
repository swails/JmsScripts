#!/usr/bin/env python


Usage= """  mutate_pdb.py -i InPdb -o OutPdb -resnum ResNum -resname ResName -seq sequence_input

InPdb:		input pdb file
OutPdb: 	output pdb file
ResNum: 	Residue number you want mutated
ResName: 	Name of the residue you want it mutated to
Sequence Input: reads in a list of mutations if there is more than one
optional cmd
format: "25 ASN"
minmize:	minimizes the mutated resiudue and everything within 10 A

"""
import protein, geometry, sys


#defaults
INPDB=False
Minimize=False
Cutoff=10
Seq=False
resnum=[]
resname=[]
#get arguments
#---------------------
args=sys.argv
if len(args) ==1:
  print Usage
  sys.exit()


for x in range(len(args)):
  if args[x]=='-i':
    Pdb1=args[x+1]
  elif args[x]=='-o':
    Pdb2=args[x+1]
  elif args[x]=='-resname':
    ResName=args[x+1]
  elif args[x]=='-resnum':
    ResNum=int(args[x+1])-1
  elif args[x]=='-seq':
    try:
      file=open(args[x+1],'r')
    except IOError:
      print args[x+1], 'does not exist'
      sys.exit()
    seqlist=file.readlines()
    file.close()
    Seq=True
  elif args[x]=='minimize':
    Minimize=True

#initiate program----
if Seq:
  p=protein.ProteinClass(Pdb=Pdb1)
  for line in seqlist:
    line=line.split()
    resnum.append(line[0])
    resname.append(line[1])
  for num,name in zip(resnum,resname):
    p=p.MutateRes(int(num)-1,name)
else:
  p=protein.ProteinClass(Pdb=Pdb1)
  p=p.MutateRes(ResNum,ResName)

if Minimize:
#create a list of all the residues that are within 10 Ang of the residue including itself
  Pos = p.ResPos()
  ResInd=[]
  for i in range(len(p)):
    if geometry.Length(Pos[i] - Pos[ResNum]) <= Cutoff:
      ResInd.append(i)
#first optimize just mutated side chain
  p.OptimizeSC(ResInd = [ResNum])

#optimize everything within 10 Ang of mutated residue
  p.OptimizeSC(ResInd = ResInd)

#check for overlap
#if p.HasOverlap():
#  print "Overlap found"
#else:
#  print "No Overlaps"

#save new file
p.WritePdb(Pdb2)





