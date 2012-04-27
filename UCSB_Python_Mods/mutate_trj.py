#!/usr/bin/env python

import os, sys, time, protein, geometry

USAGE=""" 
Usage:    mutate_trj.py -y     mdcrd {read in} 
                        -x     mdcrd {written}
                        -p     prmtop
                        -mp    mutant prmtop
                        -i     seq input
                        -int   interval
                        -end   endframe 
                        -start startframe

Seq Input     : The file which contains the sequence i.e.
				 1 ALA
                                 2 ALA
			         3 LYS

Mutant prmtop : topology file corresponding to the newly formed trajectory
Interval      : number of frames to be read in between previous frame default=1
Startframe    : frame number to begin reading default=1
Endframe      : frame number to end reading default=9999999"""                        
args =sys.argv 
if len(args)==1:
  print USAGE
  sys.exit()                        

#_____Defaults_____
startframe=1
interval=1
endframe=9999999

#parse through inputs
for input in range(len(args)):
  if args[input]== '-y':
    intrj=args[input+1]
  elif args[input]== '-x':
    outtrj=args[input+1]
  elif args[input]=='-p':
    prmtop=args[input+1]
  elif args[input]=='-int':
    interval=int(args[input+1])
  elif args[input]=='-start':
    startframe=int(args[input+1])
  elif args[input]=='-end':
    endframe=int(args[input+1])
  elif args[input]=='-mp':
    mutantprmtop=args[input+1]
  elif args[input]=='-i':
       seqfile=args[input+1]
  elif args[input].startswith('-'):
    print 'the flag ' + str(args[input]) + ' not recognized\n'
    sys.exit()

#create ptraj input to make the pdbs
ptrajfile=open('ptraj.mutate.in', 'w')
ptrajfile.write('trajin ' + str(intrj) + ' ' + str(startframe) + ' ' + str(endframe) + ' ' + str(interval))
ptrajfile.write('\ntrajout mutant.pdb pdb')
ptrajfile.close()

#run ptraj
print '_____Making pdbs_____'
os.system('ptraj ' + str(prmtop) + ' ptraj.mutate.in > ptraj.mutate.log')

print '\n\n______Mutating pdbs_______'

#get the number of pdbs made
ptrajfile=open('ptraj.mutate.log', 'r')
for line in ptrajfile:
  if 'Coordinate processing' in line:
    words=line.split()
    numpdbs=words[5]
ptrajfile.close()

print 'mutating ' + str(numpdbs) + ' pdbs'

#creates arrays which holds the res number/res name to change
resnum=[]
resname=[]
seq=open(seqfile, 'r')

for line in seq:
  line=line.split()
  resnum.append(line[0])
  resname.append(line[1])
seq.close()

timestart=time.time()
for i in range(1,int(numpdbs)+1):
  Pdb1='mutant.pdb.' + str(i)
  p=protein.ProteinClass(Pdb=Pdb1)
  for num, name in zip(resnum,resname):
    p=p.MutateRes(int(num)-1,name)
  p.WritePdb(Pdb1)
  if (i+1) % 5 ==0: print "....%d pdbs mutated" % (i+1)
timeend=time.time()
print round((timeend-timestart)/60, 2), ' min'

#Create mutant mdcrd
#print '\n_____Creating Mutant Trajectory_____'

ptrajfile=open('mutate.mdcrd.in', 'w')

for i in range(1,int(numpdbs)+1):
  ptrajfile.write('\ntrajin mutant.pdb.' + str(i))

ptrajfile.write('\ntrajout ' + str(outtrj) + ' nobox')
ptrajfile.close()

os.system('ptraj ' + str(mutantprmtop) + ' mutate.mdcrd.in > mutate.mdcrd.log')

###----clean up----
print '\n removing mutant pdbs'
os.system('rm -f mutant.pdb.*')
