#!/usr/bin/env python

import sys, os, math, time
sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

Usage=""" antebuilder.py -pdb <input.pdb> 
                -mol2 <mol2_file> 
                -frc<name_of_frcmod> 
                -nc <charge--optional default is 0> 
                -c <charge method> 
                -rn <resname--default LIG>
--------------------------------------------------------------------------------

"""

args=sys.argv
if len(args) < 5:
  print >> sys.stdout,  Usage
  sys.exit()


#--Check for the existence of antechamber
Paths= os.environ["PATH"].split(":")
antechamber=False
for executable in Paths:
  antechamber= os.path.isfile(os.path.join(executable, "antechamber"))
  if antechamber: 
    print >> sys.stdout, "\n Antechamber found: beginning the program\n"
    break
if not antechamber:
  print >> sys.stderr, """\n Antechamber not found: Check that AmberTools was properly 
installed and $AMBERHOME is defined: """
  sys.exit()

#defaults
charge=0
resname='LIG'
#-----------------------------------------------------------------------------

#parse arguements
for input in range(len(args)):
  if args[input]=='-pdb':
    inpdb=args[input+1]
  elif args[input]=='-mol2':
    mol2_file=args[input+1]
  elif args[input]=='-nc':
    charge=args[input+1]
  elif args[input]=='-rn':
    resname=args[input+1]
  elif args[input]=='-frc':
    frcmod_file=args[input+1]
  elif args[input]=='-c':
    chrg_meth=args[input+1]
  elif args[input-1] != '-nc' and args[input].startswith('-'): 
    print "'%s' flag is not recognized\n" % args[input]
    print Usage
    sys.exit() 

#----- Begin to call antechamer -----


# check to make the pdb exists
pdbexists=os.path.isfile(inpdb)
if not pdbexists:
  print >> sys.stderr, "'%s' does not exist" % inpdb
  sys.exit()

### print to user what was inputed ####
print >> sys.stdout, "==============================================="
print >> sys.stdout, "Name of the input pdb is   %s" % inpdb
print >> sys.stdout, "Name of the mol2 file is   %s" % mol2_file
print >> sys.stdout, "Name of the frcmod is      %s" % frcmod_file
print >> sys.stdout, "Charge of your molecule =  %s" % charge
print >> sys.stdout, "Name of the ligand =       %s" % resname
print >> sys.stdout, "==============================================="
#--------------

print >> sys.stdout, "\n____Running Antechamber____"

#running antechamber
startcalc=time.time()
os.system('antechamber -i %s -fi pdb -o %s -fo mol2 -c %s -s 2 -nc %s -rn %s' % (inpdb, mol2_file, chrg_meth, charge, resname))

## check for the existence of the prep file
mol2file=os.path.isfile(mol2_file)
if not mol2file:
  print >> sys.stderr, "\n Error occured mol2_file not generated"
  sys.exit()


print >> sys.stdout, "\n_____Running Parmchk_____\n"
os.system(' parmchk -i %s -f mol2 -o %s' % (mol2_file, frcmod_file) )
endcalc=time.time()
total_time=round((endcalc-startcalc)/60, 2)
print >> sys.stdout, "\nTotal time of the calculation= %s min" % total_time
