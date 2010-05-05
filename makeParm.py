#!/usr/bin/env python

#######################################################################
#                                                                     #
# This script will create an amber prmtop taking command-line options #
# using tleap. More options will be added as they're needed.          #
#                                                                     #
#    Written by Jason Swails (last updated 12-19-2009)                #
#                                                                     #
#######################################################################

import sys, utilities, os 
def printusage():
   print "Usage: makeParm.py -ff frc_fld \\ (Default leaprc.ff99SB)"
   print "                   -ci 0,1,or 2\\ (Counterions, default 0)"
   print "                   -sol dist   \\ (Solvate shell distance,"  \
       + " Default 0)"
   print "                   -cph \\ (Turns on constant pH setup)"
   print "                   -pdb pdbfile\\ (no default)"
   print "                   -run \\ (Runs tleap after script written)"
   print "                   -radii option \\ (Default \"\")"
   print "Last option asks whether to run tleap automatically, " +     \
         "or just create the leap.in script."


if len(sys.argv) == 1:
   printusage()
   sys.exit()
if "-help" in sys.argv[1].lower():
   printusage()
   sys.exit()

frc_fld = "leaprc.ff99SB"
ci      = 0
sol     = 0
cph     = "no"
run     = "no"
pdb     = ""
radii   = ""

infiles = False

try:
   for x in range(len(sys.argv)):
      if sys.argv[x].lower() == "-ff":
         frc_fld = sys.argv[x+1]
      elif sys.argv[x].lower() == "-ci":
         ci = int(sys.argv[x+1])
      elif sys.argv[x].lower() == "-sol":
         sol = float(sys.argv[x+1])
      elif sys.argv[x].lower() == "-run":
         run = "yes"
      elif sys.argv[x].lower() == "-pdb":
         pdb = sys.argv[x+1]
      elif sys.argv[x].lower() == "-cph":
         cph = "yes"
except IndexError:
   print >> sys.stderr, "Error: Command line error!"
   sys.exit()
except ValueError:
   print >> sys.stderr, "Error: Command line error!"
   sys.exit()

print "Removing py_leap.in"
os.system("rm -f py_leap.in")

leapfile = open("py_leap.in","w")
leapfile.write("source " + frc_fld + "\n")


#### BEGIN INPUT FILE TYPES ####

if pdb != "":
   if utilities.fileexists(pdb) == -1:
      leapfile.close()
      os.system("rm -f py_leap.in")
      sys.exit()
   leapfile.write("l = loadpdb " + pdb + "\n")
   infiles = True


#### END INPUT FILE TYPES ####
if not infiles:
   print "Error: You have not specified any structures!"
   leapfile.close()
   os.system("rm -f py_leap.in")
   sys.exit()

if ci == 1:
   leapfile.write("addIons l Na+\n")
elif ci == 2:
   if pdb != "":
      residue_decmp = utilities.getresdecmp_pdb(pdb)

   if len(residue_decmp) == 0:
      print "Not ready yet..."
      leapfile.close()
      os.system("rm -f py_leap.in")
      sys.exit()

   leapfile.write("addIons l Na+ " + str(residue_decmp[21]) +          \
            " Cl- " + str(residue_decmp[20]) + "\n")

if sol > 0:
   leapfile.write("solvateOct l TIP3PBOX " + str(sol) + "\n")

if radii != "" and cph == "no":
   leapfile.write("set default PBRadii " + radii + "\n")

if cph == "yes":
   if sol > 0:
      print "Error: constant pH and explicit solvent are incompatible!"
      leapfile.close()
      os.system("rm -f py_leap.in")
      sys.exit()
   if radii != "mbondi2" and radii != "":
      print "Warning: Reference energies calculated with mbondi2. " +  \
            " MBONDI2 radii will be used!"
   leapfile.write("set default PBRadii mbondi2\n")
   leapfile.write("loadOFF constph.lib\n")
   leapfile.write("loadAmberParams frcmod.constph\n")

leapfile.write("saveAmberParm l prmtop inpcrd\nquit\n")
leapfile.close()

if run == "yes":
   os.system("sleap -f py_leap.in")

