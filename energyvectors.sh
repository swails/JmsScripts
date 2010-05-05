#!/bin/sh

# The purpose of this script is to organize all of the energy vectors 
# from an MMPBSA run into two files, one that has all of the bonded
# terms (bond, angle, dihedral, 1-4 electrostatics, 1-4 vdw) and the
# other that has all of the non-bonded terms (egb/epb, electrostatics
# and vdw) in another file.

#######################################################################
#                                                                     #
#                User-defined variables. change as needed             #
#                                                                     #
#######################################################################

# _MMPBSA_ or _MMPBSA_mutant_ ?
prefix="_MMPBSA_"

# gb or pb mdouts? (lowercase only!)
solvent="gb"

# non-bonded output file
nbout="${solvent}_nonbond.dat"

# bonded output file (defaults to /dev/null, but you can put a name here)
bout=/dev/null

#######################################################################
#                                                                     #
#                   Begin the actual script                           #
#                                                                     #
#######################################################################

if [ "`whoami`" = "root" ]; then
   echo "WARNING: DO NOT RUN THIS SCRIPT AS ROOT!"
   exit 1
fi

if [ ! -f ${prefix}complex_${solvent}.mdout ]; then
   echo "Error: Can't find ${prefix}complex_${solvent}.mdout ! MMPBSA"
   echo "       files are missing! Edit energyvectors.sh to specify"
   echo "       solvent model and/or prefix (_MMPBSA_ or _MMPBSA_mutant_)"
   exit 1
fi

# set caps-ed solvent
if [ "$solvent" = "gb" ]; then
   sol="GB"
else
   sol="PB"
fi

# remove old output files
echo "All files are listed in 3 columns: Complex  Receptor  Ligand"
echo "Removing ${nbout} and ${bout}..."
rm -f $nbout $bout .*evtmp

# start with complex: store bond, angle, dihedral, 1-4 eel (eelof), 1-4 vdw
# (vdwof) vdw, eel, and esurf in .___.evtmp files

awk '$1=="BOND"{print $3}' ${prefix}complex_${solvent}.mdout > .bondc.evtmp
awk '$1=="BOND"{print $6}' ${prefix}complex_${solvent}.mdout > .anglec.evtmp
awk '$1=="BOND"{print $9}' ${prefix}complex_${solvent}.mdout > .dihedc.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}complex_${solvent}.mdout > .vdwc.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}complex_${solvent}.mdout > .eelc.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}complex_${solvent}.mdout > .egbc.evtmp
awk '$1=="1-4"{print $4}' ${prefix}complex_${solvent}.mdout > .vdwofc.evtmp
awk '$1=="1-4"{print $8}' ${prefix}complex_${solvent}.mdout > .eelofc.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}complex_${solvent}.mdout > .esurfc.evtmp

awk '$1=="BOND"{print $3}' ${prefix}receptor_${solvent}.mdout > .bondr.evtmp
awk '$1=="BOND"{print $6}' ${prefix}receptor_${solvent}.mdout > .angler.evtmp
awk '$1=="BOND"{print $9}' ${prefix}receptor_${solvent}.mdout > .dihedr.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}receptor_${solvent}.mdout > .vdwr.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}receptor_${solvent}.mdout > .eelr.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}receptor_${solvent}.mdout > .egbr.evtmp
awk '$1=="1-4"{print $4}' ${prefix}receptor_${solvent}.mdout > .vdwofr.evtmp
awk '$1=="1-4"{print $8}' ${prefix}receptor_${solvent}.mdout > .eelofr.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}receptor_${solvent}.mdout > .esurfr.evtmp

awk '$1=="BOND"{print $3}' ${prefix}ligand_${solvent}.mdout > .bondl.evtmp
awk '$1=="BOND"{print $6}' ${prefix}ligand_${solvent}.mdout > .anglel.evtmp
awk '$1=="BOND"{print $9}' ${prefix}ligand_${solvent}.mdout > .dihedl.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}ligand_${solvent}.mdout > .vdwl.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}ligand_${solvent}.mdout > .eell.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}ligand_${solvent}.mdout > .egbl.evtmp
awk '$1=="1-4"{print $4}' ${prefix}ligand_${solvent}.mdout > .vdwofl.evtmp
awk '$1=="1-4"{print $8}' ${prefix}ligand_${solvent}.mdout > .eelofl.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}ligand_${solvent}.mdout > .esurfl.evtmp

paste -d"   " .bondc.evtmp .bondr.evtmp > .evtmp
paste -d"   " .evtmp .bondl.evtmp > .bondt.evtmp
paste -d"   " .anglec.evtmp .angler.evtmp > .evtmp
paste -d"   " .evtmp .anglel.evtmp > .anglet.evtmp
paste -d"   " .dihedc.evtmp .dihedr.evtmp > .evtmp
paste -d"   " .evtmp .dihedl.evtmp > .dihedt.evtmp
paste -d"   " .vdwc.evtmp .vdwr.evtmp > .evtmp
paste -d"   " .evtmp .vdwl.evtmp > .vdwt.evtmp
paste -d"   " .eelc.evtmp .eelr.evtmp > .evtmp
paste -d"   " .evtmp .eell.evtmp > .eelt.evtmp
paste -d"   " .egbc.evtmp .egbr.evtmp > .evtmp
paste -d"   " .evtmp .egbl.evtmp > .egbt.evtmp
paste -d"   " .vdwofc.evtmp .vdwofr.evtmp > .evtmp
paste -d"   " .evtmp .vdwofl.evtmp > .vdwoft.evtmp
paste -d"   " .eelofc.evtmp .eelofr.evtmp > .evtmp
paste -d"   " .evtmp .eelofl.evtmp > .eeloft.evtmp
paste -d"   " .esurfc.evtmp .esurfr.evtmp > .evtmp
paste -d"   " .evtmp .esurfl.evtmp > .esurft.evtmp

# Time to compile everything:
echo "BOND\n" >> $bout
cat .bondt.evtmp >> $bout
echo "\nANGLE\n" >> $bout
cat .anglet.evtmp >> $bout
echo "\nDIHEDRAL\n" >> $bout
cat .dihedt.evtmp >> $bout
echo "VDWAALS\n" >> $nbout
cat .vdwt.evtmp >> $nbout
echo "\nEEL\n" >> $nbout
cat .eelt.evtmp >> $nbout
echo "\nE${sol}\n" >> $nbout
cat .egbt.evtmp >> $nbout
echo "\n1-4 VDW\n" >> $bout
cat .vdwoft.evtmp >> $bout
echo "\n1-4 EEL\n" >> $bout
cat .eeloft.evtmp >> $bout
echo "\nESURF\n" >> $nbout
cat .esurft.evtmp >> $nbout

rm -f .*evtmp

echo "Done compiling energy vectors! non-bonded energies are in ${nbout}"
echo "  bonded energies are in ${bout}!"
