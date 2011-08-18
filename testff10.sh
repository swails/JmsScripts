#! /bin/sh

# This script is supposed to automate the checking of ff10 to make sure
# that it works as advertised in the AmberTools manual.

# First thing we have to do is create the tleap script and run tleap

if [ ! -z `ls | head -n 1` ]; then
   echo "Run this script in an empty directory!"
   exit 1
fi

if [ -z `which sander` -o -z `which tleap` -o -z `which ambpdb` \
     -o -z `which ante-MMPBSA.py` ]; then
   echo "sander, tleap, ambpdb, and ante-MMPBSA.py are necessary!"
   exit 1
fi

echo "Creating the initial structures"
cat > tleap.in << EOF
source leaprc.ff10

set default pbradii mbondi2

l = sequence {ACE ASP ALA ASH GLY LYS THR VAL GLU HID TYR HIE LYN ILE HIP MET CYS PRO CYM ARG NME}

k = sequence {DA5 DC DT DG DA DC3}

m = sequence {A5 C U G A C3}

saveamberparm l protein_ff10.parm7 protein_ff10.rst7
saveamberparm k dna_ff10.parm7 dna_ff10.rst7
saveamberparm m rna_ff10.parm7 rna_ff10.rst7

quit
EOF

tleap -f tleap.in > /dev/null 2>&1

cat > min.mdin << EOF
Minimization to relax initial bad contacts, implicit solvent
 &cntrl
   imin=1,
   ncyc=100,
   maxcyc=200,
   ntpr=100,
   ntb=0,
   cut=1000,
   igb=5,
 /
EOF

echo "Relaxing protein"
sander -O -i min.mdin -p protein_ff10.parm7 -c protein_ff10.rst7 \
          -r protein_ff10_min.rst7

echo "Relaxing dna"
sander -O -i min.mdin -p dna_ff10.parm7 -c dna_ff10.rst7 \
          -r dna_ff10_min.rst7

echo "Relaxing rna"
sander -O -i min.mdin -p rna_ff10.parm7 -c rna_ff10.rst7 \
          -r rna_ff10_min.rst7

echo "Making PDBs"
ambpdb -p protein_ff10.parm7 < protein_ff10_min.rst7 > protein.pdb 2>/dev/null
ambpdb -p dna_ff10.parm7 < dna_ff10_min.rst7 > dna.pdb 2>/dev/null
ambpdb -p rna_ff10.parm7 < rna_ff10_min.rst7 > rna.pdb 2>/dev/null

# Now read in the PDB files, create all prmtops

cat > tleap.in << EOF
source leaprc.ff10

l = loadPDB protein.pdb
k = loadPDB dna.pdb
j = loadPDB rna.pdb

set default pbradii mbondi2

saveamberparm l protein_ff10.parm7 protein_ff10.rst7
saveamberparm k dna_ff10.parm7 dna_ff10.rst7
saveamberparm j rna_ff10.parm7 rna_ff10.rst7

quit
EOF

tleap -f tleap.in > /dev/null 2>&1

cat > tleap.in << EOF
source leaprc.ff99SB

l = loadPDB protein.pdb

set default pbradii mbondi2

saveamberparm l protein_ff99SB.parm7 protein_ff99SB.rst7
quit
EOF

tleap -f tleap.in > /dev/null 2>&1

cat > tleap.in << EOF
source leaprc.ff99bsc0

l = sequence {DA5 DC DT DG DA DC3}

set default pbradii mbondi2

saveamberparm l dna_ff99bsc0.parm7 dna_ff99bsc0.rst7

quit
EOF

tleap -f tleap.in > /dev/null 2>&1

cat > tleap.in << EOF
source leaprc.ff99bsc0
loadamberprep all_nuc94_ol_bsc0.in
loadamberparams frcmod.ol.dat

l = sequence {RA5 RC RU RG RA RC3}

set default pbradii mbondi2

saveamberparm l rna_ff99bsc0.parm7 rna_ff99bsc0.rst7

quit
EOF

tleap -f tleap.in > /dev/null 2>&1

cp dna_ff10_min.rst7 dna_ff99bsc0.rst7
cp rna_ff10_min.rst7 rna_ff99bsc0.rst7

# Now run sander for a short period so we can compare
cat > mdin << EOF
Implicit solvent molecular dynamics
 &cntrl
   imin=0, irest=0, ntx=1,
   ntpr=1, ntwx=2, nstlim=50,
   dt=0.002, ntt=3, tempi=300,
   temp0=300, gamma_ln=1.0, ig=1243,
   ntp=0, ntc=2, ntf=2, cut=100,
   ntb=0, igb=5,
 /
EOF

echo "Running Protein, ff10"
sander -O -o protein_ff10.mdout -p protein_ff10.parm7 -c protein_ff10.rst7 \
          -x protein_ff10.mdcrd

echo "Running Protein, ff99SB"
sander -O -o protein_ff99SB.mdout -p protein_ff99SB.parm7 -c protein_ff10.rst7 \
          -x protein_ff99SB.mdcrd

echo "Running DNA, ff10"
sander -O -o dna_ff10.mdout -p dna_ff10.parm7 -c dna_ff10.rst7 \
          -x dna_ff10.mdcrd

echo "Running DNA, ff99bsc0"
sander -O -o dna_ff99bsc0.mdout -p dna_ff99bsc0.parm7 -c dna_ff10.rst7 \
          -x dna_ff99bsc0.mdcrd

echo "Running RNA, ff10"
sander -O -o rna_ff10.mdout -p rna_ff10.parm7 -c rna_ff10.rst7 \
          -x rna_ff10.mdcrd

echo "Running RNA, ff99bsc0"
sander -O -o rna_ff99bsc0.mdout -p rna_ff99bsc0.parm7 -c rna_ff10.rst7 \
          -x rna_ff99bsc0.mdcrd

echo "Done with sander!"

echo "Extracting frcmod files to compare."

echo "Extracting protein frcmods"
ante-MMPBSA.py -p protein_ff10.parm7 -c protein_ff10.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod protein_ff10.frcmod
mv _AnteMMPBSA_.off protein_ff10.off
ante-MMPBSA.py -p protein_ff99SB.parm7 -c protein_ff99SB.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod protein_ff99SB.frcmod
mv _AnteMMPBSA_.off protein_ff99SB.off

echo "Extracting dna frcmods"
ante-MMPBSA.py -p dna_ff10.parm7 -c dna_ff10.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod dna_ff10.frcmod
mv _AnteMMPBSA_.off dna_ff10.off
ante-MMPBSA.py -p dna_ff99bsc0.parm7 -c dna_ff99bsc0.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod dna_ff99bsc0.frcmod
mv _AnteMMPBSA_.off dna_ff99bsc0.off

echo "Extracting rna frcmods"
ante-MMPBSA.py -p rna_ff10.parm7 -c rna_ff10.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod rna_ff10.frcmod
mv _AnteMMPBSA_.off rna_ff10.off
ante-MMPBSA.py -p rna_ff99bsc0.parm7 -c rna_ff99bsc0.rst7 2>/dev/null
mv _AnteMMPBSA_.frcmod rna_ff99bsc0.frcmod
mv _AnteMMPBSA_.off rna_ff99bsc0.off

echo "Done with frcmods!"

cat << EOF
To compare the results, compare the following files

protein_ff10.mdout protein_ff99SB.mdout
protein_ff10.frcmod protein_ff99SB.frcmod

dna_ff10.mdout dna_ff99bsc0.mdout
dna_ff10.frcmod dna_ff99bsc0.frcmod

rna_ff10.mdout rna_ff99bsc0.mdout
rna_ff10.frcmod rna_ff99bsc0.frcmod
EOF
