#!/bin/bash

######################################################################
#      Script to scan a dihedral angle from 0 to 360 degrees         #
#          Jason Swails                3/16/2009                     #
######################################################################

# User input Variables

NUM_STEPS=72         # Number of points in PES scan (increment is 360/NUM_STEPS degrees)

prmtop="prmtop"      # The name of the topology file for the system
inpcrd="inpcrd"      # Some initial coordinate file (it will be replaced during the scan)

# Four atoms that define the dihedral that you wish to scan
atom1=1
atom2=11
atom3=12
atom4=13

##########################################################################################
#              Begin script.  SANDER serial version must be in your path                 #
##########################################################################################

mkdir -p Restrts
mkdir -p Jarlogs
mkdir -p Mdinfos
mkdir -p Mdins
mkdir -p Mdouts
mkdir -p RSTFiles

stepNum=1
increment=`echo "360 / $NUM_STEPS" | bc`

for x in `seq 0 1 $NUM_STEPS`; do

dihedral=`echo "$x * $increment" | bc`

left=`echo "$dihedral - 100" | bc`
right=`echo "$dihedral + 100" | bc`

cat > dihedral.$x.RST << EOF
torsion restraint for polymer number $x in scan
&rst iat=${atom1},${atom2},${atom3},${atom4} r1=${left}, r2=${dihedral}, r3=${dihedral}, r4=${right}, rk2=5000., rk3=5000., /
EOF


cat > mdin.$x.scan << EOF
Basic minimization with a restraint

&cntrl
        imin=1,
        cut=1000,
        maxcyc=10000,
        ncyc=50,
        drms=1E-6,
        nmropt=1,
        igb=6,
        ntb=0,
        cut=1000,
        ntc=1,
        ntf=1,
/
&wt type='DUMPFREQ', istep1=10 /
&wt type='END'   /
DISANG=dihedral.${x}.RST
LISTIN=POUT
LISTOUT=POUT
DUMPAVE=jar.${x}.log
EOF

sander -O -i mdin.$x.scan -o mdout.$x.scan -p $prmtop -c $inpcrd -r restrt

cp restrt $inpcrd
mv restrt Restrts/restrt.$x
mv mdinfo Mdinfos/mdinfo.$x
mv mdout.$x.scan Mdouts/
mv mdin.$x.scan Mdins/
mv jar.$x.log Jarlogs/
mv dihedral.$x.RST RSTFiles/

done


# Post-process the data

# Rename files jar.0.log - jar.9.log to jar.00.log - jar.09.log so everything
# is done in order (same for mdinfo.0-mdinfo.9)

for x in `seq 0 1 9`; do
cd Mdinfos
mv mdinfo.$x mdinfo.0$x
cd ../Jarlogs
mv jar.$x.log jar.0$x.log
cd ../
done

# Parse the output files and create a file called "profile.dat" which is the desired data file

cd Jarlogs/
tail -n 1 jar* | awk '$1=="==>" {getline;print $2}' > angles.dat
cd ../Mdinfos
grep EAMBER mdinfo* | awk '{print $4}' > energies.dat
cd ../
mv Jarlogs/angles.dat Mdinfos/energies.dat .
paste -d" " angles.dat energies.dat > profile.dat

