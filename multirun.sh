#!/bin/sh

#######################################################################
#                                                                     #
# This script will create numerous calculation setups from a single   #
# set of files such that the calculation can be trivially parallel-   #
# ized across an arbitrary number of processors on a compute cluster  #
# that operates via a queue or submission system. (i.e. PBS). Note    #
# that users will need some familiarity with shell scripts to edit    #
# this script to work for their particular system. The user will need #
# to edit the INPUT FILE and SUBMISSION SCRIPT. Searching for these   #
# terms will bring you to where they are printed, and they should be  #
# edited there. The only other variables needed are the number of     #
# chunks the user wishes to break their calculation up into, the      #
# overall starting and ending frame they wish to perform analysis on, #
# and the "top-level" prmtop and mdcrd (solvated if initial_traj != 1 #
# or complex if initial_traj=1). You may also have to edit SUBMISSION #
# COMMAND for your scheduler.                                         #
#                                                                     #
#######################################################################


# BE CAREFUL TO MAKE SURE THAT THE NUMBER OF FRAMES YOU WANT TO DO MMPBSA
# FOR IS AN INTERVAL OF THE NUMBER OF CHUNKS YOU WANT TO BREAK IT UP INTO
# OR YOU MAY GET THE LAST CHUNK CALCULATING MOST OF THE FRAMES!!

# User-defined variables
numchunks=10
globalstartframe=1
globalendframe=100
interval=1
prmtop="solvated_prmtop"
mdcrd="mdcrd"

# Don't forget to edit INPUT FILE and SUBMISSION SCRIPT

#############################  BEGIN SCRIPT  ##########################

# Make sure globalendframe is not set higher than the number of frames
# in the mdcrd -- if it is, set it to the total number of frames

numatm=`awk '$2=="POINTERS"{getline; getline; print $1}' $prmtop`
numlines_per_frame=`echo "$numatm * 3 / 10" | bc`
chk=`echo "$numlines_per_frame * 10" | bc`

if [ $numatm -gt $chk ]; then
   numlines_per_frame=`echo "$numlines_per_frame + 1" | bc`
fi

numlines=`wc -l $mdcrd`
# strip the title line
numlines=`echo "$numlines - 1" | bc`
numframes=`echo "$numlines / $numlines_per_frame" | bc`
chk=`echo "$numlines_per_frame * $numframes" | bc`

if [ $numlines -ne $chk ]; then
   echo "Error: Number of frames in mdcrd is not what's expected!"
   echo "       $prmtop and $mdcrd are likely incompatible."
   exit 1
fi

if [ -d mmpbsa_chunk1 ]; then
   echo "Error: mmpbsa_chunk1 directory already exists! Remove chunk"
   echo "       directories before you rerun this script!"
fi

if [ $globalendframe -gt $numframes ]; then
   globalendframe=$numframes
fi

totframes=`echo "$globalendframe - $globalstartframe" | bc`
totframes=`echo "$totframes / $interval" | bc`

if [ $totframes -lt $chunks ]; then
   echo "Error: You're asking to break the calculation into more chunks"
   echo "       than you have frames!"
   exit 1
fi

framesperchunk=`echo "$totframes / $chunks" | bc`

startframe=$globalstartframe
endframe=`echo "$globalstartframe + $framesperchunk" | bc`

# loop through all chunks
for x in `seq 1 1 $chunks`; do

   if [ $x -eq $chunks ]; then
      endframe=$globalendframe
   fi

   mkdir mmpbsa_chunk${x}
   cd mmpbsa_chunk${x}

# INPUT FILE
   cat > mmpbsa.in << EOF
MMPBSA input file for chunk $x 
&general
   startframe=$startframe, endframe=$endframe, interval=$interval,
   verbose=1, initial_traj=1,
/
&gb
   saltcon=0.1,
/
EOF

# SUBMISSION SCRIPT
   cat > pbs.jobfile << EOF
#!/bin/bash

#PBS -N mmpbsa_chunk$x
#PBS -o pbs.out
#PBS -e pbs.err
#PBS -m abe
#PBS -M email_addy@place.com
#PBS -q queue_name
#PBS -l nodes=1:ppn=1
#PBS -l pmem=600mb

cd \$PBS_O_WORKDIR

MMPBSA.py -i mmpbsa.in -sp ../$prmtop -cp ../complex_prmtop -rp ../receptor_prmtop \\
          -lp ../ligand_prmtop -y ../$mdcrd
EOF

# SUBMISSION COMMAND
   qsub pbs.jobfile

   cd ../

done

