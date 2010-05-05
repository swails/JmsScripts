#!/bin/bash

# The purpose of this script is to process raw jar.log files generated from
# SMD simulations with AMBER.  To use this script, place it in the directory
# with each (and ONLY) of the separate jar.log files you wish to process

#Enter the name of your desired file here:
clear
echo "Enter the name of the output file."
read OUTFILE

counter=0
for x in `ls -1`; do

if [ $x != "SMDprep.sh" ]; then

if [ $counter == 0 ]; then

awk '{print $1}' $x > holder
counter=`echo "$counter + 1" | bc`
awk '{print $4}' $x > temp
paste -d" " holder temp > temp2
cat temp2 > holder

else

awk '{print $4}' $x > temp

paste -d" " holder temp > temp2

cat temp2 > holder

counter=`echo "$counter + 1" | bc`

fi 

fi

done

mv holder $OUTFILE
rm temp temp2

echo "You processed $counter SMD trajectories.  The generated file is $OUTFILE"
