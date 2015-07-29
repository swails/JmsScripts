#!/bin/bash

get_curl() {
   # get_curl will download a PDB file with the 4-letter code given as the only
   # argument using 'curl' (for when wget is unavailable)
   if [ $# -ne 1 ]; then
      echo "Bad usage: get_curl <PDB ID>"
      exit 1
   fi

   pdbcode=`python -c "print('$1'.lower())"`
   mid=`python -c "print(\"$pdbcode\"[1:3])"`

   url=ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/$mid/pdb$pdbcode.ent.gz
   curl $url > $1.pdb.gz
   retcode=$?
   if [ $retcode -eq 78 ]; then
      echo "Could not find PDB $1!"
      exit 1
   elif [ $retcode -ne 0 ]; then
      echo "Unknown error in downloading PDB $1!"
      exit 1
   else
      gunzip $1.pdb.gz
   fi
}

get_wget() {
   # get_wget will use wget to download a PDB file specified by the 4-letter
   # code
   if [ $# -ne 1 ]; then
      echo "Bad usage: wget <PDB ID>"
      exit 1
   fi

   pdbcode=`python -c "print('$1'.lower())"`
   mid=`python -c "print('$pdbcode'[1:3])"`

   url=ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/$mid/pdb$pdbcode.ent.gz
   wget $url
   retcode=$?
   if [ $retcode -eq 8 ]; then
      echo "Could not find PDB $1!"
      exit 1
   elif [ $retcode -ne 0 ]; then
      echo "Unknown error in downloading PDB $1!"
      exit 1
   else
      gunzip pdb$pdbcode.ent.gz && mv pdb$pdbcode.ent $1.pdb
   fi
}

usage() {
   echo "`basename $0` <PDB ID> [<PDB ID> [<PDB ID> [...]]]"
}

if [ $# -eq 0 ]; then
   usage
   exit 1
fi

wget_="F"
curl_="F"
# See if we have wget and/or curl
(which wget 2>&1) > /dev/null
test $? -eq 0 && wget_="T"
(which curl 2>&1) > /dev/null
test $? -eq 0 && curl_="T"

if [ "$wget_" = "F" -a "$curl_" = "F" ]; then
   echo "You must have either 'curl' or 'wget' to download PDBs!"
   exit 1
fi

# Now go through all arguments and download them all
while [ $# -gt 0 ]; do
   if [ "$wget_" = "T" ]; then
      get_wget $1
   else
      get_curl $1
   fi
   shift
done
