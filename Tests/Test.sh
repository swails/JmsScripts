#!/bin/sh

test_cleanup() {
if [ $1 -ne 0 ]; then
   echo "Failed. See $2"
else
   echo "PASSED"
   /bin/rm -f $2
fi
}
   
# This script tests all of the scripts in ScriptDevel that actually have tests
# associated with them here.

# mdout.py
echo "============================================================"
echo "Testing mdout.py"
printf "   testing minimization output:   "

python ../mdout.py --mdout trpcage.nowat.min.mdout > tmp
diff mdout_trpcage_min.check tmp > mdout_trpcage_min.diff 2>&1
test_cleanup $? mdout_trpcage_min.diff

printf "   testing heating output:        "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout > tmp
diff mdout_trpcage_heat.check tmp > mdout_trpcage_heat.diff 2>&1
test_cleanup $? mdout_trpcage_heat.diff

printf "   testing combining mdout files: "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout \
                   --next-mdout trpcage.nowat.heat1.mdout > tmp
diff mdout_trpcage_combine.check tmp > mdout_trpcage_combine.diff 2>&1
test_cleanup $? mdout_trpcage_combine.diff

# phipsigen.py
echo "============================================================"
echo "Testing phipsigen.py"

../phipsigen.py -i phipsi.in -o rama \
                -p acfca.neutral.state00.parm7 \
                -y acfca.neutral.nowat.md.nc > /dev/null

printf "   Checking residue 2 data:    "
diff rama.2.dat.check rama.2.dat > rama.2.dat.diff
test_cleanup $? rama.2.dat.diff

printf "   Checking residue 3 data:    "
diff rama.3.dat.check rama.3.dat > rama.3.dat.diff
test_cleanup $? rama.3.dat.diff

printf "   Checking residue 2 binning: "
diff rama.2.bin.dat.check rama.2.bin.dat > rama.2.bin.diff
test_cleanup $? rama.2.bin.diff

printf "   Checking residue 3 binning: "
diff rama.3.bin.dat.check rama.3.bin.dat > rama.3.bin.diff
test_cleanup $? rama.3.bin.diff

/bin/rm -f _PHIPSI_* rama.?.gnu rama.?.dat rama.?.bin.dat

# 1Dbinning.py
echo "============================================================"
echo "Testing 1Dbinning.py"
../1Dbinning.py -f rama.2.dat.check -o rama_1dbinned.dat \
                -b 30 -n -c 1 > /dev/null

printf "   Checking binning of 1D data set: "
diff rama_1dbinned.dat.check rama_1dbinned.dat > rama_1dbinned.dat.diff
test_cleanup $? rama_1dbinned.dat.diff
/bin/rm -f rama_1dbinned.dat

# 2Dbinning.py
echo "============================================================"
echo "Testing 2Dbinning.py"
../2Dbinning.py -f rama.2.dat.check -o rama_2dbinned.dat \
                -b 50x50 -n -c 1,2 -x "-176.0-180" > /dev/null

printf "   Checking binning of 2D data set: "
diff rama_2dbinned.dat.check rama_2dbinned.dat > rama_2dbinned.dat.diff
test_cleanup $? rama_2dbinned.dat.diff
/bin/rm -f rama_2dbinned.dat

# mdcrd.py
echo "============================================================"
echo "Testing mdcrd.py"
python ../mdcrd.py -p trpcage.nowat.parm7 trpcage.solv5.[1-5]_remd12.nc > tmp

printf "   Checking output of mdcrd.py: "
diff mdcrd.out.check tmp > mdcrd.out.diff
test_cleanup $? mdcrd.out.diff

printf "   Checking the RMSD data set:  "
diff AmberTraj_RMSD.dat.check AmberTraj_RMSD.dat > AmberTraj_RMSD.dat.diff
test_cleanup $? AmberTraj_RMSD.dat.diff

printf "   Checking mdcrd.py log file:  "
diff mdcrd_py.log.check mdcrd_py.log > mdcrd_py.log.diff
test_cleanup $? mdcrd_py.log.diff
/bin/rm -f AmberTraj_RMSD.dat mdcrd_py.log
echo "============================================================"

/bin/rm -f tmp
