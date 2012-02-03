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
diff -Nru mdout_trpcage_min.check tmp > mdout_trpcage_min.diff 2>&1
test_cleanup $? mdout_trpcage_min.diff

printf "   testing heating output:        "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout > tmp
diff -Nru mdout_trpcage_heat.check tmp > mdout_trpcage_heat.diff 2>&1
test_cleanup $? mdout_trpcage_heat.diff

printf "   testing combining mdout files: "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout \
                   --next-mdout trpcage.nowat.heat1.mdout > tmp
diff -Nru mdout_trpcage_combine.check tmp > mdout_trpcage_combine.diff 2>&1
test_cleanup $? mdout_trpcage_combine.diff

# phipsigen.py
echo "============================================================"
echo "Testing phipsigen.py"

../phipsigen.py -i phipsi.in -o rama \
                -p acfca.neutral.state00.parm7 \
                -y acfca.neutral.nowat.md.nc > /dev/null

printf "   Checking residue 2 data:    "
diff -Nru rama.2.dat.check rama.2.dat > rama.2.dat.diff
test_cleanup $? rama.2.dat.diff

printf "   Checking residue 3 data:    "
diff -Nru rama.3.dat.check rama.3.dat > rama.3.dat.diff
test_cleanup $? rama.3.dat.diff

printf "   Checking residue 2 binning: "
diff -Nru rama.2.bin.dat.check rama.2.bin.dat > rama.2.bin.diff
test_cleanup $? rama.2.bin.diff

printf "   Checking residue 3 binning: "
diff -Nru rama.3.bin.dat.check rama.3.bin.dat > rama.3.bin.diff
test_cleanup $? rama.3.bin.diff

/bin/rm -f _PHIPSI_* rama.?.gnu rama.?.dat rama.?.bin.dat

# 1Dbinning.py
echo "============================================================"
echo "Testing 1Dbinning.py"
../1Dbinning.py -f rama.2.dat.check -o rama_1dbinned.dat \
                -b 30 -n -c 1 > /dev/null

printf "   Checking binning of 1D data set: "
diff -Nru rama_1dbinned.dat.check rama_1dbinned.dat > rama_1dbinned.dat.diff
test_cleanup $? rama_1dbinned.dat.diff
/bin/rm -f rama_1dbinned.dat

# 2Dbinning.py
echo "============================================================"
echo "Testing 2Dbinning.py"
../2Dbinning.py -f rama.2.dat.check -o rama_2dbinned.dat \
                -b 50x50 -n -c 1,2 -x "-176.0-180" > /dev/null

printf "   Checking binning of 2D data set: "
diff -Nru rama_2dbinned.dat.check rama_2dbinned.dat > rama_2dbinned.dat.diff
test_cleanup $? rama_2dbinned.dat.diff
/bin/rm -f rama_2dbinned.dat

# mdcrd.py
echo "============================================================"
echo "Testing mdcrd.py"
python ../mdcrd.py -p trpcage.nowat.parm7 trpcage.solv5.[1-5]_remd12.nc > tmp

printf "   Checking output of mdcrd.py: "
diff -Nru mdcrd.out.check tmp > mdcrd.out.diff
test_cleanup $? mdcrd.out.diff

printf "   Checking the RMSD data set:  "
diff -Nru AmberTraj_RMSD.dat.check AmberTraj_RMSD.dat > AmberTraj_RMSD.dat.diff
test_cleanup $? AmberTraj_RMSD.dat.diff

printf "   Checking mdcrd.py log file:  "
diff -Nru mdcrd_py.log.check mdcrd_py.log > mdcrd_py.log.diff
test_cleanup $? mdcrd_py.log.diff
/bin/rm -f AmberTraj_RMSD.dat mdcrd_py.log
echo "============================================================"
echo "Testing remd.py"
python ../remd.py -l rem1.log -t TEMP -o tmp

printf "   Checking T-REM log analysis:  "
diff -Nru temp_remlog.stats.check tmp > temp_remlog.stats.diff
test_cleanup $? temp_remlog.stats.diff

python ../remd.py -l phrem.log -t pH -o tmp
printf "   Checking pH-REM log analysis: "
diff -Nru ph_remlog.stats.check tmp > ph_remlog.stats.diff
test_cleanup $? ph_remlog.stats.diff
echo "============================================================"

/bin/rm -f tmp
