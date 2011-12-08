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
echo "============================================================"

/bin/rm -f tmp
