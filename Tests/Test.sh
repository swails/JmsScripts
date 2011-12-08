#!/bin/sh

# This script tests all of the scripts in ScriptDevel that actually have tests
# associated with them here.

# mdout.py
echo "============================================================"
echo "Testing mdout.py"
printf "   testing minimization output:   "

python ../mdout.py --mdout trpcage.nowat.min.mdout > tmp
diff mdout_trpcage_min.check tmp > mdout_trpcage_min.diff 2>&1

if [ $? -ne 0 ]; then
   echo "Failed. See mdout_trpcage_min.diff"
else
   echo "PASSED"
   /bin/rm -f mdout_trpcage_min.diff
fi

printf "   testing heating output:        "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout > tmp
diff mdout_trpcage_heat.check tmp > mdout_trpcage_heat.diff 2>&1

if [ $? -ne 0 ]; then
   echo "Failed. See mdout_trpcage_heat.diff"
else
   echo "PASSED"
   /bin/rm -f mdout_trpcage_heat.diff
fi

printf "   testing combining mdout files: "

python ../mdout.py --mdout trpcage.nowat.heat0.mdout \
                   --next-mdout trpcage.nowat.heat1.mdout > tmp
diff mdout_trpcage_combine.check tmp > mdout_trpcage_combine.diff 2>&1

if [ $? -ne 0 ]; then
   echo "Failed. See mdout_trpcage_combine.diff"
else
   echo "PASSED"
   /bin/rm -f mdout_trpcage_combine.diff
fi

echo "============================================================"

/bin/rm -f tmp
