#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=6
#PBS -o LocalTest.out
#PBS -e LocalTest.err

cd $AMBERHOME

unset DO_PARALLEL TESTsander

/bin/rm -fr logs/

echo "Beginning serial tests at `date`"

(make test.serial 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_serial_test.log

export DO_PARALLEL='mpirun -np 2'

echo "Beginning parallel tests with 2 cores at `date`"

(make test.parallel 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_parallel_2_test.log

export DO_PARALLEL='mpirun -np 4'

echo "Beginning parallel tests with 4 cores at `date`"

(make test.parallel 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_parallel_4_test.log

export DO_PARALLEL='mpirun -np 8'

echo "Beginning parallel tests with 8 cores at `date`"

(make test.parallel 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_parallel_8_test.log
