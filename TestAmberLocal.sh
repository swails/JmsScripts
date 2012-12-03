#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=6
#PBS -o LocalTest.out
#PBS -e LocalTest.err

unset DO_PARALLEL TESTsander

# First, update my code
cd $AMBERHOME
git pull 2>&1 > /dev/null
git clean -f -x -d 2>&1 > /dev/null

echo "Configuring serial at `date`"
(./configure intel 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_serial_config.log
echo "Making serial at `date`"
(make -j6 install 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_serial_make.log

echo "Configuring parallel at `date`"
(./configure -mpi intel 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_parallel_config.log || \
                              error "Configuring parallel"
echo "Making parallel at `date`"
(make -j6 install 2>&1) > $HOME/TestLogs/`date +%m+%d+%y`_parallel_make.log

# Now test.
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
