#PBS -S /bin/bash
#PBS -l walltime=20:00:00
#PBS -l nodes=1:ppn=6
#PBS -o MasterTest.out
#PBS -e MasterTest.err

# Error message
error() {
   echo "Error: $1"
   exit 1
}

# Load the Amber git module
module load amber/homegit
module load intel/12.1.0
module load cuda/5.0
module load mpich2-intel/1.4.1p1_12.1.0

d=`date +%m-%d-%y`
cd $AMBERHOME
git checkout master 2>&1 > /dev/null || error "Checking out master"
git clean -f -x -d 2>&1 > /dev/null
git pull origin master  2>&1 > /dev/null || error "Pulling from origin"

echo "Configuring serial at `date`"
(./configure intel 2>&1) > $logdir/${d}_intel_serial_config.log
echo "Making serial at `date`"
(make -j6 install 2>&1) > $logdir/${d}_intel_serial_make.log || error "Making serial"

echo "Configuring parallel at `date`"
(./configure -mpi intel 2>&1) > $logdir/${d}_intel_parallel_config.log || \
                              error "Configuring parallel"
echo "Making parallel at `date`"
(make -j6 install 2>&1) > $logdir/${d}_intel_parallel_make.log || error "Making parallel"

# Now it's time to test

/bin/rm -fr logs/

echo "Beginning serial tests at `date`"

(make test.serial 2>&1) > $logdir/${d}_intel_serial_test.log

export DO_PARALLEL='mpirun -np 2'

echo "Beginning parallel tests with 2 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_intel_parallel_2_test.log

export DO_PARALLEL='mpirun -np 4'

echo "Beginning parallel tests with 4 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_intel_parallel_4_test.log

export DO_PARALLEL='mpirun -np 8'

echo "Beginning parallel tests with 8 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_intel_parallel_8_test.log

# Back up logs directory
mv logs $logdir/${d}_intel_logs

# Clean out my repo and go with the GNU compilers
echo ""
echo "Switching to the GNU compilers..."
git clean -f -x -d 2>&1 > /dev/null
module switch intel gcc
module switch mpich2-intel mpich2-gnu
echo ""

echo "Configuring serial at `date`"
(./configure gnu 2>&1) > $logdir/${d}_gnu_serial_config.log
echo "Making serial at `date`"
(make -j6 install 2>&1) > $logdir/${d}_gnu_serial_make.log || error "Making serial"

echo "Configuring parallel at `date`"
(./configure -mpi gnu 2>&1) > $logdir/${d}_gnu_parallel_config.log || \
                              error "Configuring parallel"
echo "Making parallel at `date`"
(make -j6 install 2>&1) > $logdir/${d}_gnu_parallel_make.log || error "Making parallel"

# Now it's time to test

echo "Beginning serial tests at `date`"

(make test.serial 2>&1) > $logdir/${d}_gnu_serial_test.log

export DO_PARALLEL='mpirun -np 2'

echo "Beginning parallel tests with 2 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_gnu_parallel_2_test.log

export DO_PARALLEL='mpirun -np 4'

echo "Beginning parallel tests with 4 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_gnu_parallel_4_test.log

export DO_PARALLEL='mpirun -np 8'

echo "Beginning parallel tests with 8 cores at `date`"

(make test.parallel 2>&1) > $logdir/${d}_gnu_parallel_8_test.log

mv logs $logdir/${d}_gnu_logs

echo "Cleaning everything up..."
git clean -f -x -d 2>&1 > /dev/null

echo "Done at `date`"
