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

# Set the main git repo as Amber installation to use, check out master, update
# it, clean all files, then build, test, etc.

mkdir -p /home/swails/AmberBuildTest
logdir=/home/swails/AmberBuildTest
oldamber=`cat ~/.last_amber_choice`
_set_amber_var /home/swails/amber
printf $oldamber > ~/.last_amber_choice

d=`date +%m-%d-%y`
cd $AMBERHOME
git checkout master || error "Checking out master"
git clean -f -x -d 2>&1 > /dev/null
git pull origin master || error "Pulling from origin"

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
source /usr/local/var/mpi-selector/data/mpich2-1.4.1p1-gnu-4.5.3.sh
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
