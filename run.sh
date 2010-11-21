#!/bin/bash

cat > mdin << EOF
TI calculation
&cntrl
   nstlim =1000000, nscm=2000,
   ntx=1, irest=0, ntpr=1000,
   tempi=0.0, temp0=300.0, ntt=3, gamma_ln=5.0,
   ntb=0, igb=5, cut=999.0,
   dt=0.001, 
   ntc=2, ntf=2, saltcon=0.1,
   ntwr = 10000, ntwx=1000, 
   icfe=1, clambda=0.0,
/
EOF

cat > groupfile << EOF
-O -i mdin -p hip.prmtop -c hip.inpcrd -o 0p.mdout -x 0p.mdcrd -r 0p.restrt
-O -i mdin -p hid.prmtop -c hid.inpcrd -o 0d.mdout -x 0d.mdcrd -r 0d.restrt
EOF

mpirun -np 2 sander.MPI -ng 2 -groupfile groupfile

for x in `seq 1 1 10`; do

  lambda=`echo "$x * 0.1" | bc`

cat > mdin << EOF
TI calculation
&cntrl
   nstlim =1000000, nscm=2000,
   ntx=5, irest=1, ntpr=1000,
   tempi=300.0, temp0=300.0, ntt=3, gamma_ln=5.0,
   ntb=0, igb=5, cut=999.0,
   dt=0.001, 
   ntc=2, ntf=2, saltcon=0.1,
   ntwr = 10000, ntwx=1000,
   icfe=1, clambda=$lambda,
/
EOF

  prev=`echo "$x - 1" | bc`

cat > groupfile << EOF
-O -i mdin -p hip.prmtop -c ${prev}p.restrt -o ${x}p.mdout -r ${x}p.restrt -x ${x}p.mdcrd
-O -i mdin -p hid.prmtop -c ${prev}p.restrt -o ${x}d.mdout -r ${x}d.restrt -x ${x}d.mdcrd
EOF

  mpirun -np 2 sander.MPI -ng 2 -groupfile groupfile

done

