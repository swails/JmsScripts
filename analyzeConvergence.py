#!/usr/bin/env python

import sys, math, os, utilities

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

def printusage():
   import sys
   print "analyzeConvergence.py <mdout> <mdcrd> <prmtop> {<gnuplot_script_prefix>}"
   sys.exit()

if len(sys.argv) < 2:
   printusage()
if sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--help':
   printusage()

mdcrd = '' 
prefix = 'gnu'

try:
   mdout = open(sys.argv[1],'r')
   if len(sys.argv) > 3:
      mdcrd = sys.argv[2]
      prmtop = sys.argv[3]
      tmp = open(mdcrd,'r')
      tmp2 = open(prmtop,'r')
   if len(sys.argv) > 4:
      prefix = sys.argv[4]
except IOError:
   print >> sys.stderr, "Error: mdout, prmtop, and/or mdcrd files don't exist!"
   printusage()

tmp.close()
tmp2.close()

ptraj = utilities.which('ptraj')

if ptraj == "none" and mdcrd != '':
   print "Error: ptraj is needed to analyze coordinate files!"
   printusage()

mdoutlines = mdout.readlines()
mdout.close()

# parse this block:
# NSTEP =   100000   TIME(PS) =     200.000  TEMP(K) =   308.77  PRESS =     0.0
# Etot   =    -18650.7357  EKtot   =      4532.3000  EPtot      =    -23183.0357
# BOND   =        11.0060  ANGLE   =        19.7379  DIHED      =        12.5913
# 1-4 NB =         5.9918  1-4 EEL =      -372.1093  VDWAALS    =      3276.9971
# EELEC  =    -26139.9559  EHBOND  =         0.0000  RESTRAINT  =         2.7054
# EAMBER (non-restraint)  =    -23185.7411
# Ewald error estimate:   0.1482E-03

step = []
temp = []
pres = []
etot = []
ektot = []
eptot = []
restr = []
eamber = []
density = []

print >> sys.stdout, "Parsing mdout file..."

for x in range(len(mdoutlines)):

   if "A V E R A G E S" in mdoutlines[x]: # we've reached the end
      break
      
   # determine constant pressure or not
   if mdoutlines[x].startswith("Potential function:"):
      y = x + 1
      words = mdoutlines[y].split()
      if words[5] == "2,":
         constp = True
      else:
         constp = False

   # get potential terms
   if mdoutlines[x].startswith(" NSTEP"):
      y = x
      words = mdoutlines[y].split()
      step.append(int(words[2]))
      temp.append(float(words[8]))
      pres.append(float(words[11]))
      y += 1
      words = mdoutlines[y].split()
      etot.append(float(words[2]))
      ektot.append(float(words[5]))
      eptot.append(float(words[8]))
      y += 4
      words = mdoutlines[y].split()
      restr.append(float(words[8]))
      eamber.append(eptot[len(eptot)-1]-restr[len(restr)-1])
      if constp:
         y += 2
         words = mdoutlines[y].split()
         density.append(float(words[2]))

if mdcrd != '':
   print >> sys.stdout, "Gathering RMS data with {0}...".format(ptraj)
   trajin = open('__ptraj.in','w')
   trajin.write("""trajin {0}
strip :WAT
rms first mass out __rms.dat
""".format(mdcrd))
   trajin.close()
   os.system("{0} {1} __ptraj.in 2>>__ptraj.out 1>>__ptraj.out")

print >> sys.stdout, "Printing data files..."

feamber = open('__EAMBER.dat','w')
for x in range(len(eamber)):
   feamber.write("{0} {1}\n".format(step[x],eamber[x]))
feamber.close()

ftemp = open('__TEMP.dat','w')
for x in range(len(temp)):
   ftemp.write("{0} {1}\n".format(step[x],temp[x]))
ftemp.close()

if constp:
   fpres = open('__PRESSURE.dat','w')
   for x in range(len(pres)):
      fpres.write('{0} {1}\n'.format(step[x],pres[x]))
   fpres.close()

fektot = open('__EKTOT.dat','w')
for x in range(len(ektot)):
   fektot.write("{0} {1}\n".format(step[x],ektot[x]))
fektot.close()

feptot = open('__EPTOT.dat','w')
for x in range(len(eptot)):
   feptot.write("{0} {1}\n".format(step[x],eptot[x]))
feptot.close()

fetot = open('__ETOT.dat','w')
for x in range(len(etot)):
   fetot.write("{0} {1}\n".format(step[x],etot[x]))
fetot.close()

if constp:
   fdens = open('__DENSITY.dat','w')
   for x in range(len(density)):
      fdens.write("{0} {1}\n".format(step[x],density[x]))
   fdens.close()

print >> sys.stdout, "Printing gnuplot scripts..."

gnu = open("{0}.energy".format(prefix),'w')
gnu.write("""set font 'calibri,20'
set title 'Energy equilibration'
plot '__ETOT.dat' w l lw 2 ti 'ETOT','__EKTOT.dat' w l lw 2 ti 'EKTOT', \
     '__EPTOT.dat' w l lw 2 ti 'EPTOT','__EAMBER.dat' w l lw 2 lt -1 ti 'EAMBER'
""")
gnu.close()

gnu = open('{0}.temp'.format(prefix),'w')
gnu.write("""set font 'calibri,20'
set title 'Temperature equilibration'
plot '__TEMP.dat' w l lw 2 lt -1 ti 'TEMP'
""")
gnu.close()

if constp:
   gnu = open('{0}.pressure'.format(prefix),'w')
   gnu.write("""set font 'calibri,20'
set title 'Pressure equilibration'
plot '__PRESSURE.dat' w l lw 2 lt -1 ti 'Pressure'
""")
   gnu.close()
   gnu = open('{0}.density'.format(prefix),'w')
   gnu.write("""set font 'calibri,20'
set title 'Density equilibration'
plot '__DENSITY.dat' w l lw 2 lt -1 ti 'Density'
""")
   gnu.close()

if mdcrd != '':
   gnu = open('{0}.rms'.format(prefix),'w')
   gnu.write("""set font 'calibri,20'
set title 'RMS equilibration'
plot '__RMS.dat' w l lt -1 lw 2 ti 'RMS'
""")
   gnu.close()

