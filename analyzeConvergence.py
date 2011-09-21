#!/usr/bin/env python

import sys, math, os, utilities
from optparse import OptionParser

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

clparser = OptionParser()
clparser.add_option('-o','--mdout',dest='mdout',help='MDOUT file to analyze')
clparser.add_option('-p','--prmtop',dest='prmtop',help='PRMTOP of the system')
clparser.add_option('-y','--mdcrd',dest='mdcrd',help='MDCRD of the process')
clparser.add_option('-g','--gnu-prefix',dest='prefix',
                    default='gnu', help='Prefix of gnuplot scripts')
opt, args = clparser.parse_args()

if not opt.mdout or not opt.prmtop or not opt.mdcrd:
   clparser.print_help()
   sys.exit()

try:
   mdout = open(opt.mdout,'r')
   if not os.path.exists(opt.mdcrd): raise IOError('mdcrd doesn\'t exist')
   if not os.path.exists(opt.prmtop): raise IOError('mdcrd doesn\'t exist')
except IOError:
   print >> sys.stderr, "Error: mdout, prmtop, and/or mdcrd files don't exist!"
   printusage()

ptraj = utilities.which('ptraj')

if ptraj == "none":
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
      try:
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
         y += 3
         words = mdoutlines[y].split()
         restr.append(float(words[8]))
         eamber.append(eptot[len(eptot)-1]-restr[len(restr)-1])
         if constp:
            y += 2
            words = mdoutlines[y].split()
            density.append(float(words[2]))
      except IndexError:
         print >> sys.stderr, mdoutlines[y]
         print >> sys.stderr, 'Error occurred on parsing the last line.'
         sys.exit()

if opt.mdcrd != '':
   print >> sys.stdout, "Gathering RMS data with %s..." % ptraj
   trajin = open('__ptraj.in','w')
   trajin.write("""trajin %s
strip :WAT:Na+:Cl-:Br-:Cs+:F-:I-:K+:Li+:Mg+:MG2:Rb+:IB:CIO
rms first mass out __RMS.dat
""" % opt.mdcrd)
   trajin.close()
   os.system("%s %s __ptraj.in 2>>__ptraj.out 1>>__ptraj.out" % (ptraj,opt.prmtop))

print >> sys.stdout, "Printing data files..."

feamber = open('__EAMBER.dat','w')
for x in range(len(eamber)):
   feamber.write("%s %s\n" % (x,eamber[x]))
feamber.close()

ftemp = open('__TEMP.dat','w')
for x in range(len(temp)):
   ftemp.write("%s %s\n" % (x,temp[x]))
ftemp.close()

if constp:
   fpres = open('__PRESSURE.dat','w')
   for x in range(len(pres)):
      fpres.write('%s %s\n' % (x,pres[x]))
   fpres.close()

fektot = open('__EKTOT.dat','w')
for x in range(len(ektot)):
   fektot.write("%s %s\n" % (x,ektot[x]))
fektot.close()

feptot = open('__EPTOT.dat','w')
for x in range(len(eptot)):
   feptot.write("%s %s\n" % (x,eptot[x]))
feptot.close()

fetot = open('__ETOT.dat','w')
for x in range(len(etot)):
   fetot.write("%s %s\n" % (x,etot[x]))
fetot.close()

if constp:
   fdens = open('__DENSITY.dat','w')
   for x in range(len(density)):
      fdens.write("%s %s\n" % (x,density[x]))
   fdens.close()

print >> sys.stdout, "Printing gnuplot scripts..."

gnu = open("%s.energy" % (opt.prefix),'w')
gnu.write("""set font 'calibri,20'
set title 'Energy equilibration'
plot '__ETOT.dat' w l lw 2 ti 'ETOT','__EKTOT.dat' w l lw 2 ti 'EKTOT', \
     '__EPTOT.dat' w l lw 2 ti 'EPTOT','__EAMBER.dat' w l lw 2 lt -1 ti 'EAMBER'
pause -1
""")
gnu.close()

gnu = open('%s.temp' % (opt.prefix),'w')
gnu.write("""set font 'calibri,20'
set title 'Temperature equilibration'
plot '__TEMP.dat' w l lw 2 lt -1 ti 'TEMP'
pause -1
""")
gnu.close()

if constp:
   gnu = open('%s.pressure' % (opt.prefix),'w')
   gnu.write("""set font 'calibri,20'
set title 'Pressure equilibration'
plot '__PRESSURE.dat' w l lw 2 lt -1 ti 'Pressure'
pause -1
""")
   gnu.close()
   gnu = open('%s.density' % (opt.prefix),'w')
   gnu.write("""set font 'calibri,20'
set title 'Density equilibration'
plot '__DENSITY.dat' w l lw 2 lt -1 ti 'Density'
pause -1
""")
   gnu.close()

gnu = open('%s.rms' % (opt.prefix),'w')
gnu.write("""set font 'calibri,20'
set title 'RMS equilibration'
plot '__RMS.dat' w l lt -1 lw 2 ti 'RMS'
pause -1
""")
gnu.close()

